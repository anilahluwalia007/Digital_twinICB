"""
icb_bio.py — cell biology, seed train (N-3 → N) and production-reactor sizing
NavNevAI_mAbWasteWater · CarbonTesseract 2026

Port of CHO_mAb_Process_Calculator5.xlsx (same equations, same default values):
  Bio_inputs / Bio_sim / Bio_curve  bio_run(p)          N cell model, Euler steps of bio.dt_d
  Seed_sim                          seed_sim(b, X0, d)  one seed stage at 37 °C, same biology
  Seed_train                        seed_train(p)       N-3, N-2, N-1 and N back to the vial
  Demand                            demand(p)           N size from annual demand and the cell model

Culture mode per stage (culture.preset; preset "custom" reads culture.mode_n3 .. mode_n):
  batch      no feed, nothing leaves (the workbook's seed stages)
  fedbatch   concentrated feed from bio.feed_start_d, volume grows (the workbook's N)
  perfusion  constant volume: fresh medium in, cell-free permeate out through a retention device
             whose product passage falls as it fouls (swapped below a limit); on N a cell bleed
             holds culture.perf_X_setpoint; optional loss of productivity per generation or day
  presets    workbook (default)  N-3/N-2/N-1 batch, N fed-batch: identical to the workbook
             fb_all              every stage fed-batch
             fb_n1_perf_n        N-3/N-2 batch, N-1 fed-batch, N perfusion
             perf_all            every stage perfusion
The perfusion model is not in the workbook: its defaults are typical values, not site data.

Basis switches (workbook Seed_train B6 / B7), each "biology" or "rule":
  seedtrain.duration_basis  stage durations from Seed_sim (lag, nutrients, waste) or ln(fold) / mu_rule_d
  seedtrain.density_basis   transfer where net growth falls to growth_threshold_pct of its peak
                            (late exponential phase) or the fixed xt_rule_* densities
Rules of thumb used whatever the basis (as in the workbook README): seed densities, the growth
threshold, qO2, purification yield, batch success, turnaround and every check limit.

Twin hooks
  bio_titre_points(p)        upstream.molecule = "biology": N titre curve from the cell model
  apply_seed_train_to_n1(p)  seedtrain.apply_to_n1 = True: the twin's N-1 runs as the N-1 stage
                             of a seed train sized for the twin N reactor (bioreactor.V_L)
  bio_design(p)              everything for the Outputs and Design tabs (result["bio"])

Units: day, L, g/L, mM, mOsm/kg, 1e6 cells/mL, qP in pg/cell/day.
Check against the workbook:  python3 icb_bio.py
"""
from __future__ import annotations
import math
from dataclasses import asdict

STAGES = (("n3", "N-3", "R-N3"), ("n2", "N-2", "R-N2"), ("n1", "N-1", "R-N1"))
DEFAULT_SIZES = (2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 3000, 5000, 10000, 12000, 15000, 20000, 25000)
SEED_KEYS = ("day", "Xv", "Xd", "viab", "glc", "gln", "lac", "amm", "osm", "mu", "kd", "net")
MODES = ("batch", "fedbatch", "perfusion")
MODE_LABELS = {"batch": "batch", "fedbatch": "fed-batch", "perfusion": "perfusion"}
PRESETS = {
    "workbook":     {"n3": "batch", "n2": "batch", "n1": "batch", "n": "fedbatch"},
    "fb_all":       {"n3": "fedbatch", "n2": "fedbatch", "n1": "fedbatch", "n": "fedbatch"},
    "fb_n1_perf_n": {"n3": "batch", "n2": "batch", "n1": "fedbatch", "n": "perfusion"},
    "perf_all":     {"n3": "perfusion", "n2": "perfusion", "n1": "perfusion", "n": "perfusion"},
}
LN2 = math.log(2.0)
_CACHE, _CACHE_MAX = {}, 64


# ── helpers ──────────────────────────────────────────────────────────
def _key(obj):
    return tuple((k, tuple(v) if isinstance(v, list) else v) for k, v in sorted(asdict(obj).items()))


def _cached(key, fn):
    if key not in _CACHE:
        if len(_CACHE) >= _CACHE_MAX:
            _CACHE.clear()
        _CACHE[key] = fn()
    return _CACHE[key]


def _monod(x, K):
    return x / (K + x) if K + x > 0 else 0.0


def _inhib(x, K):
    return K / (K + x) if K + x > 0 else 1.0


def _hill(x, K, n):
    """x^n / (K^n + x^n), written to avoid huge powers."""
    if x <= 0:
        return 0.0
    if K <= 0:
        return 1.0
    r = (x / K) ** n
    return r / (1.0 + r)


def _starve(G, K):
    d = K * K + G * G
    return K * K / d if d > 0 else 0.0


def _rnd(x, n=6):
    if x is None or isinstance(x, (bool, str, int)):
        return x
    if isinstance(x, float):
        return round(x, n) if math.isfinite(x) else None
    if isinstance(x, dict):
        return {k: _rnd(v, n) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_rnd(v, n) for v in x]
    return x


def _sizes(st):
    try:
        s = sorted(float(v) for v in (st.vessel_sizes_L or []) if float(v) > 0)
    except (TypeError, ValueError):
        s = []
    return s or [float(v) for v in DEFAULT_SIZES]


def _std(x, sizes):
    """Smallest standard size >= x (whole multiples of the largest above the list)."""
    for s in sizes:
        if s >= x - 1e-9:
            return s
    return math.ceil(x / sizes[-1] - 1e-9) * sizes[-1]


def _fmtL(v):
    return f"{v:,.0f}" if v >= 1 else f"{v:g}"


def _interp(a, b, f):
    return {k: (a[k] + (b[k] - a[k]) * f
                if isinstance(a[k], (int, float)) and not isinstance(a[k], bool) else a[k]) for k in a}


def _at(rows, t):
    """State at time t, interpolated between time steps."""
    if t <= rows[0]["day"]:
        return dict(rows[0])
    for i in range(1, len(rows)):
        if rows[i]["day"] >= t - 1e-12:
            a, b = rows[i - 1], rows[i]
            span = b["day"] - a["day"]
            return _interp(a, b, (t - a["day"]) / span if span > 0 else 1.0)
    return dict(rows[-1])


def basis_label(st):
    d, x = st.duration_basis == "biology", st.density_basis == "biology"
    if d and x:
        return "Biology (durations and transfer densities from Seed_sim)"
    if d:
        return "Biology + rules of thumb (durations from Seed_sim, fixed transfer densities)"
    if x:
        return "Biology + rules of thumb (transfer densities from Seed_sim, durations from fixed μ)"
    return "Rules of thumb (fixed transfer densities and fixed μ)"


def stage_modes(p):
    """Culture mode of each stage {"n3", "n2", "n1", "n"}: "batch", "fedbatch" or "perfusion".
    Without a culture block the workbook layout is used (batch seed stages, fed-batch N)."""
    c = getattr(p, "culture", None)
    base = PRESETS["workbook"]
    if c is None:
        return dict(base)
    preset = str(c.preset).lower()
    if preset in PRESETS:
        return dict(PRESETS[preset])
    out = {}
    for k, default in base.items():                       # "custom" (or any other name)
        m = str(getattr(c, "mode_" + k, default)).lower()
        out[k] = m if m in MODES else default
    return out


# ── the cell model: batch and fed-batch (Bio_sim and Seed_sim) ───────
def _simulate(b, X0, days, fed, shift):
    """Explicit Euler, exactly as the workbook: rates at step k from the state at step k."""
    dt = max(float(b.dt_d), 1e-4)
    n = max(1, int(round(days / dt)))
    V, Xv, Xd = float(b.V0_L), float(X0), 0.0
    G, Q, L, A = float(b.glc0_gL), float(b.gln0_mM), float(b.lac0_gL), float(b.amm0_mM)
    S, P, ivcd = float(b.S0_mOsm), 0.0, 0.0
    mu_max = max(b.mu_max_d, 1e-12)
    rows = []
    for k in range(n + 1):
        t = round(k * dt, 9)
        sh = shift and t >= b.shift_day - 1e-9
        osm = S + 5.55 * G + 22.2 * L + Q + A
        f_osm = b.K_osm / max(b.K_osm + max(0.0, osm - b.osm_crit), 1e-12)
        f_lag = 1.0 - math.exp(-t / b.tau_lag_d) if b.tau_lag_d > 0 else 1.0
        mu = (b.mu_max_d * (b.f_mu_shift if sh else 1.0) * f_lag
              * _monod(G, b.K_glc_gL) * _monod(Q, b.K_gln_mM)
              * _inhib(L, b.KI_lac_gL) * _inhib(A, b.KI_amm_mM) * f_osm)
        kd = ((b.f_kd_shift if sh else 1.0) * (b.kd_min_d + b.kd_amm_max_d * _hill(A, b.Kd_amm_mM, b.n_amm))
              + b.kd_osm_max_d * _hill(osm, b.Kd_osm, b.n_osm)
              + b.kd_starve_max_d * _starve(G, b.K_starve_gL))
        qP = (b.qP_pg_cell_d * (b.f_qp_shift if sh else 1.0)
              * (1.0 + b.qp_osm_rise * max(0.0, osm - b.osm_crit) / 100.0))
        fm = b.f_maint_shift if sh else 1.0
        F = b.CSFR_pL_cell_d * 1e-3 * Xv * V if (fed and t >= b.feed_start_d - 1e-9) else 0.0
        D = F / V if V > 0 else 0.0
        q_glc = mu * b.glc_per_cell + b.m_glc * fm
        q_gln = mu * b.gln_per_cell + b.m_gln * fm
        Xt = Xv + Xd
        rows.append({"day": t, "T_C": b.T_after_C if sh else b.T_before_C, "V_L": V, "F_Ld": F,
                     "Xv": Xv, "Xd": Xd, "Xt": Xt, "viab": 100.0 * Xv / Xt if Xt > 0 else 0.0,
                     "glc": G, "gln": Q, "lac": L, "amm": A, "S": S, "osm": osm, "f_osm": f_osm,
                     "mu": mu, "kd": kd, "net": mu - kd, "qP": qP, "ivcd": ivcd, "titre": P, "mass_g": P * V})
        if k == n:
            break
        dXv = (mu - kd - D) * Xv
        dXd = kd * Xv - (b.k_lysis_d + D) * Xd
        dG = D * (b.feed_glc_gL - G) - q_glc * Xv
        dQ = D * (b.feed_gln_mM - Q) - q_gln * Xv - b.k_deg_gln_d * Q
        dL = (b.Y_lac_glc * mu * b.glc_per_cell * Xv
              - b.q_lac_max * Xv * _monod(L, b.K_lac_gL) * max(0.0, 1.0 - mu / mu_max) - D * L)
        dA = b.Y_amm_gln * q_gln * Xv + b.k_deg_gln_d * Q - D * A
        dS = D * (b.S_feed_mOsm - S)
        dP = qP * Xv / 1000.0 - D * P
        Xv_new = max(0.0, Xv + dXv * dt)
        ivcd += (Xv + Xv_new) / 2.0 * dt
        Xv, Xd = Xv_new, max(0.0, Xd + dXd * dt)
        G, Q, L, A = max(0.0, G + dG * dt), max(0.0, Q + dQ * dt), max(0.0, L + dL * dt), max(0.0, A + dA * dt)
        S, P = max(0.0, S + dS * dt), max(0.0, P + dP * dt)
        V += F * dt
    return rows


# ── the cell model: perfusion ────────────────────────────────────────
def _simulate_perfusion(b, c, X0, days, bleed, shift):
    """Perfusion stage with the same cell biology as _simulate.

    Volume is constant. Fresh medium (culture.perf_*) enters at D (1/day) from culture.perf_start_d;
    live and dead cells are held back by the retention device; dissolved species leave with the
    outflow; product leaves in the permeate at its passage (sieve) and fully in the bleed.
    With bleed=True a cell bleed holds culture.perf_X_setpoint once the culture gets close to it.
    Rows every bio.dt_d (integrated in culture.perf_substeps Euler steps) with the _simulate keys:
    'titre' is the permeate titre and 'mass_g' the product harvested in the permeate so far.
    """
    dt_rec = max(float(b.dt_d), 1e-4)
    m = max(1, int(c.perf_substeps))
    dt = dt_rec / m
    n = max(1, int(round(days / dt_rec)))
    V = float(b.V0_L)
    Xv, Xd = float(X0), 0.0
    G, Q, L, A = float(b.glc0_gL), float(b.gln0_mM), float(b.lac0_gL), float(b.amm0_mM)
    S, P, ivcd = float(b.S0_mOsm), 0.0, 0.0
    mu_max = max(b.mu_max_d, 1e-12)
    Xsp = float(c.perf_X_setpoint)
    use_bleed = bool(bleed) and Xsp > 0
    tau = max(float(c.perf_bleed_tau_d), 2.0 * dt)
    keep_g = 1.0 - min(max(float(c.perf_qp_loss_per_gen), 0.0), 0.999)
    keep_d = 1.0 - min(max(float(c.perf_qp_loss_per_day), 0.0), 0.999)
    cspr = str(c.perf_rate_basis).lower() != "vvd"
    can_swap = float(c.perf_sieve0) > float(c.perf_swap_below)
    sieve, swaps, gens = float(c.perf_sieve0), 0, 0.0
    harv = lost = media = 0.0                     # per litre of working volume
    bleed_on = False
    rows = []

    def rates(t):
        """Specific rates and flows from the current state (read from the enclosing function)."""
        sh = shift and t >= b.shift_day - 1e-9
        osm = S + 5.55 * G + 22.2 * L + Q + A
        f_osm = b.K_osm / max(b.K_osm + max(0.0, osm - b.osm_crit), 1e-12)
        f_lag = 1.0 - math.exp(-t / b.tau_lag_d) if b.tau_lag_d > 0 else 1.0
        f_cap = max(0.0, 1.0 - Xv / c.perf_X_max) if c.perf_X_max > 0 else 1.0
        mu = (b.mu_max_d * (b.f_mu_shift if sh else 1.0) * f_lag * f_cap
              * _monod(G, b.K_glc_gL) * _monod(Q, b.K_gln_mM)
              * _inhib(L, b.KI_lac_gL) * _inhib(A, b.KI_amm_mM) * f_osm)
        kd = ((b.f_kd_shift if sh else 1.0) * (b.kd_min_d + b.kd_amm_max_d * _hill(A, b.Kd_amm_mM, b.n_amm))
              + b.kd_osm_max_d * _hill(osm, b.Kd_osm, b.n_osm)
              + b.kd_starve_max_d * _starve(G, b.K_starve_gL))
        qP = (b.qP_pg_cell_d * (b.f_qp_shift if sh else 1.0)
              * (1.0 + b.qp_osm_rise * max(0.0, osm - b.osm_crit) / 100.0)
              * keep_g ** gens * keep_d ** t)                       # productivity ageing
        fm = b.f_maint_shift if sh else 1.0
        D = B = 0.0
        if t >= c.perf_start_d - 1e-9:
            D = (min(max(c.perf_cspr_pL_cell_d * 1e-3 * Xv, c.perf_vvd_min), c.perf_vvd_max)
                 if cspr else max(float(c.perf_vvd), 0.0))
            if bleed_on:
                # replace what the cells grow, plus a pull back towards the setpoint
                B = min(max(mu - kd + (Xv - Xsp) / (Xsp * tau), 0.0), max(float(c.perf_bleed_max_vvd), 0.0))
                D = max(D, B)
        return {"sh": sh, "osm": osm, "f_osm": f_osm, "mu": mu, "kd": kd, "qP": qP, "D": D, "B": B,
                "q_glc": mu * b.glc_per_cell + b.m_glc * fm, "q_gln": mu * b.gln_per_cell + b.m_gln * fm}

    for k in range(n + 1):
        t = round(k * dt_rec, 9)
        r = rates(t)
        Xt = Xv + Xd
        rows.append({"day": t, "T_C": b.T_after_C if r["sh"] else b.T_before_C, "V_L": V, "F_Ld": r["D"] * V,
                     "Xv": Xv, "Xd": Xd, "Xt": Xt, "viab": 100.0 * Xv / Xt if Xt > 0 else 0.0,
                     "glc": G, "gln": Q, "lac": L, "amm": A, "S": S, "osm": r["osm"], "f_osm": r["f_osm"],
                     "mu": r["mu"], "kd": r["kd"], "net": r["mu"] - r["kd"], "qP": r["qP"], "ivcd": ivcd,
                     "titre": sieve * P, "mass_g": harv * V,
                     "titre_reactor": P, "D_vvd": r["D"], "bleed_vvd": r["B"], "sieve": sieve,
                     "gens": gens, "tff_swaps": swaps, "media_L": media * V, "bleed_loss_g": lost * V})
        if k == n:
            break
        for j in range(m):
            ts = t + j * dt
            if use_bleed and not bleed_on and Xv >= c.perf_bleed_enable_frac * Xsp:
                bleed_on = True
            r = rates(ts)
            mu, kd, D, B = r["mu"], r["kd"], r["D"], r["B"]
            Dp = D - B                                              # permeate share of the outflow
            dXv = (mu - kd - B) * Xv
            dXd = kd * Xv - (b.k_lysis_d + B) * Xd
            dG = D * (c.perf_glc_gL - G) - r["q_glc"] * Xv
            dQ = D * (c.perf_gln_mM - Q) - r["q_gln"] * Xv - b.k_deg_gln_d * Q
            dL = (b.Y_lac_glc * mu * b.glc_per_cell * Xv
                  - b.q_lac_max * Xv * _monod(L, b.K_lac_gL) * max(0.0, 1.0 - mu / mu_max) - D * L)
            dA = b.Y_amm_gln * r["q_gln"] * Xv + b.k_deg_gln_d * Q - D * A
            dS = D * (c.perf_S_mOsm - S)
            dP = r["qP"] * Xv / 1000.0 - (Dp * sieve + B) * P
            harv += Dp * sieve * P * dt
            lost += B * P * dt
            media += D * dt
            gens += max(mu, 0.0) * dt / LN2
            Xv_new = max(0.0, Xv + dXv * dt)
            ivcd += (Xv + Xv_new) / 2.0 * dt
            Xv, Xd = Xv_new, max(0.0, Xd + dXd * dt)
            G, Q, L, A = (max(0.0, G + dG * dt), max(0.0, Q + dQ * dt),
                          max(0.0, L + dL * dt), max(0.0, A + dA * dt))
            S, P = max(0.0, S + dS * dt), max(0.0, P + dP * dt)
            if ts >= c.perf_start_d - 1e-9:                         # retention device in use: it fouls
                sieve = max(0.0, sieve - (c.perf_foul_base_d + c.perf_foul_per_dead_d * Xd) * dt)
                if can_swap and sieve < c.perf_swap_below:
                    sieve, swaps = float(c.perf_sieve0), swaps + 1
    return rows


def _results(rows, b, mode="fedbatch"):
    last = rows[-1]
    pk = max(rows, key=lambda r: r["Xv"])
    lp = max(rows, key=lambda r: r["lac"])
    hit = next((r["day"] for r in rows if r["titre"] >= b.target_titre_gL - 1e-12), None)
    out = {"mode": mode,
           "final_titre_gL": last["titre"], "final_mass_kg": last["mass_g"] / 1000.0, "target_day": hit,
           "peak_vcd": pk["Xv"], "peak_vcd_day": pk["day"], "final_viability_pct": last["viab"],
           "min_glucose_gL": min(r["glc"] for r in rows), "peak_lactate_gL": lp["lac"],
           "peak_lactate_day": lp["day"], "final_lactate_gL": last["lac"],
           "peak_ammonia_mM": max(r["amm"] for r in rows), "peak_osm": max(r["osm"] for r in rows),
           "final_volume_L": last["V_L"], "final_ivcd": last["ivcd"], "final_qP": last["qP"]}
    if mode == "perfusion":
        perf = [r for r in rows if r["D_vvd"] > 0]
        out.update({"harvest_g_per_L": last["mass_g"] / max(last["V_L"], 1e-9),
                    "bleed_loss_g": last["bleed_loss_g"], "media_L": last["media_L"],
                    "mean_vvd": sum(r["D_vvd"] for r in perf) / len(perf) if perf else 0.0,
                    "generations": last["gens"], "final_sieve": last["sieve"], "tff_swaps": last["tff_swaps"]})
    return out


def bio_run(p):
    """Bio_sim for N in its culture mode (cached on the parameters that matter)."""
    b = p.bio
    mode = stage_modes(p)["n"]
    c = getattr(p, "culture", None)
    if mode == "perfusion" and c is not None:
        def go_perf():
            rows = _simulate_perfusion(b, c, b.X0, b.run_days, bleed=True, shift=bool(c.perf_temp_shift))
            return {"rows": rows, "results": _results(rows, b, "perfusion")}
        return _cached(("bio", "perfusion") + _key(b) + _key(c), go_perf)
    fed = mode != "batch"

    def go():
        rows = _simulate(b, b.X0, b.run_days, fed=fed, shift=True)
        return {"rows": rows, "results": _results(rows, b, "fedbatch" if fed else "batch")}
    return _cached(("bio", mode) + _key(b), go)


def seed_sim(b, X0, days, mode="batch", c=None):
    """Seed_sim: one seed stage at 37 °C, same biology. mode "batch" is the workbook's stage."""
    if mode == "perfusion" and c is not None:
        return _cached(("seed", "perfusion", float(X0), float(days)) + _key(b) + _key(c),
                       lambda: _simulate_perfusion(b, c, X0, days, bleed=False, shift=False))
    fed = mode == "fedbatch"
    return _cached(("seed", "fedbatch" if fed else "batch", float(X0), float(days)) + _key(b),
                   lambda: _simulate(b, X0, days, fed=fed, shift=False))


def bio_titre_points(p):
    """Titre curve for upstream.molecule = "biology" ([N day, g/L]); titre_scale is applied by titre_at."""
    run = bio_run(p)
    rows = run["rows"]
    if run["results"].get("mode") == "perfusion":
        # The twin multiplies this curve by its own permeate flow. Rescale the permeate titre so the
        # grams per day match the cell model even when the two perfusion rates differ.
        up = p.upstream
        out_frac = up.bleed_frac if up.cell_retention else up.tff_retentate_frac
        twin_perm = max(float(up.VVD) * (1.0 - out_frac), 1e-9)
        return [[r["day"], r["titre"] * (r["D_vvd"] - r["bleed_vvd"]) / twin_perm] for r in rows]
    if getattr(p.bio, "titre_link", "fedbatch") == "perfusion":
        VVD = max(float(p.upstream.VVD), 1e-9)
        return [[r["day"], r["qP"] * r["Xv"] / (1000.0 * VVD)] for r in rows]   # Inputs sheet: P = qP·Xv/D
    return [[r["day"], r["titre"]] for r in rows]


def transfer_point(rows, threshold):
    """Last moment net growth is still >= threshold × its peak (late exponential phase),
    interpolated between the time steps either side, as the workbook does."""
    net = [r["net"] for r in rows]
    ipk = max(range(len(net)), key=net.__getitem__)
    pk = net[ipk]
    lim = threshold * pk
    for i in range(ipk + 1, len(rows)):
        if net[i] < lim:
            d = net[i - 1] - net[i]
            f = (net[i - 1] - lim) / d if d > 0 else 1.0
            return dict(_interp(rows[i - 1], rows[i], min(max(f, 0.0), 1.0)), found=True, peak_net=pk)
    return dict(rows[-1], found=False, peak_net=pk)


def time_to_reach(rows, X):
    """First time the simulated viable density reaches X (None if it never does)."""
    if rows[0]["Xv"] >= X:
        return rows[0]["day"]
    for i in range(1, len(rows)):
        a, b = rows[i - 1], rows[i]
        if b["Xv"] >= X:
            d = b["Xv"] - a["Xv"]
            return a["day"] + (b["day"] - a["day"]) * ((X - a["Xv"]) / d if d > 0 else 1.0)
    return None


# ── Demand: production reactor N from annual demand ──────────────────
def demand(p, run=None):
    d, b, st = p.demand, p.bio, p.seedtrain
    R = (run or bio_run(p))["results"]
    perf = R.get("mode") == "perfusion"
    titre = R["final_titre_gL"]
    if perf:
        g_per_L = R["harvest_g_per_L"]          # grams harvested per litre of reactor over the run
        ratio = 1.0                             # constant volume
    else:
        g_per_L = titre
        ratio = R["final_volume_L"] / max(b.V0_L, 1e-9)
    kg_m3 = g_per_L * d.dsp_yield
    cycle = b.run_days + d.turnaround_d
    bpr = int(math.floor(d.operating_days / cycle + 1e-9)) if cycle > 0 else 0
    sizes, fmax = _sizes(st), max(st.max_fill_n, 1e-9)
    if perf:
        notes = [f"Perfusion basis: grams harvested in the permeate per litre of reactor over {b.run_days:g} days, "
                 "from the perfusion cell model; purification yield, batch success and turnaround are rules of thumb."]
    else:
        notes = ["Fed-batch basis: harvest titre and volume growth from Bio_sim; purification yield, batch success "
                 "and turnaround are rules of thumb. The IDB twin's own per-year totals are on the Outputs tab."]
    if bpr <= 0:
        notes.append("No whole batch fits the operating year with this run length and turnaround.")
    if kg_m3 <= 0:
        notes.append("The cell model makes no antibody, so N cannot be sized from demand.")

    def option(n):
        bpy = bpr * n
        if bpy <= 0 or kg_m3 <= 0 or d.batch_success <= 0:
            return {"reactors": n, "batches_per_year": bpy}
        kg_b = d.annual_kg / (bpy * d.batch_success)
        full = kg_b / kg_m3 * 1000.0
        ves = _std(full / fmax, sizes)
        return {"reactors": n, "batches_per_year": bpy, "kg_per_batch": kg_b, "full_L": full,
                "start_L": full / ratio if ratio > 0 else None, "vessel_L": ves, "fill_pct": 100.0 * full / ves,
                "capacity_kg_yr": bpy * d.batch_success * ves * fmax / 1000.0 * kg_m3}

    return {"basis": "perfusion" if perf else "fedbatch", "annual_kg": d.annual_kg, "titre_gL": titre,
            "g_per_L_run": g_per_L, "volume_ratio": ratio, "kg_per_m3": kg_m3,
            "cycle_d": cycle, "batches_per_reactor": bpr, "selected": option(max(1, int(d.n_reactors))),
            "options": [option(n) for n in range(1, 5)], "notes": notes}


# ── Seed_train: N-3, N-2, N-1 and N back to the vial ─────────────────
def _chk(ok, text):
    return {"ok": bool(ok), "text": text}


def seed_train(p, n_basis=None, run=None):
    b, st, br = p.bio, p.seedtrain, p.bioreactor
    c = getattr(p, "culture", None)
    modes = stage_modes(p)
    run = run or bio_run(p)
    basis = n_basis or st.n_basis
    dem = demand(p, run)
    sizes = _sizes(st)
    thr = st.growth_threshold_pct / 100.0
    dur_bio, den_bio = st.duration_basis == "biology", st.density_basis == "biology"
    notes = []

    # production reactor N
    sel = dem["selected"]
    if basis == "demand" and sel.get("start_L"):
        V_N, ves_N, src, run_N = sel["start_L"], sel["vessel_L"], "Demand-based (Demand sheet)", b.run_days
    elif basis == "manual":
        V_N, ves_N, src, run_N = st.manual_start_L, st.manual_vessel_L, "Manual", b.run_days
    else:
        if basis == "demand":
            notes.append("Demand could not size N; the twin N reactor (bioreactor.V_L) is used instead.")
        V_N = br.V_L
        ves_N = _std(V_N / max(st.max_fill_n, 1e-9), sizes)
        src, run_N, basis = "Twin N reactor (bioreactor.V_L)", p.run_days, "twin"
    X0_N = b.X0

    # each seed stage from its own simulation in its culture mode
    stages = []
    for key, name, tag in STAGES:
        mode = modes[key]
        x0 = float(getattr(st, "seed_" + key))
        sim = seed_sim(b, x0, st.sim_days, mode, c)
        tp = transfer_point(sim, thr)
        xt_rule = float(getattr(st, "xt_rule_" + key))
        xt = tp["Xv"] if den_bio else xt_rule
        t_bio = time_to_reach(sim, xt)
        t_rule = math.log(xt / x0) / st.mu_rule_d if (x0 > 0 and xt > x0 and st.mu_rule_d > 0) else 0.0
        flag = ""
        if dur_bio and t_bio is None:
            flag = (f"{name}: Seed_sim never reaches {xt:g} ×10⁶/mL in {st.sim_days:g} d; "
                    f"the rule-of-thumb duration is used")
            notes.append(flag)
        if not tp["found"]:
            notes.append(f"{name}: net growth never fell to {st.growth_threshold_pct:g} % of its peak in "
                         f"{st.sim_days:g} d; the biology transfer density is the end of the simulation.")
        dur = t_bio if (dur_bio and t_bio is not None) else t_rule
        at = _at(sim, t_bio if t_bio is not None else sim[-1]["day"])
        end = _at(sim, dur)
        perf_per_L = end.get("media_L", 0.0) / max(b.V0_L, 1e-9) if mode == "perfusion" else 0.0
        pk = tp["peak_net"]
        stages.append({"key": key, "name": name, "tag": tag, "type": getattr(st, "type_" + key), "mode": mode,
                       "max_fill": float(getattr(st, "max_fill_" + key)),
                       "min_fill": float(getattr(st, "min_fill_" + key)),
                       "x0": x0, "xt_bio": tp["Xv"], "xt_rule": xt_rule, "xt": xt,
                       "dur_bio": t_bio, "dur_rule": t_rule, "duration_d": dur,
                       "net_pct": 100.0 * at["net"] / pk if pk > 0 else None,
                       "viab_pct": at["viab"], "gln_mM": at["gln"], "glc_gL": at["glc"],
                       "peak_net": pk, "threshold_net": thr * pk,
                       "expansion": xt / x0 if x0 > 0 else None, "flag": flag,
                       "vol_growth": end["V_L"] / max(b.V0_L, 1e-9),      # fed-batch feed (1 otherwise)
                       "perf_L_per_L": perf_per_L,                       # perfused medium per L of stage
                       "perf_vvd_avg": perf_per_L / dur if dur > 0 else 0.0})

    # volumes back-calculated from N: cells needed downstream ÷ transfer density (volume at transfer)
    V_next, X0_next = V_N, X0_N
    for s in reversed(stages):
        s["V_L"] = X0_next * V_next / max(s["xt"], 1e-12)
        s["start_L"] = s["V_L"] / max(s["vol_growth"], 1e-9)            # fed-batch starts smaller
        s["media_L"] = s["V_L"] * (1.0 + s["perf_L_per_L"])              # fill + feed or + perfusion
        s["cells_1e9"] = s["xt"] * s["V_L"]              # 1e6 cells/mL × L = 1e9 cells
        s["split"] = V_next / s["V_L"] if s["V_L"] > 0 else None
        s["vessel_L"] = _std(s["V_L"] / max(s["max_fill"], 1e-9), sizes)
        s["fill_pct"] = 100.0 * s["V_L"] / s["vessel_L"]
        V_next, X0_next = s["start_L"], s["x0"]

    # oxygen: OUR at transfer density (seed) or peak VCD of the cell model (N)
    dC = st.O2_sat_air_mM * (1.0 - st.DO_pct / 100.0)
    for s in stages:
        s["OUR"] = s["xt"] * st.qO2_pmol_cell_h
        s["kLa_air"] = s["OUR"] / dC if dC > 0 else None
    peak_N = run["results"]["peak_vcd"]
    N = {"key": "n", "name": "N (production)", "tag": "R-N", "type": "Stirred tank", "source": src,
         "mode": modes["n"], "max_fill": st.max_fill_n, "min_fill": st.min_fill_n, "x0": X0_N, "V_L": V_N,
         "vessel_L": ves_N, "fill_pct": 100.0 * V_N / ves_N if ves_N else None, "duration_d": run_N,
         "peak_vcd": peak_N, "OUR": peak_N * st.qO2_pmol_cell_h,
         "kLa_air": peak_N * st.qO2_pmol_cell_h / dC if dC > 0 else None}

    # timeline from vial thaw
    s3 = stages[0]
    cells = s3["x0"] * s3["start_L"]
    doublings = max(0.0, math.log2(cells * 1e3 / max(st.vial_cells_1e6, 1e-9))) if cells > 0 else 0.0
    if dur_bio and s3["dur_bio"] and s3["dur_bio"] > 0 and s3["xt"] > s3["x0"]:
        mu_f = math.log(s3["xt"] / s3["x0"]) / s3["dur_bio"]
        mu_src = "average growth in the N-3 simulation, lag included"
    else:
        mu_f, mu_src = st.mu_rule_d, "fixed μ (rule of thumb)"
    flask = doublings * math.log(2.0) / mu_f if mu_f > 0 else 0.0
    t = flask
    for s in stages:
        s["start_d"], t = t, t + s["duration_d"]
        s["end_d"] = t
    N["start_d"], N["end_d"] = t, t + run_N
    tl = {"cells_n3_1e9": cells, "doublings": doublings, "mu_flask_d": mu_f, "mu_source": mu_src,
          "flask_d": flask, "vial_to_harvest_d": N["end_d"]}

    # checks
    cols = stages + [N]

    def c_fill(s):
        f = s.get("fill_pct")
        if f is None:
            return None
        f0 = 100.0 * s.get("start_L", s["V_L"]) / s["vessel_L"] if s.get("vessel_L") else f
        text = f"{f0:.0f}–{f:.0f} %" if abs(f - f0) > 0.5 else f"{f:.0f} %"
        return _chk(s["min_fill"] * 100 - 1e-6 <= f0 and f <= s["max_fill"] * 100 + 1e-6, text)

    def c_split(s):
        x = s.get("split")
        return None if x is None else _chk(st.split_min - 1e-9 <= x <= st.split_max + 1e-9, f"1:{x:.1f}")

    def c_growth(s):
        x = s.get("net_pct")
        return None if x is None else _chk(x >= st.growth_threshold_pct - 0.5, f"{x:.0f} % of peak")

    def c_via(s):
        x = s.get("viab_pct")
        return None if x is None else _chk(x >= st.viability_min_pct, f"{x:.1f} %")

    def c_len(s):
        if s["key"] == "n":
            return None
        x = s["duration_d"]
        return _chk(st.stage_min_d - 1e-9 <= x <= st.stage_max_d + 1e-9, f"{x:.2f} d")

    checks = [{"label": lab, "cells": [fn(s) for s in cols]} for lab, fn in (
        ("Fill within vessel limits", c_fill),
        (f"Split ratio (1:{st.split_min:g} to 1:{st.split_max:g})", c_split),
        (f"Cells still growing fast at transfer (≥ {st.growth_threshold_pct:g} % of peak)", c_growth),
        (f"Viability at transfer ≥ {st.viability_min_pct:g} %", c_via),
        (f"Stage length ({st.stage_min_d:g}–{st.stage_max_d:g} days)", c_len))]
    fails = [f"{cols[i]['name']}: {c['label']} ({cell['text']})"
             for c in checks for i, cell in enumerate(c["cells"]) if cell and not cell["ok"]]
    chain = " → ".join(_fmtL(s["vessel_L"]) + " L" for s in cols)
    overall = (f"Seed train is consistent: {chain}, {tl['vial_to_harvest_d']:.0f} days from vial to harvest"
               if not fails else
               f"Seed train needs attention ({len(fails)} check{'s' if len(fails) > 1 else ''}): "
               + "; ".join(fails[:4]))
    return {"basis": {"n": basis, "duration": st.duration_basis, "density": st.density_basis},
            "basis_label": basis_label(st), "threshold_pct": st.growth_threshold_pct, "modes": modes,
            "stages": stages, "n": N, "timeline": tl, "checks": checks, "ok": not fails,
            "overall": overall, "notes": notes, "demand": dem}


def apply_seed_train_to_n1(p):
    """seedtrain.apply_to_n1: run the twin's N-1 as the N-1 stage of a seed train sized for the
    twin N reactor. Idempotent (the train does not depend on the n1 parameters). Returns the train or None.
    A perfusion N-1 keeps its average perfusion rate (permeate to waste); batch and fed-batch run at VVD 0."""
    st = getattr(p, "seedtrain", None)
    if st is None or not st.apply_to_n1:
        return None
    tr = seed_train(p, n_basis="twin")
    s = tr["stages"][-1]
    n1 = p.n1
    n1.enabled = True
    n1.V_L = n1.transfer_L = float(s["V_L"])
    n1.VVD = float(s["perf_vvd_avg"]) if s["mode"] == "perfusion" else 0.0
    n1.X0 = float(s["x0"])
    n1.X_max = max(float(s["xt"]), n1.X0)
    days = max(float(s["duration_d"]), 1.0 / 24.0)
    n1.days, n1.target_seed_vcd = days, 0.0                # fixed duration from the chosen basis
    n1.max_days = max(n1.max_days, days)
    n1.mu_h = math.log(n1.X_max / n1.X0) / (days * 24.0) if n1.X_max > n1.X0 else 0.0
    return tr


def bio_design(p):
    """Everything the page shows: N curves and results, seed-stage curves, seed train and demand."""
    b, st = p.bio, p.seedtrain
    c = getattr(p, "culture", None)
    modes = stage_modes(p)
    run = bio_run(p)
    rows, R = run["rows"], run["results"]
    step = max(1, int(round(b.report_days / max(b.dt_d, 1e-9))))
    tr = dict(seed_train(p, "twin" if st.apply_to_n1 else None, run))
    dem = tr.pop("demand")
    seeds = {}
    for s in tr["stages"]:
        sim = seed_sim(b, s["x0"], st.sim_days, s["mode"], c)
        keys = SEED_KEYS + (("D_vvd",) if s["mode"] == "perfusion" else ())
        seeds[s["key"]] = {"name": s["name"], "mode": s["mode"], "x0": s["x0"], "xt": s["xt"],
                           "transfer_day": s["dur_bio"], "peak_net": s["peak_net"],
                           "threshold_net": s["threshold_net"],
                           "rows": [{k: r[k] for k in keys} for r in sim]}
    preset = str(getattr(c, "preset", "workbook")) if c is not None else "workbook"
    mode_txt = ", ".join(f"{name} {MODE_LABELS[modes[key]]}" for key, name, _ in STAGES)
    notes = [f"Seed train basis: {tr['basis_label']}. N from: {tr['n']['source']}.",
             f"Culture modes: {mode_txt}, N {MODE_LABELS[modes['n']]} (culture.preset = {preset})."]
    if modes["n"] == "perfusion":
        notes.append(f"N is a perfusion reactor: constant volume, cells held back by the retention device, a bleed "
                     f"holds {c.perf_X_setpoint:g} ×10⁶ cells/mL. Titre is the permeate titre; Demand sizes N from "
                     f"the grams harvested per litre of reactor over {b.run_days:g} days.")
        if c.perf_qp_loss_per_gen <= 0 and c.perf_qp_loss_per_day <= 0:
            notes.append("No productivity ageing is set (culture.perf_qp_loss_per_gen and perf_qp_loss_per_day "
                         "are 0): cells keep their productivity for the whole run. Set the loss from the cell "
                         "line stability study.")
        if c.perf_foul_per_dead_d <= 0:
            notes.append("Retention-device fouling runs at a fixed rate (culture.perf_foul_per_dead_d = 0); it "
                         "does not yet speed up when viability falls.")
    else:
        notes.append(f"N runs as {MODE_LABELS[modes['n']]} (CSFR feeding, temperature shift on day "
                     f"{b.shift_day:g}); Demand sizes N on that basis. The twin's IDB N reactor is a perfusion "
                     "reactor: set culture.preset to fb_n1_perf_n or perf_all to model N as perfusion.")
    if modes["n1"] == "perfusion":
        notes.append("A perfusion N-1 transfers at high density: raise bio.X0 (the N seed density) to suit, "
                     "or the split-ratio check will fail.")
    notes += ["Rules of thumb used whatever the basis: seed densities, the growth threshold, qO2, purification "
              "yield, batch success, turnaround and the check limits.",
              "Transfer-point values are interpolated at the transfer time (the workbook shows the next 0.1-day row)."]
    if p.upstream.molecule == "biology":
        if modes["n"] == "perfusion":
            link = "perfusion cell model, permeate titre rescaled to the twin's permeate flow (grams per day kept)"
        else:
            link = "perfusion qP·Xv/VVD" if b.titre_link == "perfusion" else "fed-batch titre"
        notes.append(f"N titre curve in this run: biology model, {link}, × titre_scale {p.upstream.titre_scale:.3g}.")
    if st.apply_to_n1:
        notes.append(f"The twin's N-1 runs as the {MODE_LABELS[modes['n1']]} N-1 stage of this train (volume, "
                     "duration and seed density from the train); the train is sized for the twin N reactor.")
    notes += tr["notes"]
    return _rnd({"bio": {"mode": modes["n"], "results": R, "series": rows, "curve": rows[::step],
                         "target_gL": b.target_titre_gL, "run_days": b.run_days, "dt_d": b.dt_d},
                 "seed_sim": seeds, "seed_train": tr, "demand": dem, "notes": notes})


if __name__ == "__main__":
    from icb_process import ICBParams
    d = bio_design(ICBParams())                       # default preset "workbook": must match the workbook
    R, T, D = d["bio"]["results"], d["seed_train"], d["demand"]
    for name, got, want in (("Final titre g/L", R["final_titre_gL"], 3.69008),
                            ("Peak VCD 1e6/mL", R["peak_vcd"], 13.02022),
                            ("Day of peak VCD", R["peak_vcd_day"], 11.6),
                            ("Final viability %", R["final_viability_pct"], 71.8422),
                            ("Day target reached", R["target_day"], 12.5),
                            ("Final volume L", R["final_volume_L"], 1777.294),
                            ("Transfer density", T["stages"][0]["xt_bio"], 3.896339),
                            ("Stage duration d", T["stages"][0]["dur_bio"], 4.369568),
                            ("Harvest volume L", D["selected"]["full_L"], 8007.294),
                            ("N starting volume L", D["selected"]["start_L"], 6307.461),
                            ("Vial to harvest d", T["timeline"]["vial_to_harvest_d"], 39.6083)):
        ok = got is not None and abs(got - want) <= 1e-3 * max(abs(want), 1.0)
        print(f"{'OK  ' if ok else 'DIFF'} {name:<22}{got!s:>14}   workbook {want}")
    print(T["overall"])

    print("\nCulture presets (typical perfusion defaults, not site data):")
    for label, preset, days, loss in (("fb_all", "fb_all", None, 0.0),
                                      ("fb_n1_perf_n", "fb_n1_perf_n", None, 0.0),
                                      ("perf_all", "perf_all", None, 0.0),
                                      ("perf_all 60 d, 0.5 %/gen", "perf_all", 60.0, 0.005)):
        q = ICBParams()
        q.culture.preset = preset
        q.culture.perf_qp_loss_per_gen = loss
        if days:
            q.bio.run_days = days
        dq = bio_design(q)
        r = dq["bio"]["results"]
        extra = ""
        if r["mode"] == "perfusion":
            extra = (f", {r['harvest_g_per_L']:.1f} g harvested per L reactor, mean {r['mean_vvd']:.2f} VVD, "
                     f"{r['generations']:.0f} generations, {r['tff_swaps']} device swap(s)")
        print(f"  {label:<26} N {r['mode']:<9} titre {r['final_titre_gL']:.2f} g/L, "
              f"viability {r['final_viability_pct']:.0f} %{extra}")
        print(f"  {'':<26} {dq['seed_train']['overall']}")
