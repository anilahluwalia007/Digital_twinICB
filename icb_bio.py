"""
icb_bio.py — cell biology, seed train (N-3 → N) and production-reactor sizing
NavNevAI_mAbWasteWater · CarbonTesseract 2026

Port of CHO_mAb_Process_Calculator5.xlsx (same equations, same default values):
  Bio_inputs / Bio_sim / Bio_curve  bio_run(p)          fed-batch cell model, Euler steps of bio.dt_d
  Seed_sim                          seed_sim(b, X0, d)  batch stage at 37 °C in fresh medium, same biology
  Seed_train                        seed_train(p)       N-3, N-2, N-1 and N back to the vial
  Demand                            demand(p)           N size from annual demand and the cell model's titre

Basis switches (workbook Seed_train B6 / B7), each "biology" or "rule":
  seedtrain.duration_basis  stage durations from Seed_sim (lag, nutrients, waste) or ln(fold) / mu_rule_d
  seedtrain.density_basis   transfer where net growth falls to growth_threshold_pct of its peak
                            (late exponential phase) or the fixed xt_rule_* densities
Both "biology" = pure biology; mixed = biology + rules of thumb; both "rule" = rules of thumb only.
Rules of thumb used whatever the basis (as in the workbook README): seed densities, the growth
threshold, qO2, purification yield, batch success, turnaround and every check limit.

Twin hooks
  bio_titre_points(p)        upstream.molecule = "biology": N titre curve from Bio_sim
  apply_seed_train_to_n1(p)  seedtrain.apply_to_n1 = True: the twin's N-1 runs as the batch N-1 stage
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


# ── the cell model (Bio_sim and Seed_sim) ────────────────────────────
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


def _results(rows, b):
    last = rows[-1]
    pk = max(rows, key=lambda r: r["Xv"])
    lp = max(rows, key=lambda r: r["lac"])
    hit = next((r["day"] for r in rows if r["titre"] >= b.target_titre_gL - 1e-12), None)
    return {"final_titre_gL": last["titre"], "final_mass_kg": last["mass_g"] / 1000.0, "target_day": hit,
            "peak_vcd": pk["Xv"], "peak_vcd_day": pk["day"], "final_viability_pct": last["viab"],
            "min_glucose_gL": min(r["glc"] for r in rows), "peak_lactate_gL": lp["lac"],
            "peak_lactate_day": lp["day"], "final_lactate_gL": last["lac"],
            "peak_ammonia_mM": max(r["amm"] for r in rows), "peak_osm": max(r["osm"] for r in rows),
            "final_volume_L": last["V_L"], "final_ivcd": last["ivcd"], "final_qP": last["qP"]}


def bio_run(p):
    """Bio_sim: fed-batch run of the biologically grounded model (cached on the bio parameters)."""
    b = p.bio

    def go():
        rows = _simulate(b, b.X0, b.run_days, fed=True, shift=True)
        return {"rows": rows, "results": _results(rows, b)}
    return _cached(("bio",) + _key(b), go)


def seed_sim(b, X0, days):
    """Seed_sim: one seed stage, batch at 37 °C in fresh medium, same biology."""
    return _cached(("seed", float(X0), float(days)) + _key(b),
                   lambda: _simulate(b, X0, days, fed=False, shift=False))


def bio_titre_points(p):
    """Titre curve for upstream.molecule = "biology" ([N day, g/L]); titre_scale is applied by titre_at."""
    rows = bio_run(p)["rows"]
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
    titre = R["final_titre_gL"]
    ratio = R["final_volume_L"] / max(b.V0_L, 1e-9)
    kg_m3 = titre * d.dsp_yield
    cycle = b.run_days + d.turnaround_d
    bpr = int(math.floor(d.operating_days / cycle + 1e-9)) if cycle > 0 else 0
    sizes, fmax = _sizes(st), max(st.max_fill_n, 1e-9)
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

    return {"annual_kg": d.annual_kg, "titre_gL": titre, "volume_ratio": ratio, "kg_per_m3": kg_m3,
            "cycle_d": cycle, "batches_per_reactor": bpr, "selected": option(max(1, int(d.n_reactors))),
            "options": [option(n) for n in range(1, 5)], "notes": notes}


# ── Seed_train: N-3, N-2, N-1 and N back to the vial ─────────────────
def _chk(ok, text):
    return {"ok": bool(ok), "text": text}


def seed_train(p, n_basis=None, run=None):
    b, st, br = p.bio, p.seedtrain, p.bioreactor
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

    # each seed stage from its own batch simulation
    stages = []
    for key, name, tag in STAGES:
        x0 = float(getattr(st, "seed_" + key))
        sim = seed_sim(b, x0, st.sim_days)
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
        pk = tp["peak_net"]
        stages.append({"key": key, "name": name, "tag": tag, "type": getattr(st, "type_" + key),
                       "max_fill": float(getattr(st, "max_fill_" + key)),
                       "min_fill": float(getattr(st, "min_fill_" + key)),
                       "x0": x0, "xt_bio": tp["Xv"], "xt_rule": xt_rule, "xt": xt,
                       "dur_bio": t_bio, "dur_rule": t_rule, "duration_d": dur,
                       "net_pct": 100.0 * at["net"] / pk if pk > 0 else None,
                       "viab_pct": at["viab"], "gln_mM": at["gln"], "glc_gL": at["glc"],
                       "peak_net": pk, "threshold_net": thr * pk,
                       "expansion": xt / x0 if x0 > 0 else None, "flag": flag})

    # volumes back-calculated from N: cells needed downstream ÷ transfer density
    V_next, X0_next = V_N, X0_N
    for s in reversed(stages):
        s["V_L"] = X0_next * V_next / max(s["xt"], 1e-12)
        s["cells_1e9"] = s["xt"] * s["V_L"]              # 1e6 cells/mL × L = 1e9 cells
        s["split"] = V_next / s["V_L"] if s["V_L"] > 0 else None
        s["vessel_L"] = _std(s["V_L"] / max(s["max_fill"], 1e-9), sizes)
        s["fill_pct"] = 100.0 * s["V_L"] / s["vessel_L"]
        V_next, X0_next = s["V_L"], s["x0"]

    # oxygen: OUR at transfer density (seed) or peak VCD of the cell model (N)
    dC = st.O2_sat_air_mM * (1.0 - st.DO_pct / 100.0)
    for s in stages:
        s["OUR"] = s["xt"] * st.qO2_pmol_cell_h
        s["kLa_air"] = s["OUR"] / dC if dC > 0 else None
    peak_N = run["results"]["peak_vcd"]
    N = {"key": "n", "name": "N (production)", "tag": "R-N", "type": "Stirred tank", "source": src,
         "max_fill": st.max_fill_n, "min_fill": st.min_fill_n, "x0": X0_N, "V_L": V_N, "vessel_L": ves_N,
         "fill_pct": 100.0 * V_N / ves_N if ves_N else None, "duration_d": run_N, "peak_vcd": peak_N,
         "OUR": peak_N * st.qO2_pmol_cell_h,
         "kLa_air": peak_N * st.qO2_pmol_cell_h / dC if dC > 0 else None}

    # timeline from vial thaw
    s3 = stages[0]
    cells = s3["x0"] * s3["V_L"]
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
        return None if f is None else _chk(s["min_fill"] * 100 - 1e-6 <= f <= s["max_fill"] * 100 + 1e-6,
                                           f"{f:.0f} %")

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
            "basis_label": basis_label(st), "threshold_pct": st.growth_threshold_pct,
            "stages": stages, "n": N, "timeline": tl, "checks": checks, "ok": not fails,
            "overall": overall, "notes": notes, "demand": dem}


def apply_seed_train_to_n1(p):
    """seedtrain.apply_to_n1: run the twin's N-1 as the batch N-1 stage of a seed train sized for the
    twin N reactor. Idempotent (the train does not depend on the n1 parameters). Returns the train or None."""
    st = getattr(p, "seedtrain", None)
    if st is None or not st.apply_to_n1:
        return None
    tr = seed_train(p, n_basis="twin")
    s = tr["stages"][-1]
    n1 = p.n1
    n1.enabled = True
    n1.V_L = n1.transfer_L = float(s["V_L"])
    n1.VVD = 0.0                                           # batch stage: no permeate to waste
    n1.X0 = float(s["x0"])
    n1.X_max = max(float(s["xt"]), n1.X0)
    days = max(float(s["duration_d"]), 1.0 / 24.0)
    n1.days, n1.target_seed_vcd = days, 0.0                # fixed duration from the chosen basis
    n1.max_days = max(n1.max_days, days)
    n1.mu_h = math.log(n1.X_max / n1.X0) / (days * 24.0) if n1.X_max > n1.X0 else 0.0
    return tr


def bio_design(p):
    """Everything the page shows: Bio_sim curves and results, Seed_sim curves, seed train and demand."""
    b, st = p.bio, p.seedtrain
    run = bio_run(p)
    rows, R = run["rows"], run["results"]
    step = max(1, int(round(b.report_days / max(b.dt_d, 1e-9))))
    tr = dict(seed_train(p, "twin" if st.apply_to_n1 else None, run))
    dem = tr.pop("demand")
    seeds = {}
    for s in tr["stages"]:
        sim = seed_sim(b, s["x0"], st.sim_days)
        seeds[s["key"]] = {"name": s["name"], "x0": s["x0"], "xt": s["xt"], "transfer_day": s["dur_bio"],
                           "peak_net": s["peak_net"], "threshold_net": s["threshold_net"],
                           "rows": [{k: r[k] for k in SEED_KEYS} for r in sim]}
    notes = [f"Seed train basis: {tr['basis_label']}. N from: {tr['n']['source']}.",
             f"Bio_sim is a fed-batch model (CSFR feeding, temperature shift on day {b.shift_day:g}); it sizes N on "
             "the Demand basis. The twin's IDB N reactor is a perfusion reactor: choose 'Biology model' as the "
             "titre curve to use its shape (scaled by the case or titre_scale).",
             "Rules of thumb used whatever the basis: seed densities, the growth threshold, qO2, purification "
             "yield, batch success, turnaround and the check limits.",
             "Transfer-point values are interpolated at the transfer time (the workbook shows the next 0.1-day row)."]
    if p.upstream.molecule == "biology":
        link = "perfusion qP·Xv/VVD" if b.titre_link == "perfusion" else "fed-batch titre"
        notes.append(f"N titre curve in this run: biology model, {link}, × titre_scale {p.upstream.titre_scale:.3g}.")
    if st.apply_to_n1:
        notes.append("The twin's N-1 runs as the batch N-1 stage of this train (VVD 0; volume, duration and seed "
                     "density from the train); the train is sized for the twin N reactor.")
    notes += tr["notes"]
    return _rnd({"bio": {"results": R, "series": rows, "curve": rows[::step], "target_gL": b.target_titre_gL,
                         "run_days": b.run_days, "dt_d": b.dt_d},
                 "seed_sim": seeds, "seed_train": tr, "demand": dem, "notes": notes})


if __name__ == "__main__":
    from icb_process import ICBParams
    d = bio_design(ICBParams())
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
