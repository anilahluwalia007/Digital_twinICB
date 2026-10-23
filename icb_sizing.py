"""
icb_sizing.py — equipment sizing for the IDB mAb line
NavNevAI_mAbWasteWater · CarbonTesseract 2026

size_equipment(p, res, history, sol_day)  sizes every unit from one hourly run (called by the simulator)
size_portfolio(p, seed)                   worst case across all titre curves incl. the biology model (/api/icb/sizing)

  Bioreactors   working volume is an input (sets the permeate flow): nominal size, media bags, pumps
  Seed train    N-3, N-2, N-1 and N vessels from icb_bio.seed_train (biology and/or rules of thumb),
                each in its culture mode (batch, fed-batch or perfusion: culture.preset)
  TFF F-101     membrane area from the peak permeate flow at a design flux
  Protein A     resin per train from the peak g/day and the cycle limit (load + recovery time)
  Vessels       peak working volume × margin / fill, next standard size
  Depth filter  area from the largest batch (L/m²) and the time available
  AEX, CEX      resin scaled so columns and load pump stay under the busy-time target
  Virus filter  flow from the busy-time target, area from flux, change volume from L/m²
  UF/DF         area from permeate flux in the time available and g/m² per lot
  Filters/pumps from the sized flows;  totes from the peak day per solution
  WFI, effluent generation / discharge for the peak day, storage for the hourly peaks

Busy-time scaling: with flows in CV/h (or flux × area), load time and the number of cycles both
scale with 1/size, so required size = configured size × busy share / target share.
Changing a column, filter or membrane moves the vessel peaks: apply and re-run until they settle.
"""
from __future__ import annotations
import copy
import math
from icb_process import (WASTE_STREAMS, BUFFER_OPS, SOLUTIONS, TITRE_PROFILES, capture_cycle, det_frac,
                         profile_names)
from icb_bio import stage_modes

VESSEL_SIZES = (1, 2, 5, 10, 20, 50, 100, 200, 300, 500, 1000, 1500, 2000, 3000, 5000, 10000, 20000, 50000)
TOTE_SIZES = (10, 20, 50, 100, 200, 500, 1000)
VESSEL_TAGS = {"sump": "T-101", "pav": "T-102", "vip": "T-VI", "dfs": "T-103", "aexv": "T-201",
               "cexv": "T-202", "vrfv": "T-301", "ufdf": "T-302", "bdsm": "M-301"}
ROUTE_TANKS = {"direct": ("T-EQ", "Direct drain / EQ tank"), "neutralisation": ("T-NEU", "Neutralisation tank"),
               "bio-kill": ("T-BK", "Bio-kill tank"), "all": ("WWTP", "All wastewater to the WWTP")}


# ── helpers ──────────────────────────────────────────────────────────
def _n(v):
    return float(v) if isinstance(v, (int, float)) and not isinstance(v, bool) and math.isfinite(v) else 0.0


def _v(x):
    if not isinstance(x, (int, float)) or isinstance(x, bool):
        return x
    if not math.isfinite(x):
        return None
    a = abs(x)
    return round(x, 3 if a < 1 else 2 if a < 100 else 1)


def _ceil(x, step):
    return round(math.ceil(x / step - 1e-9) * step, 6)


def _std(x, sizes=VESSEL_SIZES):
    for s in sizes:
        if s >= x - 1e-9:
            return s
    return _ceil(x, sizes[-1])


def _tote(L):
    if L <= 0:
        return None
    for s in TOTE_SIZES:
        if s >= L - 1e-9:
            return f"1 × {s} L"
    return f"{math.ceil(L / TOTE_SIZES[-1] - 1e-9)} × {TOTE_SIZES[-1]} L"


def _id_cm(CV_L, bed_cm):
    return math.sqrt(4.0 * CV_L * 1000.0 / (math.pi * max(bed_cm, 1e-9)))


def _row(tag, item, required, unit, recommended=None, configured=None, basis="", flag=""):
    return {"tag": tag, "item": item, "required": _v(required), "unit": unit,
            "recommended": _v(recommended), "configured": _v(configured), "basis": basis, "flag": flag}


def _pump(rows, sz, tag, item, flow, basis):
    rows.append(_row(tag, item, flow, "L/h", _ceil(flow * sz.pump_margin, 0.1), None,
                     f"{basis}; × {sz.pump_margin:g} margin"))


def _by_day(flows):
    d = {}
    for i, v in enumerate(flows):                 # history rows are hourly from h = 0
        d[i // 24] = d.get(i // 24, 0.0) + v
    return d


def _storage(flows, q):
    """Smallest tank that absorbs an hourly profile while it is emptied (or refilled) at q L/h."""
    s = peak = 0.0
    for v in flows:
        s = max(0.0, s + v - q)
        peak = max(peak, s)
    return peak


def _busy(steps, prefixes, scaling):
    """Busy hours per campaign day for lanes starting with prefixes, split into (scaling, fixed)."""
    S, F = {}, {}
    for s in steps:
        if not s["lane"].startswith(prefixes):
            continue
        t0, t1 = float(s["t0"]), float(s["t1"])
        tgt = S if scaling(s) else F
        while t1 - t0 > 1e-9:
            d = int(t0 // 24)
            e = min(t1, (d + 1) * 24.0)
            tgt[d] = tgt.get(d, 0.0) + e - t0
            t0 = e
    return S, F


def _scale(S, F, days, lanes, util):
    """Size factor k so that S/k + F <= util × 24 h × lanes on every production day; also peak busy share."""
    T = util * 24.0 * lanes
    k = u = 0.0
    for d in days:
        s, f = S.get(d, 0.0), F.get(d, 0.0)
        u = max(u, (s + f) / (24.0 * lanes))
        if s > 1e-9:
            k = max(k, s / (T - f) if T - f > 1e-6 else math.inf)
    return (k if k > 0 else None), u


def _size_columns(rows, sz, steps, prod, lane, tag, label, c, bed_cm, add_frac, pump_tag, dose_tag, dose_item):
    n = max(1, int(c.columns))
    S, _ = _busy(steps, (lane,), lambda s: True)
    SL, _ = _busy(steps, (lane,), lambda s: s["kind"] == "load")
    k1, u = _scale(S, {}, prod, n, sz.util_target)
    k2, uL = _scale(SL, {}, prod, 1, sz.util_target)
    ks = [k for k in (k1, k2) if k]
    if not ks:
        rows.append(_row(tag, f"{label} resin per column", None, "L", None, c.CV_L,
                         "did not run in the production days of this run: kept as configured"))
        return None
    CV_req = c.CV_L * max(ks)
    CV_rec = _ceil(max(CV_req, c.CV_L * sz.min_scale), 0.1)
    rows.append(_row(tag, f"{label} resin per column ({n} columns)", CV_req, "L", CV_rec, c.CV_L,
                     f"at the configured size the columns are busy {u:.0%} and the load pump {uL:.0%} of the "
                     f"peak day; sized so neither exceeds {sz.util_target:.0%}"))
    rows.append(_row(tag, f"{label} column i.d. at {bed_cm:g} cm bed", _id_cm(CV_rec, bed_cm), "cm",
                     _ceil(_id_cm(CV_rec, bed_cm), 1), _id_cm(c.CV_L, bed_cm), "per column"))
    u_lin = c.flow_CVh * bed_cm
    rows.append(_row(tag, f"{label} linear velocity", u_lin, "cm/h", None, None, f"{c.flow_CVh:g} CV/h × {bed_cm:g} cm",
                     "above max_u_cmh: lower the bed height or flow_CVh" if u_lin > sz.max_u_cmh else ""))
    q = c.flow_CVh * CV_rec
    _pump(rows, sz, pump_tag, f"{label} load / buffer pump", q, f"{c.flow_CVh:g} CV/h × {CV_rec:g} L")
    if add_frac > 0:
        _pump(rows, sz, dose_tag, dose_item, q * add_frac / (1 + add_frac),
              f"{100 * add_frac / (1 + add_frac):.1f} % of the load flow")
    return CV_rec


# ── one run ──────────────────────────────────────────────────────────
def size_equipment(p, res, history, sol_day=None):
    sz, cp, po = p.sizing, p.capture, p.polish
    steps = res.get("steps") or []
    kpi = res.get("kpi") or {}
    util, ms, am = sz.util_target, sz.min_scale, sz.area_margin
    prod = sorted({int(r["h"] // 24) for r in history if r.get("usp_on")})
    prod_set = set(prod)
    groups, sug, notes = [], {}, []

    def group(name):
        rows = []
        groups.append({"name": name, "rows": rows})
        return rows

    def suggest(sec, k, v):
        sug.setdefault(sec, {})[k] = _v(v)

    # upstream
    rows = group("Upstream: bioreactors, media and cell removal")
    n1, br = p.n1, p.bioreactor
    n1s = (res.get("design") or {}).get("n1") or {}
    if n1.enabled:
        batch_n1 = n1.VVD <= 0                                  # batch or fed-batch N-1 from the seed train
        kind = ("fed-batch" if batch_n1 and p.seedtrain.apply_to_n1 and stage_modes(p)["n1"] == "fedbatch"
                else "batch")
        rows.append(_row("R-101", (f"N-1 {kind} seed (seed-train stage)" if batch_n1 else "N-1 perfusion seed")
                         + ", working volume", n1.V_L, "L", _std(n1.V_L / sz.vessel_fill),
                         n1.V_L, f"nominal = working / {sz.vessel_fill:g} fill; seeds N at {n1s.get('n_seed_vcd')} "
                         f"×10⁶/mL after {n1s.get('days')} d ({n1s.get('reason')})"))
        m1 = n1.VVD * n1.V_L if n1.VVD > 0 else n1.V_L          # batch / fed-batch N-1: one fill per day of supply
        rows.append(_row("BAG-101", "N-1 media hold", m1 * sz.media_days, "L",
                         _std(m1 * sz.media_days * sz.vessel_margin), None,
                         (f"{kind}: {'fill plus feed' if kind == 'fed-batch' else 'one fill'} of {n1.V_L:.1f} L"
                          if batch_n1 else f"{n1.VVD:g} VVD × {n1.V_L:g} L")
                         + f" × {sz.media_days:g} d × {sz.vessel_margin:g}"))
        _pump(rows, sz, "P-100", "N-1 media fill pump" if batch_n1 else "N-1 media / permeate pump", m1 / 24.0,
              f"one fill of {n1.V_L:.1f} L in 24 h" if batch_n1 else f"{n1.VVD:g} VVD × {n1.V_L:g} L / 24 h")
    Fin = max((_n(r.get("Fin_Lh")) for r in history), default=0.0)
    Fh = max((_n(r.get("Fh_Lh")) for r in history), default=0.0)
    rows.append(_row("R-102", "N reactor, working volume", br.V_L, "L", _std(br.V_L / sz.vessel_fill), br.V_L,
                     "design input: it sets the permeate flow; nominal = working / fill"))
    rows.append(_row("BAG-102", "N media hold", Fin * 24 * sz.media_days, "L",
                     _std(Fin * 24 * sz.media_days * sz.vessel_margin), None,
                     f"peak feed {Fin:.2f} L/h × 24 h × {sz.media_days:g} d × {sz.vessel_margin:g}"))
    _pump(rows, sz, "P-101", "N outflow pump to TFF", Fin, "peak N outflow")
    a = Fh / max(sz.tff_flux_LMH, 1e-9)
    rows.append(_row("F-101", "TFF cell-removal membrane area", a, "m²", _ceil(a * am, 0.01), None,
                     f"peak permeate {Fh:.2f} L/h at {sz.tff_flux_LMH:g} LMH × {am:g}"))

    # seed train N-3 → N (icb_bio.py: biology and/or rules of thumb)
    tr = (res.get("bio") or {}).get("seed_train") or {}
    if tr.get("stages"):
        rows = group("Seed train N-3 → N (biology / rules of thumb)")
        cols = tr["stages"] + [tr["n"]]
        bad = {i: [] for i in range(len(cols))}
        for c in tr.get("checks", []):
            for i, cell in enumerate(c["cells"]):
                if cell and not cell["ok"]:
                    bad[i].append(f"{c['label']}: {cell['text']}")
        for i, s in enumerate(cols):
            mode = {"batch": "batch", "fedbatch": "fed-batch", "perfusion": "perfusion"}.get(s.get("mode"), "")
            basis = (f"{mode} stage; {tr['basis_label']}; seed {s['x0']:g} → transfer {s['xt']:.2f} ×10⁶/mL in "
                     f"{s['duration_d']:.2f} d; max fill {s['max_fill']:g}; medium {_n(s.get('media_L')):.1f} L"
                     if s["key"] != "n"
                     else f"{mode} N from {s.get('source', '')}; seeded at {s['x0']:g} ×10⁶/mL")
            rows.append(_row(s["tag"], f"{s['name']} {str(s.get('type', '')).lower()}, working volume",
                             s["V_L"], "L", s["vessel_L"], None, basis, "; ".join(bad[i])))

    # Protein A
    rows = group("Protein A capture (Brij-9 in-line)")
    cyc = capture_cycle(p)
    if cp.mode == "pcc":
        rows.append(_row("C-101", "PCC column, CV each", None, "L", None, p.pcc.CV_L, "PCC mode: kept as configured"))
        notes.append("PCC capture is reported as configured; switch Capture to Case 5 Protein A to size the columns.")
    else:
        gday = {d: g for d, g in _by_day([_n(r.get("harvest_g_h")) for r in history]).items() if d in prod_set}
        dpk = max(gday, key=gday.get) if gday else None
        g_pk = gday.get(dpk, 0.0) if dpk is not None else 0.0
        L_pk = sum(_n(r.get("Fh_Lh")) for r in history if int(r["h"] // 24) == dpk) if dpk is not None else 0.0
        C = g_pk / L_pk if L_pk > 0 else 0.0
        if C > 0:
            dose = det_frac(p)
            CV_pa = cyc["CV_L"]
            load_h = cp.DBC_gL * cp.load_util * (1 + dose) / (C * max(cp.load_flow_CVh, 1e-9))
            lim_h = max(load_h, (load_h + cyc["recovery_h"]) / cyc["trains"])
            max_cyc = 24.0 / lim_h
            CV_req = g_pk / (max_cyc * util * cp.load_util * cp.DBC_gL)
            CV_rec = _ceil(max(CV_req, CV_pa * ms), 0.1)
            rows.append(_row("C-101…104", f"Resin per train ({cyc['trains']} trains)", CV_req, "L", CV_rec, CV_pa,
                             f"peak {g_pk:.0f} g/day on campaign day {dpk + 1} at {C:.2f} g/L; a train loads "
                             f"{load_h:.1f} h and recovers {cyc['recovery_h']:.1f} h → at most {max_cyc:.1f} cycles/day; "
                             f"sized for ≤ {util:.0%} of that"))
            dia = _id_cm(CV_rec, sz.pa_bed_cm)
            rows.append(_row("C-101…104", f"Train i.d. at {sz.pa_bed_cm:g} cm bed", dia, "cm", _ceil(dia, 1),
                             _id_cm(CV_pa, sz.pa_bed_cm),
                             "if a train is two columns in series, keep this i.d. and split the bed height"))
            u_lin = cp.load_flow_CVh * sz.pa_bed_cm
            rows.append(_row("C-101…104", "Load linear velocity", u_lin, "cm/h", None, None,
                             f"{cp.load_flow_CVh:g} CV/h × {sz.pa_bed_cm:g} cm",
                             "above max_u_cmh: lower load_flow_CVh or the bed height" if u_lin > sz.max_u_cmh else ""))
            q_load = cp.load_flow_CVh * CV_rec
            _pump(rows, sz, "P-102", "Protein A load pump (permeate + Brij-9)", q_load,
                  f"{cp.load_flow_CVh:g} CV/h × {CV_rec:g} L")
            _pump(rows, sz, "P-102B", "Protein A buffer / recovery pump", cp.recovery_flow_CVh * CV_rec,
                  f"{cp.recovery_flow_CVh:g} CV/h × {CV_rec:g} L")
            if dose > 0:
                _pump(rows, sz, "P-BRIJ", "Brij-9 dosing pump", q_load * dose / (1 + dose),
                      f"{100 * dose / (1 + dose):.1f} % of the load flow")
            ncyc = kpi.get("capture_cycles") or 0
            per = ncyc / max(cyc["trains"], 1)
            rows.append(_row("C-101…104", "Cycles per train this campaign", per, "cycles", None, cp.resin_cycles_max,
                             f"{ncyc} cycles at the configured size; configured = resin cycle limit",
                             "resin limit reached" if per > cp.resin_cycles_max else ""))
            suggest("capture", "column_capacity_g", CV_rec * cp.DBC_gL)
        else:
            notes.append("Protein A not sized: no product reached the surge tank in this run.")

    # vessels
    rows = group("Vessels and surge tanks (peak working volume in this run)")
    overflow = _n((kpi.get("losses_g") or {}).get("sump_overflow")) > 0
    caps = {"sump": cp.sump_max_L if cp.mode != "pcc" else p.surge.V_max_L, "pav": p.vi.eluate_surge_max_L}
    for v in res.get("vessels") or []:
        need = _n(v.get("peak_L"))
        if need <= 0:
            continue
        nom = _std(need * sz.vessel_margin / sz.vessel_fill)
        cfg = caps.get(v["key"])
        flag = ""
        if v["key"] == "sump" and overflow:
            flag = "overflowed in this run: the real need is larger; apply and re-run"
        elif cfg is not None and need > cfg + 1e-6:
            flag = "peak above the configured capacity"
        rows.append(_row(VESSEL_TAGS.get(v["key"], ""), v["label"], need, "L", nom, cfg,
                         f"peak × {sz.vessel_margin:g} / {sz.vessel_fill:g} fill, next standard size"
                         + ("" if cfg is not None else "; no capacity limit in the model"), flag))
        if v["key"] == "sump" and cp.mode != "pcc":
            suggest("capture", "sump_max_L", nom * sz.vessel_fill)
        if v["key"] == "pav":
            suggest("vi", "eluate_surge_max_L", nom * sz.vessel_fill)

    # depth filtration
    rows = group("Depth filtration")
    d = p.depth
    loads = [s["L"] for s in steps if s["lane"] == "depth" and s["kind"] == "load"]
    if loads:
        A_cap = max(loads) / max(d.capacity_Lm2, 1e-9)
        S, F = _busy(steps, ("depth",), lambda s: s["kind"] == "load")
        k, u = _scale(S, F, prod, 1, util)
        flag = "flushes alone exceed the time available" if k == math.inf else ""
        A_time = d.area_m2 * k if k and math.isfinite(k) else 0.0
        A_req = max(A_cap, A_time)
        A_rec = _ceil(max(A_req * am, d.area_m2 * ms), 0.01)
        rows.append(_row("FH-101", "Depth filter area (per filter)", A_req, "m²", A_rec, d.area_m2,
                         f"largest batch {max(loads):.1f} L at {d.capacity_Lm2:g} L/m² → {A_cap:.2f} m²; "
                         f"busy {u:.0%} of the peak day at the configured area; × {am:g}", flag))
        _pump(rows, sz, "P-103", "Depth filtration feed pump", d.flux_LMH * A_rec, f"{d.flux_LMH:g} LMH × {A_rec:g} m²")
        _pump(rows, sz, "P-103F", "Depth filter flush pump", d.flush_flux_LMH * A_rec,
              f"{d.flush_flux_LMH:g} LMH × {A_rec:g} m²")
        rows.append(_row("FH-101", "Depth filters used this campaign", kpi.get("depth_filters"), "filters",
                         None, None, "one new filter per batch"))
        suggest("depth", "area_m2", A_rec)

    # polishing
    rows = group("AEX chromatography (flow-through)")
    cv = _size_columns(rows, sz, steps, prod, "aex", "C-201/202", "AEX", p.aex, sz.aex_bed_cm,
                       p.aex.base_frac + p.aex.nacl_frac, "P-104", "P-BASE", "Tris base + NaCl dosing pump")
    if cv:
        suggest("aex", "CV_L", cv)
    rows = group("CEX chromatography (bind and elute)")
    cv = _size_columns(rows, sz, steps, prod, "cex", "C-203/204", "CEX", p.cex, sz.cex_bed_cm,
                       p.cex.acid_frac, "P-201", "P-ACID", "Acetic acid dosing pump")
    if cv:
        suggest("cex", "CV_L", cv)
    if po.order == "cex_aex":
        notes.append("polish.order = cex_aex: the pump tags follow the AEX → CEX drawing.")

    # virus filtration
    rows = group("Virus retentive filtration")
    v = p.vrf
    S, _ = _busy(steps, ("vrf",), lambda s: True)
    k, u = _scale(S, {}, prod, 1, util)
    vol = sum(_n(r.get("L_vf_in")) for r in history)
    if k:
        q_req = v.flow_Lh * k
        q_rec = _ceil(max(q_req, v.flow_Lh * ms), 0.1)
        A_rec = _ceil(q_rec / max(sz.vf_flux_LMH, 1e-9) * am, 0.01)
        cap = A_rec * sz.vf_Lm2
        rows.append(_row("F-301 A/B", "Filtration flow", q_req, "L/h", q_rec, v.flow_Lh,
                         f"busy {u:.0%} of the peak day at {v.flow_Lh:g} L/h; sized for ≤ {util:.0%}"))
        rows.append(_row("F-301 A/B", "Virus filter area (per holder)", q_rec / sz.vf_flux_LMH, "m²", A_rec, None,
                         f"{q_rec:g} L/h at {sz.vf_flux_LMH:g} LMH × {am:g}"))
        rows.append(_row("F-301 A/B", "Volume per filter before change", cap, "L", None, po.vf_capacity_L,
                         f"{A_rec:g} m² × {sz.vf_Lm2:g} L/m²"))
        rows.append(_row("F-301 A/B", "Filters this campaign", math.ceil(vol / cap) if cap > 0 else None, "filters",
                         None, None, f"{vol:.0f} L filtered"))
        _pump(rows, sz, "P-202", "Virus filtration feed pump", q_rec, "sized filtration flow")
        suggest("vrf", "flow_Lh", q_rec)
        suggest("polish", "vf_capacity_L", cap)

    # UF/DF
    rows = group("UF/DF and formulation")
    s = p.sptff
    S, F = _busy(steps, ("ufdf",), lambda st: st["kind"] in ("perm", "df"))
    k, u = _scale(S, F, prod, 1, util)
    lots = res.get("lots") or []
    g_lot = max((_n(L.get("g")) for L in lots), default=0.0) / max(s.sptff_yield, 1e-9)
    if k:
        fin = math.isfinite(k)
        flux_req = s.perm_flux_Lh * (k if fin else 1.0)
        A_flux = flux_req / max(sz.ufdf_flux_LMH, 1e-9)
        A_load = g_lot / max(sz.ufdf_load_gm2, 1e-9)
        A_req = max(A_flux, A_load)
        A_cfg = s.perm_flux_Lh / max(sz.ufdf_flux_LMH, 1e-9)
        A_rec = _ceil(max(A_req * am, A_cfg * ms), 0.01)
        rows.append(_row("UF-301", "UF/DF membrane area", A_req, "m²", A_rec, A_cfg,
                         f"permeate {flux_req:.1f} L/h at {sz.ufdf_flux_LMH:g} LMH → {A_flux:.2f} m²; largest lot "
                         f"{g_lot:.0f} g at {sz.ufdf_load_gm2:g} g/m² → {A_load:.2f} m²; busy {u:.0%} of the peak day",
                         "" if fin else "flushes, CIP and transfers alone exceed the time available"))
        _pump(rows, sz, "P-UF", "UF/DF feed (crossflow) pump", A_rec * sz.ufdf_crossflow_LMH,
              f"{sz.ufdf_crossflow_LMH:g} LMH × {A_rec:g} m²")
        suggest("sptff", "perm_flux_Lh", A_rec * sz.ufdf_flux_LMH)
    _pump(rows, sz, "P-301", "Transfer VRF vessel → UF/DF", s.transfer_Lh, "transfer_Lh")
    _pump(rows, sz, "P-302", "Transfer UF/DF → BDS hold mixer", s.transfer_Lh, "transfer_Lh")

    # BDS and cryovault
    rows = group("BDS hold, fill and cryovault")
    b = p.bds
    _pump(rows, sz, "P-303", "BDS fill pump", b.fill_Lh, "fill_Lh")
    per_lot = max((L.get("containers", 0) for L in lots), default=0)
    rows.append(_row("CRYO", f"Cryovault containers ({b.container_L:g} L)", kpi.get("containers"), "containers",
                     None, None, f"{len(lots)} lots, up to {per_lot} containers per lot: freezer positions for the campaign"))

    # in-line filters
    rows = group("In-line 0.2 µm filters")
    CVpa = sug.get("capture", {}).get("column_capacity_g", cp.column_capacity_g) / max(cp.DBC_gL, 1e-9)
    for tag, item, q in (("FL-101", "Protein A eluate line", cp.recovery_flow_CVh * CVpa),
                         ("FL-201", "AEX outlet", p.aex.flow_CVh * sug.get("aex", {}).get("CV_L", p.aex.CV_L)),
                         ("FL-202", "CEX outlet", p.cex.flow_CVh * sug.get("cex", {}).get("CV_L", p.cex.CV_L)),
                         ("FL-301", "Virus filtrate line", sug.get("vrf", {}).get("flow_Lh", p.vrf.flow_Lh)),
                         ("FL-BDS", "BDS fill filter", b.fill_Lh)):
        a = q / max(sz.filter_flux_LMH, 1e-9)
        rows.append(_row(tag, item, a, "m²", _ceil(a * am, 0.01), None,
                         f"{q:.1f} L/h at {sz.filter_flux_LMH:g} LMH × {am:g}"))

    # solutions
    rows = group("Solutions: totes and buffer preparation (peak day of use)")
    peak, tot = {}, {}
    for day, use in (sol_day or {}).items():
        for k2, L in use.items():
            if k2 == "wfi":
                continue
            if L > peak.get(k2, (0.0, 0))[0]:
                peak[k2] = (L, day)
            tot[day] = tot.get(day, 0.0) + L
    for k2, (L, day) in sorted(peak.items(), key=lambda kv: -kv[1][0]):
        rows.append(_row("TOTE", SOLUTIONS.get(k2, {}).get("label", k2), L, "L/day",
                         _tote(L * sz.tote_days * sz.vessel_margin), None,
                         f"peak on campaign day {day + 1}; {sz.tote_days:g} d supply × {sz.vessel_margin:g}"))
    if tot:
        top = max(tot.values())
        rows.append(_row("BUF-PREP", "Buffer preparation capacity (all solutions)", top, "L/day",
                         _ceil(top * sz.vessel_margin, 1), None, f"peak day × {sz.vessel_margin:g}"))

    # WFI and wastewater
    rows = group("Utilities: WFI and wastewater")
    vm, vf, pm = sz.vessel_margin, sz.vessel_fill, sz.pump_margin
    wfi = [sum(_n(r.get("wfi_" + k3)) for k3 in BUFFER_OPS) for r in history]
    if any(wfi):
        q = max(_by_day(wfi).values()) / 24.0
        st = _storage(wfi, q)
        rows.append(_row("WFI", "WFI generation", q, "L/h", _ceil(q * pm, 0.1), None, "peak day / 24 h"))
        rows.append(_row("WFI", "WFI storage", st, "L", _std(st * vm / vf), None,
                         "covers the hourly draw peaks at that generation rate"))
        rows.append(_row("WFI", "WFI loop, peak hourly draw", max(wfi), "L/h", _ceil(max(wfi) * pm, 0.1), None, ""))
    routes = {}
    for k3, (_, route, _, _) in WASTE_STREAMS.items():
        routes.setdefault(route, []).append(k3)
    routes["all"] = list(WASTE_STREAMS)
    for route, keys in routes.items():
        flow = [sum(_n(r.get("ww_" + k3)) for k3 in keys) for r in history]
        if not any(flow):
            continue
        pk_day = max(_by_day(flow).values())
        q = pk_day / 24.0
        st = _storage(flow, q)
        basis = f"discharged at {q:.1f} L/h (peak day {pk_day:.0f} L / 24 h)"
        if route == "bio-kill":
            st, basis = max(st, pk_day), "one day of collection for batch kill"
        tag, label = ROUTE_TANKS.get(route, ("WW", route))
        rows.append(_row(tag, label + ": storage", st, "L", _std(st * vm / vf), None, basis))
        _pump(rows, sz, tag, label + ": discharge pump", q, "peak day / 24 h")
        rows.append(_row(tag, label + ": peak inflow (drain capacity)", max(flow), "L/h",
                         _ceil(max(flow) * pm, 0.1), None, "largest hourly flow"))

    notes.append(f"Sized from this run's hourly profiles: units busy ≤ {util:.0%} of the peak day, vessels × "
                 f"{vm:g} / {vf:g} fill, pumps × {pm:g}, filters and membranes × {am:g}. Bigger columns, filters or "
                 "membranes move the vessel peaks: press 'Apply sized values and re-run' and check again (one or two rounds).")
    return {"groups": groups, "suggested_params": sug, "notes": notes,
            "scope": (res.get("design") or {}).get("molecule_label")}


# ── worst case across molecules ──────────────────────────────────────
def size_portfolio(p, seed=0):
    from icb_dynamic import ICBSimulator          # here to avoid a circular import
    names = profile_names(p)
    best, order, gnames, sug, notes, labels = {}, [], [], {}, [], []
    for name in names:
        q = copy.deepcopy(p)
        q.upstream.mode, q.upstream.molecule = "profile", name
        res = ICBSimulator(q, forecaster=None, seed=seed).run()
        lab = (res.get("design") or {}).get("molecule_label", name)
        labels.append(lab)
        sz = res.get("sizing") or {}
        if sz.get("error"):
            notes.append(f"{lab}: {sz['error']}")
            continue
        for g in sz.get("groups", []):
            if g["name"] not in gnames:
                gnames.append(g["name"])
            for r in g["rows"]:
                key = (g["name"], r["tag"], r["item"])
                cur = best.get(key)
                if cur is None:
                    order.append(key)
                if cur is None or _n(r["required"]) > _n(cur["required"]):
                    best[key] = dict(r, set_by=lab)
        for sec, kv in (sz.get("suggested_params") or {}).items():
            for k, v in kv.items():
                if _n(v) > sug.setdefault(sec, {}).get(k, -math.inf):
                    sug[sec][k] = v
        for n in sz.get("notes", []):
            if n not in notes:
                notes.append(n)
    notes.insert(0, "Worst case of " + ", ".join(labels) + " with the same settings; 'Set by' names the molecule "
                    "that needs the most.")
    groups = [{"name": g, "rows": [best[k] for k in order if k[0] == g]} for g in gnames]
    return {"groups": groups, "suggested_params": sug, "notes": notes, "molecules": labels}
