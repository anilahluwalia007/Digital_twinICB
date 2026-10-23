"""
icb_cases.py — basis-of-design Cases 1–5 for the IDB mAb digital twin
NavNevAI_mAbWasteWater · CarbonTesseract 2026

  #  Case                              Avg titre   Batches/yr  Yield  Capture load
                                       g/L/day                        g/L resin
  1  Design Case: High Yield/Volume      3.0          20       71 %       80
  2  Design Case: Turndown Volume        3.0           3        –          –
  3  Design Case: Low Yield/Volume       0.4           3       25 %       35
  4  Rating Case: Rilveg                 0.6          20       71 %        –
  5  Scenario Case: High-High Yield      4.8          20       71 %        –

Cases 1–3 bracket the operating range every item of equipment must handle; the design case is their
envelope (worst case per item). Case 4 rates the lead molecule Rilveg for cost of goods. Case 5 is a
higher-titre scenario: simulate it, scope the incremental capital against Case 1, and SteerCo decides
whether it replaces Case 1 as the high design case.

How a case changes a run (applied on top of the Inputs, at run time):
  titre    the chosen titre-curve shape is scaled (upstream.titre_scale) so its run average equals the
           case value. Basis "productivity" (default): g made per L of N working volume per day
           = mean permeate titre × VVD. Basis "titre": mean permeate titre in g/L.
  shape    "auto" picks the Figure 2 curve whose own average is nearest the target (least scaling);
           "biology" uses the Bio_sim curve from icb_bio.py (fed-batch titre or perfusion qP·Xv/VVD).
  yield    harvest-to-DS train yield (Protein A → fill): every step yield is scaled by one factor.
  load     Protein A load in g per L resin: DBC_gL = load / load_util; the column volume
           (column_capacity_g / DBC_gL from the Inputs) is kept, so the grams per cycle follow.
  batches  per-year totals = one simulated batch × batches per year; reactor and downstream lines
           from occupancy (days per batch + turnaround) against the operating days per year.
Blank cells: Cases 2 and 5 inherit yield and load from Case 1 (the same process at another volume or
titre); Case 4 keeps the load of its recipe set. Every value can be edited on the Cases tab.
  recipe   Case 4 first loads the Rilveg 2936 recipe set (icb_rilveg.py: reactor, TFF, LambdaFabSelect
           capture, polishing and UF/DF solutions); the titre, yield and load locks then act on it.
  titrants every case sizes the in-line Tris base before AEX and acetic acid before CEX from the pH
           each step has to reach (icb_process.titrant_plan); confirm both by titration.
Case 4 yield 71 % = the three yields the map states (0.80 × 0.97 × 0.92 = 0.713); lower it once the
depth-filter, virus-filter and UF/DF losses are measured.
"""
from __future__ import annotations
import copy
import math
from icb_process import (TITRE_PROFILES, EXTRA_PROFILES, profile_points, profile_label, titre_at, train_yield,
                         capture_cycle, process_design, titrant_plan)

VALUE_KEYS = ("avg_titre_gLd", "batches_per_year", "yield_est", "capture_load_gL")
CASES = {
    "1": {"name": "Design Case: High Yield/Volume", "short": "High Yield/Volume", "kind": "design",
          "avg_titre_gLd": 3.0, "batches_per_year": 20.0, "yield_est": 0.71, "capture_load_gL": 80.0},
    "2": {"name": "Design Case: Turndown Volume", "short": "Turndown Volume", "kind": "design",
          "avg_titre_gLd": 3.0, "batches_per_year": 3.0, "yield_est": None, "capture_load_gL": None},
    "3": {"name": "Design Case: Low Yield/Volume", "short": "Low Yield/Volume", "kind": "design",
          "avg_titre_gLd": 0.4, "batches_per_year": 3.0, "yield_est": 0.25, "capture_load_gL": 35.0},
    "4": {"name": "Rating Case: Rilveg", "short": "Rilveg", "kind": "rating",
          "avg_titre_gLd": 0.6, "batches_per_year": 20.0, "yield_est": 0.71, "capture_load_gL": None},
    "5": {"name": "Scenario Case: High-High Yield", "short": "High-High Yield", "kind": "scenario",
          "avg_titre_gLd": 4.8, "batches_per_year": 20.0, "yield_est": 0.71, "capture_load_gL": None},
}
INHERIT = {"2": "1", "5": "1"}            # blank yield / load taken from this case
RECIPES = {"4": "rilveg"}                 # case id -> molecule recipe set
DESIGN_CASES = ("1", "2", "3")
KIND_LABEL = {"design": "Design case", "rating": "Rating case", "scenario": "Scenario case"}
PURPOSE = {
    "1": "Design case: high titre at high volume. With Cases 2 and 3 it brackets the range every item must handle.",
    "2": "Design case: same titre as Case 1 at turndown volume (3 batches a year).",
    "3": "Design case: low titre, low yield and low volume; the low end of the operating range.",
    "4": "Rating case: the lead pipeline molecule Rilveg, simulated for detailed cost of goods "
         "with its own recipe set (Rilveg 2936 map).",
    "5": "Scenario case above the design case: simulate, scope the incremental capital against Case 1, "
         "and SteerCo decides whether Case 5 replaces Case 1 as the high design case.",
}
NOTE_1 = ("Note 1: Case 5 is simulated and a feasibility-level scope developed to compare the incremental "
          "capital needed for Case 5 with the capital for Case 1. The increment goes to SteerCo early in "
          "Concept Design to decide whether Case 5 is kept as the design basis or dropped.")
SETTINGS = {
    "titre_basis": "productivity",   # "productivity" (g/L reactor/day) or "titre" (mean permeate g/L)
    "shape": "auto",                 # "auto", a TITRE_PROFILES key, "biology" or "custom"
    "lock_titre": True,
    "lock_yield": True,
    "lock_load": True,
    "turnaround_days": 3.0,          # clean, SIP and set-up between batches (assumption)
    "operating_days": 330.0,         # available days per year (assumption)
    "envelope_with_case5": False,    # add Case 5 to the design envelope
}
BASIS_UNIT = {"productivity": "g/L/day", "titre": "g/L"}
YMIN, YMAX = 0.01, 0.999
KPI_KEYS = ("drug_substance_kg", "product_made_kg", "overall_yield_pct", "final_titre_gL", "first_cycle_n_day",
            "capture_cycles", "peak_cycles_day", "capture_wait_h", "buffer_m3", "buffer_L_per_g_ds",
            "buffer_daily_peak_L", "wfi_m3", "wastewater_m3", "wastewater_m3_per_kg_ds", "ww_daily_peak_m3",
            "depth_filters", "vf_changes", "tff_swaps", "lots", "containers", "in_process_end_g", "n1_days")
ANNUAL_ROWS = (("drug_substance_kg", "Drug substance", "kg"),
               ("product_made_kg", "Product made in the N reactor", "kg"),
               ("media_used_m3", "Media (N-1 + N)", "m³"),
               ("buffer_m3", "Solutions incl. WFI and Brij-9", "m³"),
               ("wfi_m3", "WFI used directly", "m³"),
               ("detergent_L", "Brij-9 concentrate", "L"),
               ("wastewater_m3", "Wastewater, all routes", "m³"),
               ("capture_cycles", "Protein A cycles", "cycles"),
               ("depth_filters", "Depth filters", "filters"),
               ("vf_changes", "Virus filter changes", "changes"),
               ("prefilters", "0.2 um prefilters before the virus filter", "filters"),
               ("tff_swaps", "TFF set swaps", "swaps"),
               ("lots", "Lots", "lots"),
               ("containers", "Cryovault containers", "containers"))
SKIP_UNITS = {"cm", "cm/h", "cycles", "filters", "containers"}   # not capacity items in the increment


# ── helpers ──────────────────────────────────────────────────────────
def _f(v):
    if v is None or v == "":
        return None
    try:
        x = float(v)
    except (TypeError, ValueError):
        return None
    return x if math.isfinite(x) else None


def _num(v):
    return float(v) if isinstance(v, (int, float)) and not isinstance(v, bool) and math.isfinite(v) else 0.0


def _r(x, n=4):
    return round(x, n) if isinstance(x, float) and math.isfinite(x) else x


def _rr(x):
    if not isinstance(x, (int, float)) or isinstance(x, bool) or not math.isfinite(x):
        return None
    a = abs(x)
    return round(x, 4 if a < 1 else 3 if a < 10 else 2 if a < 1000 else 1)


def settings(d=None):
    s = dict(SETTINGS)
    for k, v in (d or {}).items():
        if k not in s or v is None:
            continue
        cur = s[k]
        try:
            if isinstance(cur, bool):
                s[k] = v if isinstance(v, bool) else str(v).lower() in ("1", "true", "yes", "on")
            elif isinstance(cur, float):
                s[k] = float(v)
            else:
                s[k] = str(v)
        except (TypeError, ValueError):
            pass
    s["operating_days"] = min(max(s["operating_days"], 1.0), 366.0)
    s["turnaround_days"] = max(s["turnaround_days"], 0.0)
    if s["titre_basis"] not in BASIS_UNIT:
        s["titre_basis"] = "productivity"
    return s


def case_spec(cid, values=None):
    """Case values after edits and inheritance, with where each value came from."""
    values = values or {}

    def raw(c):
        d = {k: CASES[c][k] for k in VALUE_KEYS}
        edited = set()
        for k, v in (values.get(c) or {}).items():
            if k not in VALUE_KEYS:
                continue
            x = _f(v)
            if k == "yield_est" and x is not None and x > 1.0:
                x /= 100.0                                  # given in %
            same = (x is None and d[k] is None) or (x is not None and d[k] is not None and abs(x - d[k]) < 1e-9)
            if not same:
                edited.add(k)
            d[k] = x
        return d, edited

    d, edited = raw(cid)
    src = {}
    for k in VALUE_KEYS:
        if d[k] is not None:
            src[k] = "edited" if k in edited else "case table"
            continue
        parent = INHERIT.get(cid)
        if parent and k in ("yield_est", "capture_load_gL"):
            pd, _ = raw(parent)
            if pd[k] is not None:
                d[k] = pd[k]
                src[k] = f"inherited from Case {parent}"
                continue
        src[k] = "not in the table: Inputs value kept"
    spec = dict(CASES[cid])
    spec.update(d)
    spec["id"], spec["source"] = cid, src
    return spec


def _basis_value(p, basis):
    """Noise-free run average of the titre curve at titre_scale = 1, on the chosen basis."""
    q = copy.deepcopy(p)
    q.upstream.titre_scale = 1.0
    pts = profile_points(q)
    n = max(1, int(round(q.run_days * 24)))
    c = sum(titre_at(q, (k + 0.5) / 24.0, pts) for k in range(n)) / n
    return c * q.upstream.VVD if basis == "productivity" else c


def _pick_shape(p, target, basis):
    best = None
    for name in TITRE_PROFILES:
        q = copy.deepcopy(p)
        q.upstream.molecule = name
        v = _basis_value(q, basis)
        if v > 0:
            score = abs(math.log(target / v))
            if best is None or score < best[0]:
                best = (score, name)
    return best[1] if best else "purple"


def _yield_keys(p):
    keys = [("pcc", "capture_yield") if p.capture.mode == "pcc" else ("capture", "capture_yield"),
            ("depth", "depth_yield"), ("polish", "aex_yield"), ("polish", "cex_yield"),
            ("polish", "vf_yield"), ("sptff", "sptff_yield")]
    if p.vi.mode == "low_ph":
        keys.append(("vi", "vi_yield"))
    return keys


def _calibrate_yield(p, target):
    keys = _yield_keys(p)
    vals = {k: float(getattr(getattr(p, k[0]), k[1])) for k in keys}
    for _ in range(40):
        cur = math.prod(vals.values())
        if cur <= 0 or abs(cur / target - 1.0) < 1e-7:
            break
        up = target > cur
        free = [k for k in keys if not (up and vals[k] >= YMAX - 1e-12)]
        if not free:
            break
        f = (target / cur) ** (1.0 / len(free))
        for k in free:
            vals[k] = min(max(vals[k] * f, YMIN), YMAX)
    return vals


def _load_now(p):
    cyc = capture_cycle(p)
    return cyc["target_g"] / max(cyc["CV_L"], 1e-9), cyc["CV_L"]


# ── apply one case to a parameter set ────────────────────────────────
def apply_case(p, case):
    """Change p in place for the case; return what was done (None when no case)."""
    if not case:
        return None
    cid = str(case.get("id") or "")
    if cid not in CASES:
        return None
    st = settings(case.get("settings"))
    spec = case_spec(cid, case.get("values"))
    changes = []

    def setp(sec, key, val, why):
        obj = getattr(p, sec)
        old = getattr(obj, key)
        if isinstance(old, (int, float)) and not isinstance(old, bool) and isinstance(val, (int, float)):
            if abs(old - val) <= 1e-9 * max(1.0, abs(old)):
                return
        elif old == val:
            return
        setattr(obj, key, val)
        changes.append({"param": f"{sec}.{key}", "inputs": _r(old), "case": _r(val), "why": why})

    recipe = RECIPES.get(cid)
    if recipe == "rilveg":
        from icb_rilveg import apply_rilveg          # here to avoid a circular import
        apply_rilveg(p, verbose=False)
        changes.append({"param": "recipe", "inputs": "platform", "case": "rilveg",
                        "why": "Case 4 runs the Rilveg 2936 recipe set (icb_rilveg.py)"})

    if p.upstream.mode != "profile":
        setp("upstream", "mode", "profile", "the cases are IDB titre-curve runs")
    basis, target = st["titre_basis"], spec["avg_titre_gLd"]
    shape = st["shape"]
    if shape == "custom" and not p.upstream.custom_points:
        shape = "auto"
    if shape == "auto":
        shape = _pick_shape(p, target, basis) if target else p.upstream.molecule
    if shape not in ("custom", "biology") and shape not in TITRE_PROFILES:
        shape = "purple"
    setp("upstream", "molecule", shape, "titre-curve shape for this case")
    curve = _basis_value(p, basis)
    if st["lock_titre"] and target and curve > 0:
        setp("upstream", "titre_scale", target / curve,
             f"scales the {profile_label(p).lower()} ({curve:.3g} {BASIS_UNIT[basis]} as drawn) "
             f"to a run average of {target:g} {BASIS_UNIT[basis]}")

    y0 = train_yield(p)
    ty = spec["yield_est"]
    if st["lock_yield"] and ty:
        vals = _calibrate_yield(p, min(max(ty, YMIN), YMAX))
        for (sec, key), v in vals.items():
            setp(sec, key, v, f"step yields scaled by one factor so Protein A → fill is {ty:.0%}")

    load = spec["capture_load_gL"]
    if st["lock_load"] and load:
        if p.capture.mode == "pcc":
            setp("pcc", "DBC_gL", load / max(p.pcc.load_factor, 1e-9),
                 f"PCC load {load:g} g/L resin = DBC × load_factor")
        else:
            cv = p.capture.column_capacity_g / max(p.capture.DBC_gL, 1e-9)
            dbc = load / max(p.capture.load_util, 1e-9)
            setp("capture", "DBC_gL", dbc, f"load {load:g} g/L resin = DBC × load_util ({p.capture.load_util:g})")
            setp("capture", "column_capacity_g", dbc * cv, f"column volume kept at {cv:.2f} L per train")
    tt = titrant_plan(p)
    setp("aex", "base_frac", tt["base_frac"],
         f"2 M Tris base to take the eluate from pH {p.capture.elute_pH:g} to 7.4 (buffer only; confirm by titration)")
    setp("cex", "acid_frac", tt["acid_frac"],
         f"1 M acetic acid to take the AEX pool from pH 7.4 to {p.cex.load_pH:g} (buffer only; confirm by titration)")
    lnow, cv_now = _load_now(p)
    return {"id": cid, "name": spec["name"], "short": spec["short"], "kind": spec["kind"],
            "kind_label": KIND_LABEL[spec["kind"]], "purpose": PURPOSE[cid],
            "spec": {k: spec[k] for k in VALUE_KEYS}, "source": spec["source"], "settings": st,
            "basis": basis, "basis_unit": BASIS_UNIT[basis],
            "shape": p.upstream.molecule, "shape_label": profile_label(p),
            "curve_value": _rr(curve), "titre_scale": _r(float(p.upstream.titre_scale)),
            "plan_value": _rr(curve * p.upstream.titre_scale),
            "yield_inputs": _rr(y0), "yield_plan": _rr(train_yield(p)),
            "capture_load_gL": _rr(lnow), "capture_CV_L": _rr(cv_now), "changes": changes}


def rebase_suggestions(sug, base, q):
    """Sizing suggests column_capacity_g at the case DBC; restate it at the Inputs DBC so that
    'apply and re-run' keeps the sized column volume when the case sets DBC again."""
    if not sug or "capture" not in sug or q.capture.mode == "pcc":
        return sug
    cap = sug["capture"].get("column_capacity_g")
    if cap is None:
        return sug
    cv = cap / max(q.capture.DBC_gL, 1e-9)
    sug["capture"]["column_capacity_g"] = round(cv * base.capture.DBC_gL, 3)
    return sug


# ── per-year results, occupancy and checks ───────────────────────────
def occupancy(p, n1_days, st, batches):
    turn, op = st["turnaround_days"], st["operating_days"]
    tail = max(p.dsp.tail_h, 0.0) / 24.0
    occ = {"n1": (n1_days + turn) if p.n1.enabled else 0.0,
           "n": p.run_days + turn,
           "dsp": p.run_days + tail + turn}
    out = {"campaign_days": _rr((n1_days if p.n1.enabled else 0.0) + p.run_days + tail),
           "turnaround_days": turn, "operating_days": op, "units": [],
           "max_batches_one_line": int(op // occ["n"]) if occ["n"] > 0 else None}
    for key, label in (("n1", "N-1 seed reactor"), ("n", "N production reactor"), ("dsp", "Downstream train")):
        d = occ[key]
        if d <= 0:
            continue
        need = batches * d if batches else None
        lines = int(math.ceil(need / op - 1e-9)) if need else None
        out["units"].append({"key": key, "label": label, "days_per_batch": _rr(d),
                             "days_per_year": _rr(need) if need else None, "lines": lines,
                             "utilisation_pct": round(100 * need / (lines * op), 1) if lines else None})
    return out


def _checks(p, res, info):
    k = res.get("kpi") or {}
    loss = k.get("losses_g") or {}
    made = (k.get("product_made_kg") or 0.0) * 1000.0
    ds = (k.get("drug_substance_kg") or 0.0) * 1000.0
    to_cap = made - sum(loss.get(x, 0.0) for x in ("tff_retentate", "cell_bleed", "ramp_harvest"))
    rows = [r for r in (res.get("series") or []) if r.get("usp_on")]
    tit = [r["titre_harvest"] for r in rows if isinstance(r.get("titre_harvest"), (int, float))]
    prod = made / max(p.bioreactor.V_L, 1e-9) / max(p.run_days, 1e-9)
    mean_t = sum(tit) / len(tit) if tit else None
    spec = (info or {}).get("spec") or {}
    basis = (info or {}).get("basis", "productivity")
    tgt, ty, load = spec.get("avg_titre_gLd"), spec.get("yield_est"), spec.get("capture_load_gL")
    out = []

    def add(key, label, unit, target, value, tol=None, note=""):
        d = {"key": key, "label": label, "unit": unit, "target": _rr(target), "value": _rr(value), "note": note}
        if target and value is not None and tol:
            d["delta_pct"] = _rr(100 * (value / target - 1))
            d["ok"] = abs(value / target - 1) <= tol
        out.append(d)

    add("productivity", "Average productivity (made per L of N volume per day)", "g/L/day",
        tgt if basis == "productivity" else None, prod, 0.05)
    add("titre", "Average permeate titre", "g/L", tgt if basis == "titre" else None, mean_t, 0.05)
    add("yield_plan", "Train yield, plan (product of the step yields)", "%", 100 * ty if ty else None,
        100 * train_yield(p), 0.01)
    add("yield_sim", "Train yield, simulated (DS ÷ permeate sent to capture)", "%", 100 * ty if ty else None,
        100 * ds / to_cap if to_cap > 0 else None, 0.05,
        "lower than plan when product is still in process at the end or columns are overloaded")
    add("yield_overall", "Overall yield, simulated (DS ÷ product made)", "%", None, k.get("overall_yield_pct"),
        note="includes the TFF cell concentrate and any permeate diverted before capture")
    lnow, _ = _load_now(p)
    add("load", "Protein A load", "g/L resin", load, lnow, 0.01)
    return out


def annual(p, res, info):
    st = info["settings"] if info else settings(None)
    b = (info or {}).get("spec", {}).get("batches_per_year")
    k = res.get("kpi") or {}
    rows = []
    for key, label, unit in ANNUAL_ROWS:
        v = k.get(key)
        if v is None:
            continue
        rows.append({"label": label, "unit": unit, "per_batch": v, "per_year": _rr(v * b) if b else None})
    for route, m3 in ((res.get("wastewater") or {}).get("by_route_m3") or {}).items():
        rows.append({"label": f"Wastewater to {route}", "unit": "m³", "per_batch": m3,
                     "per_year": _rr(m3 * b) if b else None})
    cyc = k.get("capture_cycles") or 0
    if cyc:
        pcc = p.capture.mode == "pcc"
        n = max(1, int(p.pcc.n_columns if pcc else p.capture.trains))
        rmax = max(1, int(p.pcc.resin_cycles_max if pcc else p.capture.resin_cycles_max))
        rows.append({"label": "Protein A resin lifetimes used per column / train", "unit": "lifetimes",
                     "per_batch": _rr(cyc / n / rmax), "per_year": _rr(cyc * b / n / rmax) if b else None})
    n1_days = ((res.get("design") or {}).get("n1") or {}).get("days") or 0.0
    return {"batches_per_year": b, "rows": rows, "occupancy": occupancy(p, n1_days, st, b),
            "checks": _checks(p, res, info)}


# ── all cases ────────────────────────────────────────────────────────
def case_api_info():
    return {"cases": {cid: dict(c, id=cid, purpose=PURPOSE[cid], inherit=INHERIT.get(cid),
                                recipe=RECIPES.get(cid))
                      for cid, c in CASES.items()},
            "order": list(CASES), "design_cases": list(DESIGN_CASES), "settings": SETTINGS,
            "shapes": dict({k: v["label"] for k, v in TITRE_PROFILES.items()}, **EXTRA_PROFILES), "note1": NOTE_1,
            "kind_labels": KIND_LABEL, "basis_units": BASIS_UNIT}


def preview_cases(p_base, case=None):
    """Noise-free plan for every case (fast: no simulation)."""
    st = settings((case or {}).get("settings"))
    values = (case or {}).get("values")
    out = []
    for cid in CASES:
        q = copy.deepcopy(p_base)
        info = apply_case(q, {"id": cid, "settings": st, "values": values})
        d = process_design(q)
        n1d = (d.get("n1") or {}).get("days") or 0.0
        out.append({"id": cid, "info": info,
                    "design": {k: d.get(k) for k in ("peak_g_day", "peak_day", "peak_cycles_day",
                                                     "max_cycles_day_at_peak", "first_cycle_n_day",
                                                     "first_cycle_sump_L", "peak_buffer_L_day",
                                                     "total_g_harvested", "ds_kg_est")},
                    "notes": (d.get("notes") or [])[1:],
                    "occupancy": occupancy(q, n1d, st, info["spec"].get("batches_per_year"))})
    return {"cases": out, "settings": st}


def merge_sizing(by_case):
    """Worst case per sizing row across cases ('set_by' names the case)."""
    best, order, gnames, sug, notes = {}, [], [], {}, []
    for cid, sz in by_case.items():
        if not sz or sz.get("error"):
            notes.append(f"Case {cid}: {(sz or {}).get('error', 'no sizing')}")
            continue
        for g in sz.get("groups", []):
            if g["name"] not in gnames:
                gnames.append(g["name"])
            for r in g["rows"]:
                key = (g["name"], r["tag"], r["item"])
                cur = best.get(key)
                if cur is None:
                    order.append(key)
                if cur is None or _num(r.get("required")) > _num(cur.get("required")):
                    best[key] = dict(r, set_by="Case " + cid)
        for sec, kv in (sz.get("suggested_params") or {}).items():
            for k, v in kv.items():
                if _num(v) > sug.setdefault(sec, {}).get(k, -math.inf):
                    sug[sec][k] = v
    groups = [{"name": g, "rows": [best[k] for k in order if k[0] == g]} for g in gnames]
    notes.insert(0, "Design envelope: each item is sized for the case that needs the most ('Set by'). "
                    "Apply the sized values and re-run the cases: bigger columns or filters move the vessel peaks.")
    return {"groups": groups, "suggested_params": sug, "notes": notes}


def case_increment(sz1, sz5, a1, a5):
    """Case 5 against Case 1: capacity ratio and six-tenths scaling factor per item."""
    if not sz1 or not sz5 or sz1.get("error") or sz5.get("error"):
        return None
    idx = {(g["name"], r["tag"], r["item"]): r for g in sz1.get("groups", []) for r in g["rows"]}
    rows = []
    for g in sz5.get("groups", []):
        for r in g["rows"]:
            r1 = idx.get((g["name"], r["tag"], r["item"]))
            if r1 is None or r.get("unit") in SKIP_UNITS:
                continue
            q1, q5 = _num(r1.get("required")), _num(r.get("required"))
            if q1 <= 0 or q5 <= 0:
                continue
            ratio = q5 / q1
            rows.append({"group": g["name"], "tag": r["tag"], "item": r["item"], "unit": r["unit"],
                         "case1": r1.get("required"), "case5": r.get("required"),
                         "rec1": r1.get("recommended"), "rec5": r.get("recommended"),
                         "ratio": round(ratio, 3), "factor": round(ratio ** 0.6, 3)})
    u1 = {u["key"]: u for u in (a1 or {}).get("occupancy", {}).get("units", [])}
    for u in (a5 or {}).get("occupancy", {}).get("units", []):
        o = u1.get(u["key"])
        if o and o.get("lines") and u.get("lines"):
            ratio = u["lines"] / o["lines"]
            rows.append({"group": "Production lines", "tag": "", "item": u["label"] + " lines", "unit": "lines",
                         "case1": o["lines"], "case5": u["lines"], "rec1": o["lines"], "rec5": u["lines"],
                         "ratio": round(ratio, 3), "factor": round(ratio, 3)})
    rows.sort(key=lambda x: -x["ratio"])
    grow = sum(1 for r in rows if r["ratio"] > 1.05)
    same = sum(1 for r in rows if 0.95 <= r["ratio"] <= 1.05)
    return {"rows": rows, "grow": grow, "same": same, "shrink": len(rows) - grow - same,
            "notes": [f"{grow} items need more capacity for Case 5 than for Case 1 (ratio > 1.05), "
                      f"{same} are about the same.",
                      "Ratio = required capacity Case 5 ÷ Case 1 from the same Inputs equipment. "
                      "Factor = ratio^0.6 (six-tenths rule) for equipment, whole-line ratio for production "
                      "lines: a scaling index for the estimate, not a cost.",
                      "Consumables (resin cycles, filters, containers) are in the per-year table, not here."]}


def compare_cases(p_base, case=None, seed=0, ids=None):
    from icb_dynamic import ICBSimulator           # here to avoid a circular import
    st = settings((case or {}).get("settings"))
    values = (case or {}).get("values")
    ids = [c for c in (ids or list(CASES)) if c in CASES]
    out, sizing, ann = [], {}, {}
    for cid in ids:
        q = copy.deepcopy(p_base)
        info = apply_case(q, {"id": cid, "settings": st, "values": values})
        res = ICBSimulator(q, forecaster=None, seed=seed).run()
        a = annual(q, res, info)
        sz = res.get("sizing") or {}
        rebase_suggestions(sz.get("suggested_params"), p_base, q)
        sizing[cid], ann[cid] = sz, a
        k = res.get("kpi") or {}
        out.append({"id": cid, "info": info, "kpi": {key: k.get(key) for key in KPI_KEYS},
                    "alarms": k.get("alarms") or {}, "annual": a,
                    "design_notes": (res.get("design") or {}).get("notes", [])[1:]})
    env_ids = [c for c in DESIGN_CASES if c in sizing] + (["5"] if st["envelope_with_case5"] and "5" in sizing else [])
    return {"ok": True, "cases": out, "settings": st, "envelope_cases": env_ids,
            "envelope": merge_sizing({c: sizing[c] for c in env_ids}),
            "increment": case_increment(sizing.get("1"), sizing.get("5"), ann.get("1"), ann.get("5")),
            "note1": NOTE_1}
