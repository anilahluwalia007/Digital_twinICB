"""
icb_process.py — Integrated Dynamic BioProcessing (IDB) for mAbs:
flowsheet, parameters, titre curves and design checks.
NavNevAI_mAbWasteWater · CarbonTesseract 2026

IDB (site definition)
  N-1   perfusion seed: media in = cell-free permeate out to waste; cells stay,
        so the seed reaches a high density.
  N     continuous media feed; outflow = cells + broth + product (in = out).
  TFF   removes all cells from the N outflow (cell concentrate -> bio-kill);
        permeate (broth + product) -> harvest sump.
  DSP   standard Protein A from the sump, one cycle each time the sump holds
        one column load -> low-pH VI -> AEX + CEX (polish.order) -> virus filter -> SPTFF.
        3-column PCC is kept as an option.

The permeate titre curve of each molecule (Figure 2) governs the process:
titre changes every hour, so buffer and wastewater are planned in litres per
gram recovered, not per day.

Upstream modes
  profile      titre(t) from the molecule's permeate curve (IDB, default)
  mechanistic  ATF perfusion with cell return and bleed control (legacy,
               fully continuous ICB, not IDB; kept for comparison)

Units: time h, volume L, mass g, cell density 1e6 cells/mL (= 1e9 cells/L),
specific rates pg/cell/h. Defaults are illustrative, not site data.
"""
from __future__ import annotations
import bisect
import copy
import math
from dataclasses import dataclass, field, asdict, is_dataclass

# Permeate titre (g/L) against N day, digitised from Figure 2. Replace with the
# raw data. The Excel polynomial fits are not used: they only hold to day 14
# (the purple cubic peaks at day 13.9 and reaches zero near day 21).
# After the last point the titre is held at the last value.
TITRE_PROFILES = {
    "purple": {"label": "Purple molecule", "points": [
        [0.0, 0.40], [0.25, 0.39], [0.75, 0.39], [1.25, 0.43], [1.75, 0.51], [2.25, 0.63],
        [2.75, 0.77], [3.25, 0.96], [3.75, 1.15], [4.25, 1.36], [4.75, 1.59], [5.25, 1.83],
        [5.75, 2.07], [6.25, 2.32], [6.75, 2.57], [7.25, 2.82], [7.75, 3.07], [8.25, 3.31],
        [8.75, 3.53], [9.25, 3.75], [9.75, 3.96], [10.25, 4.15], [10.75, 4.31], [11.25, 4.48],
        [11.75, 4.61], [12.25, 4.72], [12.75, 4.82], [13.25, 4.90], [13.75, 4.97]]},
    "green": {"label": "Green molecule", "points": [
        [0.0, 0.23], [0.25, 0.23], [1.25, 0.26], [1.75, 0.30], [2.25, 0.37], [2.75, 0.46],
        [3.25, 0.57], [3.75, 0.68], [4.25, 0.82], [4.75, 0.96], [5.25, 1.10], [5.75, 1.25],
        [6.25, 1.40], [6.75, 1.55], [7.25, 1.69], [7.75, 1.84], [8.25, 1.98], [8.75, 2.12],
        [9.25, 2.25], [9.75, 2.37], [10.25, 2.48], [10.75, 2.59], [11.25, 2.69], [11.75, 2.77],
        [12.25, 2.84], [12.75, 2.89], [13.25, 2.95], [13.75, 2.98]]},
    "blue": {"label": "Blue molecule", "points": [
        [0.0, 0.01], [0.25, 0.02], [0.75, 0.05], [1.25, 0.09], [1.75, 0.12], [2.25, 0.16],
        [2.75, 0.20], [3.25, 0.23], [3.75, 0.27], [4.25, 0.32], [4.75, 0.39], [5.25, 0.47],
        [5.75, 0.56], [6.25, 0.67], [6.75, 0.75], [7.25, 0.93], [7.75, 0.96], [8.25, 1.02],
        [8.75, 1.08], [9.25, 1.19], [9.75, 1.27], [10.25, 1.39], [10.75, 1.48], [11.25, 1.51],
        [11.75, 1.46], [12.25, 1.44], [12.75, 1.42], [13.25, 1.42], [13.75, 1.41]]},
    "red": {"label": "Red molecule", "points": [
        [0.0, 0.04], [1.0, 0.10], [2.0, 0.20], [3.0, 0.41], [4.0, 0.67], [5.0, 0.64],
        [6.0, 0.68], [7.0, 0.67], [8.0, 0.67], [9.0, 0.76], [10.0, 0.77], [11.0, 0.79],
        [12.0, 0.76], [13.0, 0.75], [14.0, 0.74]]},
}


@dataclass
class N1Params:
    enabled: bool = True
    V_L: float = 15.0                # N-1 working volume
    VVD: float = 2.0                 # media in = cell-free permeate out to waste
    days: float = 5.0
    X0: float = 0.5                  # N-1 inoculation VCD
    mu_h: float = 0.025              # net growth with cells retained
    X_max: float = 40.0              # density cap
    transfer_L: float = 15.0         # volume moved into N at inoculation


@dataclass
class UpstreamParams:
    mode: str = "profile"            # "profile" (IDB) or "mechanistic" (legacy ATF model)
    molecule: str = "purple"         # key in TITRE_PROFILES, or "custom"
    custom_points: list = field(default_factory=list)   # [[N day, g/L], ...]
    VVD: float = 1.2                 # N feed = N outflow (cells + broth + product)
    tff_retentate_frac: float = 0.02 # share of N outflow leaving the TFF as cell concentrate
    titre_scale: float = 1.0         # what-if multiplier on the curve
    time_stretch: float = 1.0        # >1 slows the curve down
    variation_sigma: float = 0.04    # slow random titre drift (log scale)
    variation_tau_h: float = 48.0


@dataclass
class BioreactorParams:
    V_L: float = 83.0                # N working volume (both modes)
    # --- legacy ATF perfusion model only ---
    X0: float = 0.5                  # used when N-1 is disabled
    X_setpoint: float = 50.0
    X_max: float = 90.0
    mu_max_h: float = 0.028
    Ks_gL: float = 0.5
    kd_h: float = 0.0015
    k_lysis_h: float = 0.005
    qP_pg_cell_h: float = 1.2
    qp_ou_sigma: float = 0.05
    qp_ou_tau_h: float = 72.0
    qS_pg_cell_h: float = 1.5
    S_feed_gL: float = 8.0
    S0_gL: float = 4.0
    VVD_min: float = 0.3
    VVD_max: float = 2.5
    bleed_max_VVD: float = 0.5
    sieving0: float = 0.98
    sieve_loss_per_day: float = 0.006
    sieve_replace_below: float = 0.85
    atf_swap_h: int = 2
    atf_swap_flush_L: float = 20.0


@dataclass
class CaptureParams:
    mode: str = "batch"              # "batch" (standard Protein A from the sump) or "pcc"
    column_capacity_g: float = 100.0 # mass the column binds (DBC x CV)
    DBC_gL: float = 40.0             # -> CV = 2.5 L
    load_util: float = 0.90          # load target as a share of capacity (margin >= 2x titre noise)
    load_flow_CVh: float = 15.0      # 4 min residence
    equil_CV: float = 3.0
    wash_CV: float = 6.0
    elute_CV: float = 3.0
    eluate_CV: float = 2.0           # collected product peak
    strip_CV: float = 2.0
    cip_CV: float = 3.0
    reequil_CV: float = 5.0
    recovery_flow_CVh: float = 10.0
    cip_contact_h: float = 0.25
    capture_yield: float = 0.97
    sump_max_L: float = 600.0
    max_sump_hold_h: float = 96.0    # product stability limit in the sump
    capture_start_day: float = 0.0   # N day when permeate first goes to the sump
    resin_cycles_max: int = 200


@dataclass
class SurgeParams:                   # PCC mode: level-controlled harvest surge
    V_max_L: float = 40.0
    level_sp: float = 0.4


@dataclass
class PCCParams:
    n_columns: int = 3
    CV_L: float = 0.5
    DBC_gL: float = 45.0
    load_factor: float = 1.1
    sat_factor: float = 1.3
    bt_frac: float = 0.5
    capture_yield: float = 0.97
    wash_CV: float = 6.0
    elute_CV: float = 3.0
    eluate_CV: float = 2.0
    strip_CV: float = 2.0
    cip_CV: float = 3.0
    reequil_CV: float = 5.0
    recovery_flow_CVh: float = 20.0
    cip_contact_h: float = 0.25
    load_flow_max_Lh: float = 12.0
    min_residence_min: float = 2.0
    resin_cycles_max: int = 200

    def recovery_time_h(self) -> float:
        cvs = self.wash_CV + self.elute_CV + self.strip_CV + self.cip_CV + self.reequil_CV
        return cvs / self.recovery_flow_CVh + self.cip_contact_h


@dataclass
class VIParams:
    eluate_surge_max_L: float = 20.0
    heel_L: float = 0.5
    tube_V_L: float = 2.0
    t_hold_min_h: float = 1.0
    eluate_pH_mean: float = 3.9
    eluate_pH_sd: float = 0.05
    pH_target: float = 3.6
    pH_max: float = 3.8
    ph_drop_per_acid_frac: float = 10.0
    acid_bias: float = 0.03
    acid_max: float = 0.15
    base_frac: float = 0.06
    vi_yield: float = 0.99


@dataclass
class PolishParams:
    cex_yield: float = 0.92
    cex_buffer_L_per_g: float = 0.25
    cex_pool_gL: float = 10.0
    aex_yield: float = 0.98
    aex_buffer_L_per_g: float = 0.05
    aex_flush_frac: float = 0.05
    vf_yield: float = 0.99
    vf_capacity_L: float = 2000.0
    vf_flush_L: float = 2.0
    order: str = "aex_cex"           # "aex_cex" (AEX then CEX, as on the site PFD) or "cex_aex" (legacy)


@dataclass
class SPTFFParams:
    target_gL: float = 50.0
    DV: float = 7.0
    sptff_yield: float = 0.98
    excipient_frac: float = 0.05


@dataclass
class PATParams:
    vcd_noise: float = 0.03
    titre_noise: float = 0.05
    ph_noise: float = 0.02


@dataclass
class ControlParams:
    CSPR_pL_cell_day: float = 30.0   # legacy model
    bleed_kc: float = 0.07           # legacy model, L/h per 1e6 cells/mL
    bleed_ti_h: float = 24.0
    bleed_enable_frac: float = 0.95
    surge_tau_h: float = 4.0         # PCC mode
    eluate_tau_h: float = 3.0
    acid_kc: float = 0.05
    acid_ti_h: float = 6.0
    viability_alarm: float = 0.80
    glucose_alarm_gL: float = 0.3
    titre_drop_alarm: float = 0.15   # NevAI early-warning threshold


def _cast(old, new):
    if isinstance(old, bool):
        return new if isinstance(new, bool) else str(new).lower() in ("1", "true", "yes", "on")
    if isinstance(old, int):
        return int(float(new))
    if isinstance(old, float):
        return float(new)
    if isinstance(old, str):
        return str(new)
    if isinstance(old, list):
        return list(new) if isinstance(new, (list, tuple)) else old
    return new


@dataclass
class ICBParams:
    n1: N1Params = field(default_factory=N1Params)
    upstream: UpstreamParams = field(default_factory=UpstreamParams)
    bioreactor: BioreactorParams = field(default_factory=BioreactorParams)
    capture: CaptureParams = field(default_factory=CaptureParams)
    surge: SurgeParams = field(default_factory=SurgeParams)
    pcc: PCCParams = field(default_factory=PCCParams)
    vi: VIParams = field(default_factory=VIParams)
    polish: PolishParams = field(default_factory=PolishParams)
    sptff: SPTFFParams = field(default_factory=SPTFFParams)
    pat: PATParams = field(default_factory=PATParams)
    control: ControlParams = field(default_factory=ControlParams)
    run_days: float = 14.0           # N days (N-1 days come before)
    lot_window_h: float = 24.0
    disturbances: list = field(default_factory=list)

    @classmethod
    def from_dict(cls, d=None):
        p = cls()
        for key, val in (d or {}).items():
            if not hasattr(p, key) or val is None:
                continue
            cur = getattr(p, key)
            if is_dataclass(cur) and isinstance(val, dict):
                for k2, v2 in val.items():
                    if hasattr(cur, k2) and v2 is not None and not callable(getattr(cur, k2)):
                        try:
                            setattr(cur, k2, _cast(getattr(cur, k2), v2))
                        except (TypeError, ValueError):
                            pass
            elif key == "disturbances" and isinstance(val, list):
                p.disturbances = [dict(x) for x in val if isinstance(x, dict)]
            else:
                try:
                    setattr(p, key, _cast(cur, val))
                except (TypeError, ValueError):
                    pass
        return p

    def to_dict(self):
        return asdict(self)


# name: (area, route, driver, description)
# driver: volume = follows the feed rate; mass = follows grams recovered;
#         event = start-up, filter changes and upsets
WASTE_STREAMS = {
    "n1_permeate":         ("USP", "direct",         "volume", "N-1 permeate: spent media removed to intensify the seed"),
    "tff_retentate":       ("USP", "bio-kill",       "volume", "TFF cell concentrate from the N outflow"),
    "cell_bleed":          ("USP", "bio-kill",       "volume", "Cell bleed (legacy ATF model)"),
    "atf_flush":           ("USP", "bio-kill",       "event",  "ATF filter swap flush (legacy ATF model)"),
    "ramp_harvest":        ("USP", "direct",         "volume", "Permeate sent to drain before capture starts"),
    "sump_overflow":       ("USP", "direct",         "event",  "Harvest sump overflow (cell-free)"),
    "capture_flowthrough": ("DSP", "direct",         "volume", "Protein A load flow-through: spent permeate, high organic load"),
    "capture_buffers":     ("DSP", "direct",         "mass",   "Protein A equilibration, wash, re-equilibration, elution pre/post-peak"),
    "capture_strip_cip":   ("DSP", "neutralisation", "mass",   "Protein A acid strip and caustic CIP"),
    "vi_flush":            ("DSP", "neutralisation", "event",  "VI tube start-up buffer"),
    "vi_divert":           ("DSP", "neutralisation", "event",  "Out-of-spec VI material diverted"),
    "cex_flowthrough":     ("DSP", "direct",         "mass",   "CEX load flow-through (AEX pool; VI pool in legacy order)"),
    "cex_buffers":         ("DSP", "neutralisation", "mass",   "CEX buffers incl. strip and CIP"),
    "aex_buffers":         ("DSP", "direct",         "mass",   "AEX membrane equilibration and flush"),
    "vf_flush":            ("DSP", "direct",         "mass",   "Virus filter pre-use flush"),
    "sptff_permeate":      ("DSP", "direct",         "mass",   "SPTFF permeate incl. inline DF buffer"),
}
BUFFER_OPS = ("capture", "vi", "cex", "aex", "vf", "sptff")


def _r(x, n=3):
    return round(x, n) if isinstance(x, (int, float)) and math.isfinite(x) else None


# ── titre curves ──────────────────────────────────────────────────────
def profile_points(p):
    up = p.upstream
    if up.molecule == "custom" and up.custom_points:
        try:
            pts = sorted([float(a), float(b)] for a, b in up.custom_points)
            if pts:
                return pts
        except (TypeError, ValueError):
            pass
    return TITRE_PROFILES.get(up.molecule, TITRE_PROFILES["purple"])["points"]


def profile_label(p):
    if p.upstream.molecule == "custom" and p.upstream.custom_points:
        return "Custom curve"
    return TITRE_PROFILES.get(p.upstream.molecule, TITRE_PROFILES["purple"])["label"]


def titre_at(p, n_day, pts=None):
    """Permeate titre (g/L) at N day n_day: linear interpolation, held after the last point."""
    pts = pts or profile_points(p)
    up = p.upstream
    t = n_day / max(up.time_stretch, 1e-6)
    xs = [x for x, _ in pts]
    if t <= xs[0]:
        v = pts[0][1]
    elif t >= xs[-1]:
        v = pts[-1][1]
    else:
        i = bisect.bisect_right(xs, t)
        (x0, y0), (x1, y1) = pts[i - 1], pts[i]
        v = y0 + (y1 - y0) * (t - x0) / (x1 - x0) if x1 > x0 else y1
    return max(0.0, v * up.titre_scale)


# ── unit summaries ────────────────────────────────────────────────────
def capture_cycle(p):
    """Per-cycle volumes and times for the selected capture mode."""
    if p.capture.mode == "pcc":
        pc = p.pcc
        cv = pc.CV_L
        buf = (pc.wash_CV + pc.elute_CV + pc.strip_CV + pc.cip_CV + pc.reequil_CV) * cv
        strip_cip = (pc.strip_CV + pc.cip_CV) * cv
        return {"mode": "pcc", "CV_L": cv, "target_g": pc.load_factor * pc.DBC_gL * cv,
                "capacity_g": pc.sat_factor * pc.DBC_gL * cv, "buffer_L": buf,
                "eluate_L": pc.eluate_CV * cv, "strip_cip_L": strip_cip,
                "other_buffer_waste_L": buf - pc.eluate_CV * cv - strip_cip,
                "load_flow_Lh": pc.load_flow_max_Lh, "recovery_h": pc.recovery_time_h(),
                "yield": pc.capture_yield, "resin_cycles_max": pc.resin_cycles_max}
    cp = p.capture
    cv = cp.column_capacity_g / max(cp.DBC_gL, 1e-9)
    cvs = cp.equil_CV + cp.wash_CV + cp.elute_CV + cp.strip_CV + cp.cip_CV + cp.reequil_CV
    buf = cvs * cv
    strip_cip = (cp.strip_CV + cp.cip_CV) * cv
    return {"mode": "batch", "CV_L": cv, "target_g": cp.column_capacity_g * cp.load_util,
            "capacity_g": cp.column_capacity_g, "buffer_L": buf,
            "eluate_L": cp.eluate_CV * cv, "strip_cip_L": strip_cip,
            "other_buffer_waste_L": buf - cp.eluate_CV * cv - strip_cip,
            "load_flow_Lh": max(cp.load_flow_CVh * cv, 1e-6),
            "recovery_h": cvs / max(cp.recovery_flow_CVh, 1e-9) + cp.cip_contact_h,
            "yield": cp.capture_yield, "resin_cycles_max": cp.resin_cycles_max}


def train_yield(p, cyc=None):
    cyc = cyc or capture_cycle(p)
    vi, po, sp = p.vi, p.polish, p.sptff
    return (cyc["yield"] * vi.vi_yield * po.cex_yield * po.aex_yield
            * po.vf_yield * sp.sptff_yield)


def n1_summary(p):
    n1, br = p.n1, p.bioreactor
    if not n1.enabled:
        return {"enabled": False, "n_seed_vcd": round(br.X0, 3), "waste_L": 0.0, "media_L": 0.0}
    hours = int(round(n1.days * 24))
    X = min(n1.X0 * math.exp(n1.mu_h * hours), n1.X_max)
    seed = X * min(n1.transfer_L, n1.V_L) / br.V_L
    media = n1.VVD * n1.V_L / 24.0 * hours
    return {"enabled": True, "days": n1.days, "V_L": n1.V_L, "VVD": n1.VVD,
            "final_vcd": round(X, 2), "n_seed_vcd": round(seed, 3),
            "media_L": round(media, 1), "waste_L": round(media, 1)}


def buffer_L_per_g(p):
    """Buffer litres per gram of product entering capture: the IDB planning constants."""
    vi, po, sp = p.vi, p.polish, p.sptff
    cyc = capture_cycle(p)
    tg = max(cyc["target_g"], 1e-9)
    elu = cyc["eluate_L"] / tg                       # eluate L per g loaded
    g1 = cyc["yield"] * vi.vi_yield                  # g after capture + VI
    if po.order != "cex_aex":                        # AEX flow-through, then CEX bind/elute
        aex = g1 * po.aex_buffer_L_per_g
        g2 = g1 * po.aex_yield
        cex = g2 * po.cex_buffer_L_per_g
        g3 = g2 * po.cex_yield
        v3 = g3 / po.cex_pool_gL                     # CEX pool feeds the virus filter
    else:                                            # legacy: CEX, then AEX
        cex = g1 * po.cex_buffer_L_per_g
        g2 = g1 * po.cex_yield
        aex = g2 * po.aex_buffer_L_per_g
        g3 = g2 * po.aex_yield
        v3 = g2 / po.cex_pool_gL * (1 + po.aex_flush_frac)
    g_ds = g3 * po.vf_yield * sp.sptff_yield
    v_ret = g_ds / sp.target_gL
    d = {
        "capture": cyc["buffer_L"] / tg,
        "vi": elu * vi.acid_bias + elu * (1 + vi.acid_bias) * vi.base_frac,
        "cex": cex,
        "aex": aex,
        "vf": v3 * po.vf_flush_L / max(po.vf_capacity_L, 1e-9),
        "sptff": v_ret * (sp.DV + sp.excipient_frac),
    }
    d["total"] = sum(d.values())
    d["ds_g_per_g"] = g_ds
    return d


# ── design checks ─────────────────────────────────────────────────────
def steady_state_design(p):
    """Legacy ATF perfusion model: algebraic steady state at the VCD setpoint."""
    br, vi, ct = p.bioreactor, p.vi, p.control
    cyc = capture_cycle(p)
    V, X = br.V_L, br.X_setpoint
    notes = ["Legacy ATF perfusion model (cells returned to the reactor): fully continuous ICB, not IDB."]

    Fin = ct.CSPR_pL_cell_day * X * V * 1e-3 / 24.0
    Fin = min(max(Fin, br.VVD_min * V / 24.0), br.VVD_max * V / 24.0)
    S = br.S_feed_gL - br.qS_pg_cell_h * X * 1e-3 * V / Fin
    if S <= 0:
        notes.append("Medium cannot supply enough glucose at this VCD; raise CSPR or feed glucose.")
        S = 1e-6
    mu = br.mu_max_h * S / (br.Ks_gL + S) * max(0.0, 1.0 - X / br.X_max)
    net = mu - br.kd_h
    if net <= 0:
        notes.append("Net growth is zero or negative at the setpoint; VCD cannot be held there.")
    Fb = min(max(net, 0.0) * V, Fin)
    if Fb > br.bleed_max_VVD * V / 24.0:
        notes.append("Required bleed exceeds the bleed pump limit (bleed_max_VVD).")
    Fh = Fin - Fb
    Xd = br.kd_h * X / (Fb / V + br.k_lysis_h)
    viab = X / (X + Xd)
    P = br.qP_pg_cell_h * X * 1e-3 * V / max(Fh * br.sieving0 + Fb, 1e-9)
    titre = P * br.sieving0
    harvest_g_h = Fh * titre

    target = cyc["target_g"]
    t_fill = target / harvest_g_h if harvest_g_h > 0 else math.inf
    if cyc["mode"] == "pcc":
        pc = p.pcc
        if pc.n_columns < 3:
            notes.append("PCC needs at least 3 columns for uninterrupted loading.")
        elif cyc["recovery_h"] > (pc.n_columns - 2) * t_fill:
            notes.append(f"PCC recovery ({cyc['recovery_h']:.1f} h) exceeds the time available "
                         f"({(pc.n_columns - 2) * t_fill:.1f} h); loading will stall.")
        res_min = pc.CV_L / Fh * 60.0 if Fh > 0 else math.inf
        if res_min < pc.min_residence_min:
            notes.append("Load residence time per column is below the minimum; use larger columns.")
        if pc.load_factor > pc.sat_factor:
            notes.append("load_factor above sat_factor: product will be lost in breakthrough.")
    else:
        load_h = target / titre / cyc["load_flow_Lh"] if titre > 0 else math.inf
        if load_h + cyc["recovery_h"] > t_fill:
            notes.append(f"Batch capture cycle ({load_h + cyc['recovery_h']:.1f} h) is longer than the "
                         f"time to collect one load ({t_fill:.1f} h); the sump will keep growing.")

    eluate_Lh = cyc["eluate_L"] / t_fill if math.isfinite(t_fill) and t_fill > 0 else 0.0
    vi_res = vi.tube_V_L / (eluate_Lh * (1 + vi.acid_bias)) if eluate_Lh > 0 else math.inf
    if vi_res < vi.t_hold_min_h:
        notes.append("VI tube too small for the eluate flow; hold time below minimum.")
    if len(notes) == 1:
        notes.append("Design is consistent at steady state.")
    ds_g_h = harvest_g_h * train_yield(p, cyc)
    return {
        "mode": "mechanistic",
        "perfusion_VVD": _r(Fin * 24 / V), "feed_Lh": _r(Fin, 2), "bleed_Lh": _r(Fb, 2),
        "bleed_pct_of_feed": _r(100 * Fb / Fin, 1), "harvest_Lh": _r(Fh, 2),
        "glucose_gL": _r(S, 2), "growth_rate_h": _r(mu, 4),
        "doubling_h": _r(math.log(2) / mu, 1) if mu > 0 else None,
        "viability_pct": _r(100 * viab, 1), "harvest_titre_gL": _r(titre),
        "harvest_g_day": _r(harvest_g_h * 24, 0),
        "capture_fill_h": _r(t_fill, 2), "capture_recovery_h": _r(cyc["recovery_h"], 2),
        "capture_cycles_per_day": _r(24 / t_fill, 2) if math.isfinite(t_fill) else 0,
        "eluate_conc_gL": _r(target * cyc["yield"] / cyc["eluate_L"], 1),
        "vi_residence_h": _r(vi_res, 2), "ds_kg_day": _r(ds_g_h * 24 / 1000),
        "notes": notes,
    }


def profile_design(p):
    """IDB planning view from the titre curve (noise-free): daily grams, cycles and buffer."""
    up, br, cp = p.upstream, p.bioreactor, p.capture
    cyc = capture_cycle(p)
    per_g = buffer_L_per_g(p)
    pts = profile_points(p)
    target = cyc["target_g"]
    F = up.VVD * br.V_L / 24.0
    Fr = F * up.tff_retentate_frac
    Fp = F - Fr
    notes = []

    daily, total_g = [], 0.0
    sump_g = sump_L = 0.0
    first = None
    for d in range(int(math.ceil(p.run_days))):
        g_d = 0.0
        for k in range(24):
            t = d + k / 24.0
            if t >= p.run_days:
                break
            if t < cp.capture_start_day:
                continue
            g = Fp * titre_at(p, t + 0.5 / 24.0, pts)
            g_d += g
            if first is None:
                sump_g += g
                sump_L += Fp
                if sump_g >= target:
                    first = {"n_day": round(t + 1.0 / 24.0, 2), "sump_L": round(sump_L, 1)}
        total_g += g_d
        daily.append({"n_day": d + 1, "titre_gL": round(titre_at(p, d + 0.5, pts), 3),
                      "g": round(g_d, 1), "cycles": round(g_d / target, 2),
                      "capture_buffer_L": round(g_d * per_g["capture"], 1),
                      "buffer_L": round(g_d * per_g["total"], 1)})

    pk = max(daily, key=lambda r: r["g"]) if daily else {"g": 0.0, "n_day": None}
    peak_g = pk["g"]
    C_pk = peak_g / (Fp * 24.0) if Fp > 0 else 0.0
    if cyc["mode"] == "batch":
        net_g = target
        load_h = (target / C_pk) / cyc["load_flow_Lh"] if C_pk > 0 else math.inf
        cycle_h = load_h + cyc["recovery_h"]
        if 1.0 - cp.load_util < 2.0 * p.pat.titre_noise - 1e-9:
            notes.append(f"Load margin ({100 * (1 - cp.load_util):.0f} %) is less than twice the titre "
                         f"measurement noise ({100 * p.pat.titre_noise:.0f} %); expect occasional column "
                         "overloads. Lower load_util or improve the titre measurement.")
    else:
        pc = p.pcc
        # mass carried over on column 2 is eluted in the next cycle, not this one
        net_g = target - pc.bt_frac * max(0.0, target - pc.DBC_gL * pc.CV_L)
        cycle_h = cyc["recovery_h"] / max(pc.n_columns - 2, 1)
        if pc.n_columns < 3:
            notes.append("PCC needs at least 3 columns for uninterrupted loading.")
        if pc.load_factor > pc.sat_factor:
            notes.append("load_factor above sat_factor: product will be lost in breakthrough.")
        if Fp > cyc["load_flow_Lh"]:
            notes.append("Permeate flow exceeds the PCC load pump limit.")
    peak_cycles = peak_g / max(net_g, 1e-9)
    max_cyc = 24.0 / cycle_h if cycle_h > 0 and math.isfinite(cycle_h) else 0.0
    if peak_cycles > max_cyc:
        notes.append(f"On N day {pk['n_day']} the column needs {peak_cycles:.1f} cycles but can run "
                     f"about {max_cyc:.1f} a day; the sump will grow. Use a larger or second column.")
    elif max_cyc > 0 and peak_cycles > 0.9 * max_cyc:
        notes.append(f"On N day {pk['n_day']} the column runs at {100 * peak_cycles / max_cyc:.0f} % of its "
                     f"cycle limit ({peak_cycles:.1f} of {max_cyc:.1f} a day); titre swings will cause "
                     "waits and sump build-up. Allow at least 10 % spare.")

    sump_max = cp.sump_max_L if cyc["mode"] == "batch" else p.surge.V_max_L
    if first is None:
        notes.append("The sump never collects one column load in this run; use a smaller column or a longer run.")
    elif cyc["mode"] == "batch":
        hold_h = (first["n_day"] - cp.capture_start_day) * 24.0
        if first["sump_L"] > sump_max:
            notes.append(f"Before the first cycle the sump must hold {first['sump_L']:.0f} L, "
                         f"more than its {sump_max:.0f} L capacity.")
        if hold_h > cp.max_sump_hold_h:
            notes.append(f"The first product waits about {hold_h:.0f} h in the sump (limit "
                         f"{cp.max_sump_hold_h:.0f} h); check product stability or load a partial first cycle.")
    last_day = pts[-1][0] * up.time_stretch
    if p.run_days > last_day + 0.5:
        notes.append(f"Run is longer than the titre curve ({last_day:.1f} days); titre is held at "
                     f"{pts[-1][1] * up.titre_scale:.2f} g/L after that.")

    first_buf = None
    if first and daily:
        first_buf = daily[min(int(first["n_day"]), len(daily) - 1)]["buffer_L"]
    peak_buf = max((r["buffer_L"] for r in daily), default=0.0)
    return {
        "mode": "profile", "molecule": up.molecule, "molecule_label": profile_label(p),
        "feed_L_day": _r(F * 24, 1), "permeate_L_day": _r(Fp * 24, 1),
        "tff_retentate_L_day": _r(Fr * 24, 2),
        "capture_mode": cyc["mode"], "column_CV_L": _r(cyc["CV_L"], 2),
        "load_target_g": _r(target, 1), "capture_buffer_L_per_cycle": _r(cyc["buffer_L"], 1),
        "capture_cycle_h_at_peak": _r(cycle_h, 2), "max_cycles_day_at_peak": _r(max_cyc, 1),
        "first_cycle_n_day": first["n_day"] if first else None,
        "first_cycle_sump_L": first["sump_L"] if first else None,
        "peak_day": pk["n_day"], "peak_g_day": _r(peak_g, 1),
        "peak_cycles_day": _r(peak_cycles, 2), "peak_buffer_L_day": _r(peak_buf, 1),
        "buffer_swing": _r(peak_buf / first_buf, 1) if first_buf else None,
        "buffer_L_per_g": {k: round(v, 4) for k, v in per_g.items()},
        "total_g_harvested": _r(total_g, 1),
        "ds_kg_est": _r(total_g * train_yield(p, cyc) / 1000),
        "daily": daily,
        "notes": notes or ["Design is consistent with the titre curve."],
    }


def process_design(p):
    d = steady_state_design(p) if p.upstream.mode == "mechanistic" else profile_design(p)
    d["n1"] = n1_summary(p)
    return d


def portfolio_sizing(p):
    """Run the planning view for every molecule and name the sizing cases."""
    names = list(TITRE_PROFILES) + (["custom"] if p.upstream.custom_points else [])
    rows = []
    for name in names:
        q = copy.deepcopy(p)
        q.upstream.mode = "profile"
        q.upstream.molecule = name
        d = profile_design(q)
        rows.append({"molecule": name, "label": d["molecule_label"],
                     "peak_g_day": d["peak_g_day"], "peak_day": d["peak_day"],
                     "peak_cycles_day": d["peak_cycles_day"],
                     "max_cycles_day_at_peak": d["max_cycles_day_at_peak"],
                     "first_cycle_n_day": d["first_cycle_n_day"],
                     "first_cycle_sump_L": d["first_cycle_sump_L"],
                     "peak_buffer_L_day": d["peak_buffer_L_day"],
                     "total_g": d["total_g_harvested"]})
    buf_case = max(rows, key=lambda r: r["peak_buffer_L_day"] or 0.0)
    top = buf_case["peak_buffer_L_day"] or 0.0
    for r in rows:
        r["pct_of_buffer_sizing"] = round(100 * (r["peak_buffer_L_day"] or 0.0) / top) if top > 0 else None
    sump_rows = [r for r in rows if r["first_cycle_sump_L"] is not None]
    sump_case = max(sump_rows, key=lambda r: r["first_cycle_sump_L"]) if sump_rows else None
    notes = [f"Size buffer supply and column cycling for the {buf_case['label'].lower()}: "
             f"{buf_case['peak_buffer_L_day']} L/day at {buf_case['peak_cycles_day']} cycles/day. "
             "The other molecules run inside that envelope; only their timing differs."]
    if sump_case:
        notes.append(f"Size the harvest sump for the {sump_case['label'].lower()}: "
                     f"{sump_case['first_cycle_sump_L']} L before its first cycle (N day "
                     f"{sump_case['first_cycle_n_day']}). Slow starters set sump size and hold time.")
    missing = [r["label"] for r in rows if r["first_cycle_sump_L"] is None]
    if missing:
        notes.append("Never collects one column load in this run: " + ", ".join(missing) + ".")
    return {"rows": rows, "buffer_sizing_case": buf_case["molecule"],
            "sump_sizing_case": sump_case["molecule"] if sump_case else None,
            "peak_buffer_L_day": top, "notes": notes}
