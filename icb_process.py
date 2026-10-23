"""
icb_process.py — Integrated Dynamic BioProcessing (IDB) for mAbs:
flowsheet, parameters, titre curves, solutions and design checks.
NavNevAI_mAbWasteWater · CarbonTesseract 2026

IDB (site definition)
  N-1   perfusion seed: media in = cell-free permeate out to waste.
  N     continuous media feed; outflow = cells + broth + product (in = out).
  TFF   cell-free permeate -> conditioned media surge tank (harvest sump).
        upstream.cell_retention = True (default): cells go back to the reactor and
        only the cell bleed (bleed_frac of the outflow) goes to bio-kill; product
        passage (sieve) falls with use and the TFF set is swapped at its capacity.
        upstream.cell_retention = False: all cells leave in the TFF cell concentrate
        (tff_retentate_frac of the outflow) to bio-kill.
  DSP   Case 5 flowsheet, run step by step in icb_dsp.py:
        Brij-9 in-line -> Protein A (columns 1&2, 3&4) -> Protein A vessel -> depth
        filtration -> surge tank -> Tris base + NaCl in-line -> AEX (2 columns) -> AEX
        vessel -> acetic acid in-line -> CEX (2 columns) -> CEX vessel -> virus
        filtration (2 holders) -> VRF vessel -> UF/DF -> BDS hold mixer -> 0.2 um ->
        cryovault containers.
        Options kept: low-pH VI after Protein A (vi.mode = "low_ph"), CEX before AEX
        (polish.order = "cex_aex"), 3-column PCC capture (capture.mode = "pcc").

Every solution on the flowsheet has a key in SOLUTIONS with its composition, so a
run reports litres per solution and kilograms per raw material. Step volumes and
times are typical placeholders, not site data. The flowsheet names the depth
filtration solutions but not their order; the order used here is assumed.

Biology, seed train N-3 -> N and production-reactor sizing from demand live in icb_bio.py
(parameters bio, seedtrain, demand below; titre curve "biology"; seedtrain.apply_to_n1).

Units: time h, volume L, mass g, cell density 1e6 cells/mL, rates pg/cell/h.
"""
from __future__ import annotations
import bisect
import copy
import math
from dataclasses import dataclass, field, asdict, is_dataclass
from icb_bio import bio_titre_points, apply_seed_train_to_n1

# Permeate titre (g/L) against N day, digitised from Figure 2. Replace with raw data.
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

# titre curves not digitised from Figure 2
EXTRA_PROFILES = {"biology": "Biology model (Bio_sim)"}


def profile_names(p):
    """Every titre curve the portfolio and sizing loops run: Figure 2, the biology model, any custom curve."""
    return list(TITRE_PROFILES) + list(EXTRA_PROFILES) + (["custom"] if p.upstream.custom_points else [])

# ── solutions on the Case 5 flowsheet ─────────────────────────────────
PKA_ACETIC, PKA_TRIS = 4.76, 8.07          # at 25 °C


def _acetate(mM, pH):
    """Sodium acetate + acetic acid at the stated total and pH (Henderson-Hasselbalch)."""
    base = 1.0 / (1.0 + 10 ** (PKA_ACETIC - pH))
    M = mM / 1000.0
    return {"Sodium acetate (anhydrous)": M * base * 82.03, "Acetic acid (glacial)": M * (1 - base) * 60.05}


def _tris(mM, pH):
    """Tris base titrated with HCl to the stated pH."""
    acid = 1.0 / (1.0 + 10 ** (pH - PKA_TRIS))
    M = mM / 1000.0
    return {"Tris base": M * 121.14, "Hydrochloric acid (as 100 % HCl)": M * acid * 36.46}


def _nacl(mM):
    return {"Sodium chloride": mM / 1000.0 * 58.44}


def titrant_frac(buffer_mM, pKa_buffer, pH_in, pH_out, titrant_M, titrant_pKa=None, base=True):
    """Litres of titrant per litre of pool to move one buffer from pH_in to pH_out.

    Buffer only: the protein and any other species are not included, so this is a
    starting point for the in-line addition, not a substitute for a titration.
    titrant_pKa accounts for a weak titrant (Tris base, acetic acid); leave it None
    for a strong one (NaOH, HCl).
    """
    def f(pH):
        return 1.0 / (1.0 + 10 ** (pKa_buffer - pH))
    need = (buffer_mM / 1000.0) * (f(pH_out) - f(pH_in))       # mol of strong base per L
    if titrant_pKa is None:
        eff = 1.0
    elif base:
        eff = 1.0 / (1.0 + 10 ** (pH_out - titrant_pKa))        # protonated share accepts H+
    else:
        eff = 1.0 / (1.0 + 10 ** (titrant_pKa - pH_out))        # dissociated share gives H+
    return abs(need) / max(titrant_M * eff, 1e-12)


def titrant_plan(p, capture_buffer_mM=25.0, aex_buffer_mM=50.0):
    """In-line base before AEX and acid before CEX from the pH each step has to reach."""
    return {"base_frac": titrant_frac(capture_buffer_mM, PKA_ACETIC, p.capture.elute_pH,
                                      7.4, 2.0, PKA_TRIS, base=True),
            "acid_frac": titrant_frac(aex_buffer_mM, PKA_TRIS, 7.4, p.cex.load_pH,
                                      1.0, PKA_ACETIC, base=False)}


def set_titrants(p, capture_buffer_mM=25.0, aex_buffer_mM=50.0):
    """Size the in-line base and acid from the pH each step has to reach.
    Call it after the case is applied; confirm both by titration."""
    t = titrant_plan(p, capture_buffer_mM, aex_buffer_mM)
    p.aex.base_frac = t["base_frac"]
    p.cex.acid_frac = t["acid_frac"]
    return p


def _mix(*parts):
    out = {}
    for part in parts:
        for k, v in part.items():
            out[k] = out.get(k, 0.0) + v
    return out


def _sol(label, where, g_per_L):
    return {"label": label, "where": where, "g_per_L": g_per_L}


_PRO = 115.13
# grams of each raw material per litre (approximate: pH trim of arginine/phosphate not included)
SOLUTIONS = {
    "wfi": _sol("WFI", "All units: flushes, rinses, in-line dilution", {}),
    "brij_10": _sol("10 %wt Brij-9 concentrate", "Detergent viral inactivation (stream 8)",
                    {"Brij-9 (PEG-9 lauryl ether)": 100.0}),
    "tris_74": _sol("50 mM Tris pH 7.4", "Protein A, depth filtration, AEX", _tris(50, 7.4)),
    "pa_wash1": _sol("25 mM Tris, 0.42 M arginine, 0.5 M NaCl, 0.1 % Virodex TXR-1, pH 9.0", "Protein A wash 1",
                     _mix(_tris(25, 9.0), {"L-arginine": 0.42 * 174.20}, _nacl(500), {"Virodex TXR-1": 1.0})),
    "pa_wash2": _sol("25 mM sodium acetate pH 6.5", "Protein A wash 2", _acetate(25, 6.5)),
    "pa_elute": _sol("25 mM sodium acetate pH 4.45", "Protein A elution", _acetate(25, 4.45)),
    "pa_strip": _sol("100 mM acetic acid", "Protein A strip", {"Acetic acid (glacial)": 0.1 * 60.05}),
    "naoh_1n": _sol("1 N sodium hydroxide (stock)",
                    "Protein A (0.5 N in-line), depth filtration, AEX, CEX, UF/DF (and 0.1 N storage)",
                    {"Sodium hydroxide": 40.00}),
    "storage": _sol("2 % benzyl alcohol, 100 mM sodium acetate pH 5.0", "Protein A, AEX and CEX storage",
                    _mix({"Benzyl alcohol": 20.0}, _acetate(100, 5.0))),
    "df_ac50": _sol("50 mM sodium acetate pH 5.0", "Depth filtration equilibration and chase", _acetate(50, 5.0)),
    "nacl_2m": _sol("2 M sodium chloride (stock)", "Depth filtration (0.5 M in-line), base addition before AEX",
                    _nacl(2000)),
    "df_arg": _sol("1 M arginine, 50 mM sodium phosphate pH 7", "Depth filtration",
                   {"L-arginine": 174.20, "Sodium phosphate (as NaH2PO4)": 0.05 * 119.98}),
    "tris_base_2m": _sol("2 M Tris base", "On-line base addition before AEX", {"Tris base": 2 * 121.14}),
    "tris_nacl_74": _sol("50 mM Tris, 0.5 M NaCl, pH 7.4", "AEX and CEX strip", _mix(_tris(50, 7.4), _nacl(500))),
    "acetic_1m": _sol("1 M acetic acid", "On-line acid addition before CEX", {"Acetic acid (glacial)": 60.05}),
    "cex_eq": _sol("50 mM sodium acetate, 120 mM NaCl, pH 5.5", "CEX equilibration and wash",
                   _mix(_acetate(50, 5.5), _nacl(120))),
    "cex_elute": _sol("50 mM sodium acetate, 280 mM NaCl, pH 6.0", "CEX elution, virus filtration",
                      _mix(_acetate(50, 6.0), _nacl(280))),
    "uf_df1": _sol("15 mM acetate, 130 mM proline, pH 5.1", "UF/DF equilibration and diafiltration",
                   _mix(_acetate(15, 5.1), {"L-proline": 0.130 * _PRO})),
    "uf_df2": _sol("15 mM acetate, 113 mM proline, pH 5.30", "UF/DF recovery flush",
                   _mix(_acetate(15, 5.30), {"L-proline": 0.113 * _PRO})),
    "uf_spike": _sol("25 mM acetate, 113 mM proline, 1.17 M sucrose, 6 mg/mL poloxamer 188, pH 5.30",
                     "Formulation spike",
                     _mix(_acetate(25, 5.30), {"L-proline": 0.113 * _PRO, "Sucrose": 1.17 * 342.30,
                                               "Poloxamer 188": 6.0})),
    "vi_acid": _sol("Low-pH VI acid", "Low-pH viral inactivation (legacy mode)", {}),
    "vi_base": _sol("Low-pH VI neutralising base", "Low-pH viral inactivation (legacy mode)", {}),
}
# solutions made in-line from a stock: (stock, share) per litre delivered
DILUTIONS = {
    "naoh_05n": (("naoh_1n", 0.5), ("wfi", 0.5)),
    "naoh_01n": (("naoh_1n", 0.1), ("wfi", 0.9)),
    "nacl_05m": (("nacl_2m", 0.25), ("wfi", 0.75)),
}
SOLUTION_ALIASES = {
    "naoh_05n": "0.5 N NaOH (1 N diluted in-line with WFI)",
    "naoh_01n": "0.1 N NaOH storage (1 N diluted with WFI)",
    "nacl_05m": "0.5 M NaCl (2 M diluted in-line with WFI)",
}


# ── parameters ────────────────────────────────────────────────────────
@dataclass
class N1Params:
    enabled: bool = True
    V_L: float = 15.0
    VVD: float = 2.0
    days: float = 5.0                # fixed N-1 duration (used only when target_seed_vcd = 0)
    target_seed_vcd: float = 2.0     # inoculate N when the transfer gives this VCD (1e6/mL); 0 = fixed days
    max_days: float = 8.0            # transfer anyway after this long if the target is not reached
    X0: float = 0.5
    mu_h: float = 0.025
    X_max: float = 40.0
    transfer_L: float = 15.0


@dataclass
class UpstreamParams:
    mode: str = "profile"            # "profile" (IDB) or "mechanistic" (legacy ATF model)
    molecule: str = "purple"
    custom_points: list = field(default_factory=list)
    VVD: float = 1.2
    titre_scale: float = 1.0
    time_stretch: float = 1.0
    variation_sigma: float = 0.04
    variation_tau_h: float = 48.0
    # --- cell handling at the TFF ---
    cell_retention: bool = True       # True: cells go back to the reactor, only the bleed leaves
    bleed_frac: float = 0.02          # share of the reactor outflow taken as cell bleed
    tff_retentate_frac: float = 0.02  # used only when cell_retention is False
    sieve0: float = 0.98              # product passage through a fresh TFF set
    sieve_loss_per_day: float = 0.006
    tff_area_m2: float = 7.2
    tff_capacity_Lm2: float = 2000.0  # permeate per m² before a set is swapped
    tff_swap_flush_L: float = 240.0   # media and caustic flush per swap


def cell_out_frac(up):
    """Share of the N outflow that leaves with the cells (bleed or TFF cell concentrate)."""
    return up.bleed_frac if up.cell_retention else up.tff_retentate_frac


@dataclass
class BioreactorParams:
    V_L: float = 83.0
    X0: float = 0.5
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
    mode: str = "batch"              # "batch" (Case 5 recipe) or "pcc"
    trains: int = 2                  # columns 1&2 and 3&4: one loads while the other recovers
    column_capacity_g: float = 100.0 # per train (DBC x CV)
    DBC_gL: float = 40.0             # -> CV = 2.5 L per train
    load_util: float = 0.90
    load_flow_CVh: float = 15.0      # load incl. in-line detergent
    postload_CV: float = 3.0         # 50 mM Tris pH 7.4
    wash1_CV: float = 3.0            # arginine / NaCl / Virodex pH 9.0
    wash2_CV: float = 3.0            # 25 mM acetate pH 6.5
    elute_CV: float = 3.0            # 25 mM acetate pH 4.45
    eluate_CV: float = 2.0           # collected peak
    strip_CV: float = 2.0            # 100 mM acetic acid
    cip_CV: float = 3.0              # 0.5 N NaOH (1 N diluted in-line)
    cip_contact_h: float = 0.25
    reequil_CV: float = 5.0
    recovery_flow_CVh: float = 10.0
    startup_wfi_CV: float = 3.0
    startup_naoh_CV: float = 3.0
    equil_CV: float = 3.0            # pre-use equilibration
    storage_wfi_CV: float = 2.0
    storage_CV: float = 3.0          # 2 % benzyl alcohol
    capture_yield: float = 0.97
    sump_max_L: float = 600.0
    max_sump_hold_h: float = 96.0
    capture_start_day: float = 0.0
    resin_cycles_max: int = 200
    elute_pH: float = 4.45
    sol_postload: str = "tris_74"
    sol_wash1: str = "pa_wash1"
    sol_wash2: str = "pa_wash2"
    sol_elute: str = "pa_elute"
    sol_strip: str = "pa_strip"
    sol_cip: str = "naoh_05n"
    sol_equil: str = "tris_74"


@dataclass
class SurgeParams:
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
    mode: str = "detergent"          # "detergent" (Brij-9 before Protein A, Case 5) or "low_ph"
    det_target_pct: float = 0.5      # Brij-9 in the load (placeholder)
    det_conc_pct: float = 10.0       # concentrate (stream 8)
    det_min_contact_min: float = 1.0
    det_line_hold_min: float = 2.0
    eluate_surge_max_L: float = 40.0 # Protein A vessel SUM capacity
    heel_L: float = 0.5              # --- low-pH mode only below ---
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
class DepthParams:                   # new filter each batch; solution order assumed
    area_m2: float = 0.2
    capacity_Lm2: float = 150.0
    flux_LMH: float = 150.0
    flush_flux_LMH: float = 300.0
    wfi_flush_Lm2: float = 50.0
    equil_Lm2: float = 20.0          # 50 mM acetate pH 5.0
    chase_Lm2: float = 10.0          # same, kept with product
    salt_Lm2: float = 10.0           # 0.5 M NaCl (2 M in-line)
    arginine_Lm2: float = 10.0
    naoh_Lm2: float = 10.0
    naoh_contact_h: float = 0.25
    tris_Lm2: float = 10.0
    trigger_elutions: float = 3.0    # flowsheet: minimum of 3 elution volumes
    depth_yield: float = 0.98


@dataclass
class AEXParams:
    columns: int = 2
    CV_L: float = 2.0
    load_gL: float = 150.0
    flow_CVh: float = 10.0
    base_frac: float = 0.04          # 2 M Tris base, L per L load
    nacl_frac: float = 0.03          # 2 M NaCl, L per L load
    chase_CV: float = 3.0
    strip_CV: float = 3.0
    cip_CV: float = 3.0
    cip_contact_h: float = 0.5
    reequil_CV: float = 5.0
    startup_wfi_CV: float = 3.0
    startup_naoh_CV: float = 3.0
    equil_CV: float = 5.0
    storage_wfi_CV: float = 2.0
    storage_CV: float = 3.0
    min_load_frac: float = 0.5


@dataclass
class CEXParams:
    columns: int = 2
    CV_L: float = 3.0
    DBC_gL: float = 50.0
    flow_CVh: float = 10.0
    acid_frac: float = 0.03          # 1 M acetic acid, L per L load
    wash_CV: float = 3.0
    elute_CV: float = 5.0
    eluate_CV: float = 3.0
    strip_CV: float = 3.0
    cip_CV: float = 3.0
    cip_contact_h: float = 0.5
    reequil_CV: float = 5.0
    startup_wfi_CV: float = 3.0
    startup_naoh_CV: float = 3.0
    equil_CV: float = 5.0
    storage_wfi_CV: float = 2.0
    storage_CV: float = 3.0
    min_load_frac: float = 0.5
    load_util: float = 0.90          # share of binding capacity actually loaded
    load_pH: float = 5.5
    sol_eq: str = "cex_eq"
    sol_elute: str = "cex_elute"


@dataclass
class PolishParams:
    order: str = "aex_cex"
    aex_yield: float = 0.98
    cex_yield: float = 0.92
    vf_yield: float = 0.99
    vf_capacity_L: float = 300.0
    vf_flush_L: float = 2.0
    prefilter_per_day: float = 0.0   # 0.2 um prefilters changed per day before the virus filter


@dataclass
class VRFParams:
    flow_Lh: float = 20.0
    equil_L: float = 2.0
    chase_L: float = 1.0
    trigger_elutions: float = 3.0
    sol_equil: str = "cex_elute"
    integrity_pre_h: float = 2.0
    integrity_post_h: float = 2.0
    max_hold_h: float = 24.0         # most material one filter may see


@dataclass
class SPTFFParams:                   # UF/DF and formulation (name kept for compatibility)
    target_gL: float = 50.0
    DV: float = 7.0
    sptff_yield: float = 0.98
    excipient_frac: float = 0.10     # formulation spike, share of final volume
    uf1_gL: float = 25.0
    perm_flux_Lh: float = 15.0
    transfer_Lh: float = 30.0
    flush_Lh: float = 30.0
    recover_L: float = 2.0
    min_V_L: float = 1.0
    wfi_flush_L: float = 20.0
    equil_L: float = 5.0
    cip_L: float = 10.0
    cip_contact_h: float = 1.0
    storage_L: float = 10.0
    min_lot_g: float = 1.0
    sol_df: str = "uf_df1"
    sol_recover: str = "uf_df2"
    sol_spike: str = "uf_spike"


@dataclass
class BDSParams:
    mix_h: float = 0.5
    filter_flush_L: float = 1.0
    fill_Lh: float = 10.0
    container_L: float = 2.0


@dataclass
class DSPParams:
    tail_h: float = 72.0             # downstream run-out after the N run ends


@dataclass
class PATParams:
    vcd_noise: float = 0.03
    titre_noise: float = 0.05
    ph_noise: float = 0.02


@dataclass
class ControlParams:
    CSPR_pL_cell_day: float = 30.0
    bleed_kc: float = 0.07
    bleed_ti_h: float = 24.0
    bleed_enable_frac: float = 0.95
    surge_tau_h: float = 4.0
    eluate_tau_h: float = 3.0
    acid_kc: float = 0.05
    acid_ti_h: float = 6.0
    viability_alarm: float = 0.80
    glucose_alarm_gL: float = 0.3
    titre_drop_alarm: float = 0.15
    # --- NevAI in operation (icb_nevai_ops.py) ---
    nevai_act: bool = False           # let forecasts start an early partial Protein A load
    nevai_partial_min_frac: float = 0.3   # smallest partial load, share of a full load
    nevai_hold_margin_h: float = 6.0  # act when product would pass the hold limit minus this
    nevai_watch_lead_h: float = 6.0   # compare titre with the forecast issued this long before
    nevai_watch_window_h: float = 6.0
    nevai_watch_baseline_h: float = 48.0
    nevai_watch_threshold: float = 0.12   # flag a forecast-error shift above this share
    buffer_prep_Lh: float = 0.0       # buffer preparation capacity, L/h (0 = peak day / 20 h)
    buffer_hold_days: float = 1.0     # prepared hold target, days of forecast use


@dataclass
class SizingParams:                  # basis for icb_sizing.py: typical values, replace with site standards
    util_target: float = 0.80        # max share of the peak day a unit (or its load pump) may be busy
    vessel_margin: float = 1.20      # on the peak working volume
    vessel_fill: float = 0.80        # max working fill of a vessel
    pump_margin: float = 1.25
    area_margin: float = 1.25        # filters and membranes (fouling)
    min_scale: float = 0.25          # never shrink a unit below this share of its configured size
    tote_days: float = 1.0           # supply held per tote, days of peak use
    media_days: float = 1.0
    pa_bed_cm: float = 20.0
    aex_bed_cm: float = 15.0
    cex_bed_cm: float = 20.0
    max_u_cmh: float = 600.0         # max linear velocity
    tff_flux_LMH: float = 20.0       # cell-removal TFF permeate flux
    vf_flux_LMH: float = 60.0
    vf_Lm2: float = 200.0            # virus filter throughput before change
    ufdf_flux_LMH: float = 40.0
    ufdf_load_gm2: float = 500.0
    ufdf_crossflow_LMH: float = 300.0
    filter_flux_LMH: float = 1000.0  # in-line 0.2 um filters


@dataclass
class BioParams:                     # Bio_inputs: typical CHO values until fitted to the cell line
    mu_max_d: float = 0.8
    tau_lag_d: float = 0.6
    K_glc_gL: float = 0.5
    K_gln_mM: float = 0.3
    KI_lac_gL: float = 10.0
    KI_amm_mM: float = 15.0
    osm_crit: float = 340.0
    K_osm: float = 110.0
    shift_day: float = 7.0           # beyond run_days = 37 °C throughout
    T_before_C: float = 36.8         # labels only
    T_after_C: float = 33.0
    f_mu_shift: float = 0.45
    f_qp_shift: float = 1.35
    f_maint_shift: float = 0.75
    f_kd_shift: float = 0.7
    kd_min_d: float = 0.01
    kd_amm_max_d: float = 0.2
    Kd_amm_mM: float = 10.0
    n_amm: float = 3.0
    kd_osm_max_d: float = 0.25
    Kd_osm: float = 560.0
    n_osm: float = 8.0
    kd_starve_max_d: float = 0.5
    K_starve_gL: float = 0.1
    k_lysis_d: float = 0.05
    glc_per_cell: float = 0.3        # 1/Y_X/Glc, g/L per 1e6 cells/mL
    m_glc: float = 0.06
    gln_per_cell: float = 0.4        # 1/Y_X/Gln, mM per 1e6 cells/mL
    m_gln: float = 0.015
    Y_lac_glc: float = 0.9
    q_lac_max: float = 0.06
    K_lac_gL: float = 0.5
    Y_amm_gln: float = 0.7
    k_deg_gln_d: float = 0.07
    S0_mOsm: float = 262.0
    S_feed_mOsm: float = 900.0
    qp_osm_rise: float = 0.35
    qP_pg_cell_d: float = 22.0
    X0: float = 0.8                  # also the N seed density in the seed train
    glc0_gL: float = 6.0
    gln0_mM: float = 4.0
    lac0_gL: float = 0.1
    amm0_mM: float = 0.5
    V0_L: float = 1400.0
    feed_start_d: float = 2.0
    CSFR_pL_cell_d: float = 2.0
    feed_glc_gL: float = 50.0
    feed_gln_mM: float = 50.0
    dt_d: float = 0.1
    run_days: float = 14.0
    target_titre_gL: float = 3.0
    report_days: float = 0.5
    titre_link: str = "fedbatch"     # titre curve "biology": "fedbatch" (Bio_sim titre) or "perfusion" (qP·Xv/VVD)


@dataclass
class SeedTrainParams:               # Seed_train: N-3, N-2, N-1 and N back to the vial
    n_basis: str = "demand"          # N starting volume: "demand", "twin" (bioreactor.V_L), "manual"
    manual_start_L: float = 1400.0
    manual_vessel_L: float = 2000.0
    duration_basis: str = "biology"  # "biology" (Seed_sim) or "rule" (fixed mu_rule_d)
    density_basis: str = "biology"   # "biology" (net growth >= threshold) or "rule" (xt_rule_*)
    growth_threshold_pct: float = 80.0
    mu_rule_d: float = 0.65
    vial_cells_1e6: float = 30.0
    qO2_pmol_cell_h: float = 0.15
    O2_sat_air_mM: float = 0.2
    DO_pct: float = 40.0
    sim_days: float = 10.0
    type_n3: str = "Rocking bag"
    type_n2: str = "Stirred tank"
    type_n1: str = "Stirred tank"
    max_fill_n3: float = 0.5
    max_fill_n2: float = 1.0
    max_fill_n1: float = 1.0
    max_fill_n: float = 1.0
    min_fill_n3: float = 0.1
    min_fill_n2: float = 0.2
    min_fill_n1: float = 0.2
    min_fill_n: float = 0.2
    seed_n3: float = 0.5
    seed_n2: float = 0.5
    seed_n1: float = 0.5
    xt_rule_n3: float = 3.0
    xt_rule_n2: float = 4.0
    xt_rule_n1: float = 6.0
    split_min: float = 4.0
    split_max: float = 12.0
    stage_min_d: float = 2.0
    stage_max_d: float = 5.0
    viability_min_pct: float = 95.0
    apply_to_n1: bool = False        # run the twin's N-1 as this batch N-1 stage (train sized for the twin N)
    vessel_sizes_L: list = field(default_factory=lambda: [2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 3000,
                                                          5000, 10000, 12000, 15000, 20000, 25000])


@dataclass
class DemandParams:                  # Demand: production reactor N from annual demand
    annual_kg: float = 400.0
    dsp_yield: float = 0.75
    batch_success: float = 0.95
    turnaround_d: float = 3.0
    operating_days: float = 330.0
    n_reactors: int = 1


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
    depth: DepthParams = field(default_factory=DepthParams)
    aex: AEXParams = field(default_factory=AEXParams)
    cex: CEXParams = field(default_factory=CEXParams)
    polish: PolishParams = field(default_factory=PolishParams)
    vrf: VRFParams = field(default_factory=VRFParams)
    sptff: SPTFFParams = field(default_factory=SPTFFParams)
    bds: BDSParams = field(default_factory=BDSParams)
    dsp: DSPParams = field(default_factory=DSPParams)
    pat: PATParams = field(default_factory=PATParams)
    control: ControlParams = field(default_factory=ControlParams)
    sizing: SizingParams = field(default_factory=SizingParams)
    bio: BioParams = field(default_factory=BioParams)
    seedtrain: SeedTrainParams = field(default_factory=SeedTrainParams)
    demand: DemandParams = field(default_factory=DemandParams)
    run_days: float = 14.0
    lot_window_h: float = 24.0       # UF/DF collects this long per lot
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
WASTE_STREAMS = {
    "n1_permeate":         ("USP", "direct",         "volume", "N-1 permeate: spent media removed to intensify the seed"),
    "tff_retentate":       ("USP", "bio-kill",       "volume", "TFF cell concentrate from the N outflow"),
    "cell_bleed":          ("USP", "bio-kill",       "volume", "Cell bleed (TFF cell retention or legacy ATF model)"),
    "atf_flush":           ("USP", "bio-kill",       "event",  "ATF filter swap flush (legacy ATF model)"),
    "tff_swap_flush":      ("USP", "bio-kill",       "event",  "TFF set swap: media push, caustic and rinse"),
    "ramp_harvest":        ("USP", "direct",         "volume", "Permeate to drain before capture starts (stream 2)"),
    "sump_overflow":       ("USP", "direct",         "event",  "Surge tank overflow (cell-free)"),
    "capture_flowthrough": ("DSP", "direct",         "volume", "Protein A load flow-through with Brij-9: high organic load"),
    "capture_buffers":     ("DSP", "direct",         "mass",   "Protein A post-load wash, washes 1-2, elution pre/post-peak, re-equilibration"),
    "capture_strip_cip":   ("DSP", "neutralisation", "mass",   "Protein A acetic acid strip and 0.5 N NaOH CIP"),
    "unit_flushes":        ("DSP", "direct",         "event",  "Pre-use and end-of-campaign WFI flushes and equilibration of columns"),
    "sanitisation":        ("DSP", "neutralisation", "event",  "Pre-use caustic sanitisation of columns"),
    "vi_flush":            ("DSP", "neutralisation", "event",  "Low-pH VI tube start-up buffer (low-pH mode)"),
    "vi_divert":           ("DSP", "neutralisation", "event",  "Out-of-spec VI material diverted (low-pH mode)"),
    "depth_flush":         ("DSP", "direct",         "mass",   "Depth filter WFI flush and equilibration (new filter per batch)"),
    "depth_post":          ("DSP", "neutralisation", "mass",   "Depth filtration salt, arginine, caustic and Tris flushes (order assumed)"),
    "aex_buffers":         ("DSP", "direct",         "mass",   "AEX strip and re-equilibration"),
    "aex_cip":             ("DSP", "neutralisation", "mass",   "AEX 1 N NaOH CIP"),
    "cex_flowthrough":     ("DSP", "direct",         "mass",   "CEX load flow-through (acidified load)"),
    "cex_buffers":         ("DSP", "direct",         "mass",   "CEX wash, elution pre/post-peak, strip and re-equilibration"),
    "cex_cip":             ("DSP", "neutralisation", "mass",   "CEX 1 N NaOH CIP"),
    "vf_flush":            ("DSP", "direct",         "mass",   "Virus filter WFI flush and equilibration (new filter)"),
    "sptff_permeate":      ("DSP", "direct",         "mass",   "UF/DF permeate incl. diafiltration buffer"),
    "ufdf_cip":            ("DSP", "neutralisation", "mass",   "UF/DF membrane flushes, equilibration and 1 N NaOH CIP"),
    "bds_flush":           ("DSP", "direct",         "mass",   "0.2 um fill filter WFI flush"),
}
BUFFER_OPS = ("capture", "vi", "depth", "adjust", "aex", "cex", "vf", "sptff", "bds")


def _r(x, n=3):
    return round(x, n) if isinstance(x, (int, float)) and math.isfinite(x) else None


# ── titre curves ──────────────────────────────────────────────────────
def profile_points(p):
    up = p.upstream
    if up.molecule == "biology":
        return bio_titre_points(p)
    if up.molecule == "custom" and up.custom_points:
        try:
            pts = sorted([float(a), float(b)] for a, b in up.custom_points)
            if pts:
                return pts
        except (TypeError, ValueError):
            pass
    return TITRE_PROFILES.get(up.molecule, TITRE_PROFILES["purple"])["points"]


def profile_label(p):
    m = p.upstream.molecule
    if m == "custom" and p.upstream.custom_points:
        return "Custom curve"
    if m == "biology":
        return EXTRA_PROFILES["biology"] + (" — perfusion qP·Xv/VVD" if p.bio.titre_link == "perfusion"
                                            else " — fed-batch titre")
    return TITRE_PROFILES.get(m, TITRE_PROFILES["purple"])["label"]


def titre_at(p, n_day, pts=None):
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
def det_frac(p):
    """Litres of Brij-9 concentrate per litre of permeate loaded (0 in low-pH mode)."""
    vi = p.vi
    if vi.mode == "low_ph" or vi.det_target_pct <= 0:
        return 0.0
    c = max(vi.det_conc_pct, vi.det_target_pct + 1e-6)
    return vi.det_target_pct / (c - vi.det_target_pct)


def capture_cycle(p):
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
                "yield": pc.capture_yield, "resin_cycles_max": pc.resin_cycles_max, "trains": 1}
    cp = p.capture
    cv = cp.column_capacity_g / max(cp.DBC_gL, 1e-9)
    cvs = cp.postload_CV + cp.wash1_CV + cp.wash2_CV + cp.elute_CV + cp.strip_CV + cp.cip_CV + cp.reequil_CV
    buf = cvs * cv
    strip_cip = (cp.strip_CV + cp.cip_CV) * cv
    return {"mode": "batch", "CV_L": cv, "target_g": cp.column_capacity_g * cp.load_util,
            "capacity_g": cp.column_capacity_g, "buffer_L": buf,
            "eluate_L": cp.eluate_CV * cv, "strip_cip_L": strip_cip,
            "other_buffer_waste_L": buf - cp.eluate_CV * cv - strip_cip,
            "load_flow_Lh": max(cp.load_flow_CVh * cv, 1e-6),
            "recovery_h": cvs / max(cp.recovery_flow_CVh, 1e-9) + cp.cip_contact_h,
            "yield": cp.capture_yield, "resin_cycles_max": cp.resin_cycles_max,
            "trains": max(1, int(cp.trains))}


def train_yield(p, cyc=None):
    cyc = cyc or capture_cycle(p)
    po = p.polish
    y = cyc["yield"] * p.depth.depth_yield * po.aex_yield * po.cex_yield * po.vf_yield * p.sptff.sptff_yield
    return y * (p.vi.vi_yield if p.vi.mode == "low_ph" else 1.0)


def n1_plan(p):
    """When the N-1 seed is transferred to N and the VCD it delivers."""
    n1, br = p.n1, p.bioreactor
    if not n1.enabled:
        return {"enabled": False, "hours": 0}
    share = min(n1.transfer_L, n1.V_L) / max(br.V_L, 1e-9)
    if n1.target_seed_vcd > 0:
        hmax = max(0, int(round(n1.max_days * 24)))
        need = n1.target_seed_vcd / max(share, 1e-12)        # N-1 density needed at transfer
        if need <= n1.X0:
            h, reason = 0, "target reached"
        elif n1.mu_h <= 0:
            h, reason = hmax, "N-1 does not grow (mu_h <= 0); transferred at max_days"
        elif need > n1.X_max:
            h = min(hmax, int(math.ceil(math.log(n1.X_max / n1.X0) / n1.mu_h - 1e-9)))
            reason = "target needs more than the N-1 maximum density; transferred at maximum"
        else:
            h = int(math.ceil(math.log(need / n1.X0) / n1.mu_h - 1e-9))
            reason = "target reached"
            if h > hmax:
                h, reason = hmax, "max_days reached before the target"
    else:
        h, reason = max(0, int(round(n1.days * 24))), "fixed N-1 duration"
    X = min(n1.X0 * math.exp(n1.mu_h * h), n1.X_max)
    return {"enabled": True, "hours": h, "X_end": X, "seed_vcd": X * share, "reason": reason}


def n1_summary(p):
    n1, br = p.n1, p.bioreactor
    pl = n1_plan(p)
    if not pl["enabled"]:
        return {"enabled": False, "hours": 0, "days": 0.0, "n_seed_vcd": round(br.X0, 3),
                "waste_L": 0.0, "media_L": 0.0, "reason": "N-1 off"}
    perf = n1.VVD * n1.V_L / 24.0 * pl["hours"]
    media = perf if n1.VVD > 0 else n1.V_L          # batch N-1 (seed-train mode): one fill of medium
    return {"enabled": True, "hours": pl["hours"], "days": round(pl["hours"] / 24.0, 2),
            "V_L": n1.V_L, "VVD": n1.VVD, "final_vcd": round(pl["X_end"], 2),
            "n_seed_vcd": round(pl["seed_vcd"], 3), "target_seed_vcd": n1.target_seed_vcd or None,
            "reason": pl["reason"], "media_L": round(media, 1), "waste_L": round(perf, 1),
            "batch": n1.VVD <= 0}


def _nominal_lot_g(p):
    """Grams entering capture per lot window, averaged over the run (planning only)."""
    up, br = p.upstream, p.bioreactor
    if up.mode == "mechanistic":
        return 500.0
    Fp = up.VVD * br.V_L / 24.0 * (1 - cell_out_frac(up))
    pts = profile_points(p)
    n = max(1, int(p.run_days * 24))
    g = sum(Fp * titre_at(p, k / 24.0, pts) for k in range(n))
    return max(g / n * p.lot_window_h, 1.0)


def buffer_L_per_g(p, lot_g=None):
    """Planning litres of solution (incl. WFI) per gram loaded on Protein A, by unit.
    Brij-9 follows the load volume, so it is returned separately per litre of load."""
    vi, dp, ae, ce, vr, po, sp, bd = p.vi, p.depth, p.aex, p.cex, p.vrf, p.polish, p.sptff, p.bds
    cyc = capture_cycle(p)
    tg = max(cyc["target_g"], 1e-9)
    d = dict.fromkeys(BUFFER_OPS, 0.0)
    d["capture"] = cyc["buffer_L"] / tg
    g = cyc["yield"]
    v = cyc["eluate_L"] / tg
    if vi.mode == "low_ph":
        d["vi"] = v * vi.acid_bias + v * (1 + vi.acid_bias) * vi.base_frac
        v *= (1 + vi.acid_bias) * (1 + vi.base_frac)
        g *= vi.vi_yield
    batch = max(dp.trigger_elutions, 1.0) * tg
    A = dp.area_m2
    d["depth"] = A * (dp.wfi_flush_Lm2 + dp.equil_Lm2 + dp.chase_Lm2 + dp.salt_Lm2
                      + dp.arginine_Lm2 + dp.naoh_Lm2 + dp.tris_Lm2) / batch
    v += A * dp.chase_Lm2 / batch
    g *= dp.depth_yield
    for unit in (("cex", "aex") if po.order == "cex_aex" else ("aex", "cex")):
        if unit == "aex":
            cyc_g = max(min(ae.CV_L * ae.load_gL, batch * g), 1e-9)
            d["aex"] += ae.CV_L * (ae.chase_CV + ae.strip_CV + ae.cip_CV + ae.reequil_CV) * g / cyc_g
            d["adjust"] += v * (ae.base_frac + ae.nacl_frac)
            v = v * (1 + ae.base_frac + ae.nacl_frac) + ae.CV_L * ae.chase_CV * g / cyc_g
            g *= po.aex_yield
        else:
            cyc_g = max(min(ce.CV_L * ce.DBC_gL * ce.load_util, batch * g), 1e-9)
            d["cex"] += ce.CV_L * (ce.wash_CV + ce.elute_CV + ce.strip_CV + ce.cip_CV + ce.reequil_CV) * g / cyc_g
            d["adjust"] += v * ce.acid_frac
            v = ce.CV_L * ce.eluate_CV * g / cyc_g
            g *= po.cex_yield
    vrf_batch = max(vr.trigger_elutions, 1.0) * ce.CV_L * ce.eluate_CV
    d["vf"] = (v * (po.vf_flush_L + vr.equil_L) / max(po.vf_capacity_L, 1e-9)
               + v * vr.chase_L / max(vrf_batch, 1e-9))
    g *= po.vf_yield
    lot = max(lot_g if lot_g else _nominal_lot_g(p), 1e-9)
    g_ds = g * sp.sptff_yield
    d["sptff"] = (g * sp.DV / max(sp.uf1_gL, 1e-9) + g_ds * sp.excipient_frac / max(sp.target_gL, 1e-9)
                  + (2 * sp.wfi_flush_L + sp.equil_L + sp.cip_L + sp.storage_L + sp.recover_L) / lot)
    d["bds"] = bd.filter_flush_L / lot
    d["total"] = sum(d[k] for k in BUFFER_OPS)
    d["ds_g_per_g"] = g_ds
    d["detergent_L_per_L_load"] = det_frac(p)
    return d


# ── design checks ─────────────────────────────────────────────────────
def steady_state_design(p):
    """Legacy ATF perfusion model: algebraic steady state at the VCD setpoint."""
    br, vi, ct = p.bioreactor, p.vi, p.control
    cyc = capture_cycle(p)
    dose = det_frac(p)
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
        if (pc.CV_L / Fh * 60.0 if Fh > 0 else math.inf) < pc.min_residence_min:
            notes.append("Load residence time per column is below the minimum; use larger columns.")
        if pc.load_factor > pc.sat_factor:
            notes.append("load_factor above sat_factor: product will be lost in breakthrough.")
    else:
        load_h = target / titre * (1 + dose) / cyc["load_flow_Lh"] if titre > 0 else math.inf
        lim = max(load_h, (load_h + cyc["recovery_h"]) / cyc["trains"])
        if lim > t_fill:
            notes.append(f"Protein A needs {lim:.1f} h per load across {cyc['trains']} train(s) but a load "
                         f"arrives every {t_fill:.1f} h; the sump will keep growing.")
    vi_res = None
    if vi.mode == "low_ph":
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
        "vi_residence_h": _r(vi_res, 2) if vi_res is not None else None,
        "ds_kg_day": _r(ds_g_h * 24 / 1000), "notes": notes,
    }


def profile_design(p):
    """IDB planning view from the titre curve (noise-free): daily grams, cycles and solutions."""
    up, br, cp, vi, dp = p.upstream, p.bioreactor, p.capture, p.vi, p.depth
    cyc = capture_cycle(p)
    per_g = buffer_L_per_g(p)
    dose = det_frac(p)
    pts = profile_points(p)
    target = cyc["target_g"]
    F = up.VVD * br.V_L / 24.0
    Fr = F * cell_out_frac(up)
    Fp = F - Fr
    notes = []

    daily, total_g = [], 0.0
    sump_g = sump_L = 0.0
    first = None
    for d in range(int(math.ceil(p.run_days))):
        g_d = L_d = 0.0
        for k in range(24):
            t = d + k / 24.0
            if t >= p.run_days:
                break
            if t < cp.capture_start_day:
                continue
            g = Fp * titre_at(p, t + 0.5 / 24.0, pts)
            g_d += g
            L_d += Fp
            if first is None:
                sump_g += g
                sump_L += Fp
                if sump_g >= target:
                    first = {"n_day": round(t + 1.0 / 24.0, 2), "sump_L": round(sump_L, 1)}
        total_g += g_d
        det_L = L_d * dose
        daily.append({"n_day": d + 1, "titre_gL": round(titre_at(p, d + 0.5, pts), 3),
                      "g": round(g_d, 1), "cycles": round(g_d / target, 2),
                      "capture_buffer_L": round(g_d * per_g["capture"], 1),
                      "detergent_L": round(det_L, 1),
                      "buffer_L": round(g_d * per_g["total"] + det_L, 1)})

    pk = max(daily, key=lambda r: r["g"]) if daily else {"g": 0.0, "n_day": None}
    peak_g = pk["g"]
    C_pk = peak_g / (Fp * 24.0) if Fp > 0 else 0.0
    if cyc["mode"] == "batch":
        net_g = target
        load_h = (target / C_pk) * (1 + dose) / cyc["load_flow_Lh"] if C_pk > 0 else math.inf
        cycle_h = load_h + cyc["recovery_h"]
        lim_h = max(load_h, cycle_h / cyc["trains"])
        if 1.0 - cp.load_util < 2.0 * p.pat.titre_noise - 1e-9:
            notes.append(f"Load margin ({100 * (1 - cp.load_util):.0f} %) is less than twice the titre "
                         f"measurement noise ({100 * p.pat.titre_noise:.0f} %); expect occasional overloads.")
    else:
        pc = p.pcc
        net_g = target - pc.bt_frac * max(0.0, target - pc.DBC_gL * pc.CV_L)
        cycle_h = cyc["recovery_h"] / max(pc.n_columns - 2, 1)
        lim_h = cycle_h
        if pc.n_columns < 3:
            notes.append("PCC needs at least 3 columns for uninterrupted loading.")
        if pc.load_factor > pc.sat_factor:
            notes.append("load_factor above sat_factor: product will be lost in breakthrough.")
        if Fp > cyc["load_flow_Lh"]:
            notes.append("Permeate flow exceeds the PCC load pump limit.")
    peak_cycles = peak_g / max(net_g, 1e-9)
    max_cyc = 24.0 / lim_h if lim_h > 0 and math.isfinite(lim_h) else 0.0
    if peak_cycles > max_cyc:
        notes.append(f"On N day {pk['n_day']} Protein A needs {peak_cycles:.1f} cycles but can run about "
                     f"{max_cyc:.1f} a day; the sump will grow. Use larger columns or another train.")
    elif max_cyc > 0 and peak_cycles > 0.9 * max_cyc:
        notes.append(f"On N day {pk['n_day']} Protein A runs at {100 * peak_cycles / max_cyc:.0f} % of its "
                     f"cycle limit ({peak_cycles:.1f} of {max_cyc:.1f} a day). Allow at least 10 % spare.")
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
                         f"{cp.max_sump_hold_h:.0f} h); check stability or load a partial first cycle.")
    last_day = pts[-1][0] * up.time_stretch
    if p.run_days > last_day + 0.5:
        notes.append(f"Run is longer than the titre curve ({last_day:.1f} days); titre is held at "
                     f"{pts[-1][1] * up.titre_scale:.2f} g/L after that.")
    if vi.mode != "low_ph":
        res_min = cyc["CV_L"] / max(cyc["load_flow_Lh"], 1e-9) * 60.0
        if vi.det_line_hold_min + res_min < vi.det_min_contact_min:
            notes.append(f"Brij-9 contact before and on the column ({vi.det_line_hold_min + res_min:.1f} min) is "
                         f"shorter than the {vi.det_min_contact_min:.0f} min required; add a hold loop.")
    pav_need = max(dp.trigger_elutions, 1.0) * cyc["eluate_L"]
    if pav_need > vi.eluate_surge_max_L:
        notes.append(f"The Protein A vessel must hold {pav_need:.0f} L ({dp.trigger_elutions:g} eluates) before "
                     f"depth filtration starts, more than its {vi.eluate_surge_max_L:.0f} L.")
    if pav_need / max(dp.area_m2, 1e-9) > dp.capacity_Lm2:
        notes.append(f"Depth filter loading {pav_need / max(dp.area_m2, 1e-9):.0f} L/m² exceeds its "
                     f"{dp.capacity_Lm2:.0f} L/m² capacity; use a larger filter.")
    mu_need = (up.bleed_frac if up.cell_retention else 1.0) * up.VVD / 24.0
    if mu_need > br.mu_max_h:
        notes.append(
            f"Cells must grow at {mu_need:.4f} /h to hold density, above mu_max_h "
            f"({br.mu_max_h:.4f} /h): the culture would wash out. "
            + ("Lower the bleed or the perfusion rate." if up.cell_retention else
               "Without cell retention the whole outflow carries cells out; "
               "set upstream.cell_retention = True if the TFF returns them."))
    tff_cap = up.tff_area_m2 * up.tff_capacity_Lm2
    if tff_cap > 0 and Fp > 0:
        life_d = tff_cap / (Fp * 24.0)
        notes.append(f"TFF set life {life_d:.1f} days at {Fp * 24:.0f} L/day "
                     f"({up.tff_area_m2:g} m² at {up.tff_capacity_Lm2:g} L/m²): "
                     f"about {max(0, math.ceil(p.run_days / life_d) - 1)} swap(s) in the run.")
    first_buf = daily[min(int(first["n_day"]), len(daily) - 1)]["buffer_L"] if first and daily else None
    peak_buf = max((r["buffer_L"] for r in daily), default=0.0)
    info = ("Downstream runs the Case 5 flowsheet step by step (Brij-9, Protein A, depth filtration, AEX, CEX, "
            "virus filtration, UF/DF, BDS fill). Step volumes are typical placeholders and the depth-filtration "
            "solution order is assumed: replace them with the process description.")
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
        "notes": [info] + (notes or ["Design is consistent with the titre curve."]),
    }


def process_design(p):
    tr = apply_seed_train_to_n1(p)            # seedtrain.apply_to_n1: the N-1 comes from the seed-train design
    d = steady_state_design(p) if p.upstream.mode == "mechanistic" else profile_design(p)
    d["n1"] = n1 = n1_summary(p)
    if tr is not None:
        s = tr["stages"][-1]
        d["notes"].append(f"N-1 is the batch N-1 stage of the seed train ({tr['basis_label']}): {s['V_L']:.1f} L, "
                          f"{s['duration_d']:.2f} d, {s['x0']:g} → {s['xt']:.2f} ×10⁶/mL; seeds N at "
                          f"{n1['n_seed_vcd']} ×10⁶/mL.")
    if n1.get("enabled") and n1["reason"] not in ("target reached", "fixed N-1 duration"):
        d["notes"].append(f"N-1: {n1['reason']}. N is seeded at {n1['n_seed_vcd']} ×10⁶/mL after {n1['days']} d "
                          f"(target {p.n1.target_seed_vcd:g}). Raise transfer_L, X_max or max_days.")
    return d


def portfolio_sizing(p):
    names = profile_names(p)
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
    notes = [f"Size solution supply and column cycling for the {buf_case['label'].lower()}: "
             f"{buf_case['peak_buffer_L_day']} L/day at {buf_case['peak_cycles_day']} cycles/day."]
    if sump_case:
        notes.append(f"Size the harvest sump for the {sump_case['label'].lower()}: "
                     f"{sump_case['first_cycle_sump_L']} L before its first cycle (N day "
                     f"{sump_case['first_cycle_n_day']}).")
    missing = [r["label"] for r in rows if r["first_cycle_sump_L"] is None]
    if missing:
        notes.append("Never collects one column load in this run: " + ", ".join(missing) + ".")
    return {"rows": rows, "buffer_sizing_case": buf_case["molecule"],
            "sump_sizing_case": sump_case["molecule"] if sump_case else None,
            "peak_buffer_L_day": top, "notes": notes}
