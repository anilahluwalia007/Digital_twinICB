"""
icb_rilveg.py — Rilveg (2936) recipe set for the IDB mAb digital twin
NavNevAI_mAbWasteWater · CarbonTesseract 2026

Values are read from the Rilveg 2936 basis-of-design map. Anything the map does
not state is left at the platform default and marked below.

Use:
    from icb_rilveg import apply_rilveg
    apply_rilveg(p)                 # p is an ICBParams
or post RILVEG_OVERRIDES as "params" to /api/icb/simulate (call
register_rilveg_solutions() once at start-up so the solution keys resolve).

Case 4 in icb_cases.py applies this set automatically (RECIPES = {"4": "rilveg"}).

Not in the map, so not set here: Protein A, AEX and CEX column volumes, depth-filter
and virus-filter areas, UF/DF membrane area. Size them from the Sizing tab
(apply the suggested values and re-run) before reading Case 4 cost of goods.
"""
from __future__ import annotations
from icb_process import (SOLUTIONS, DILUTIONS, SOLUTION_ALIASES,
                         _acetate, _nacl, _mix, _sol)

_HIS, _ARG_HCL = 155.15, 210.66          # histidine, arginine-HCl (g/mol)

RILVEG_SOLUTIONS = {
    "lfs_elute": _sol("25 mM sodium acetate pH 3.6", "LambdaFabSelect elution",
                      _acetate(25, 3.6)),
    "lfs_strip2": _sol("120 mM phosphoric acid, 167 mM acetic acid, 2.2 % benzyl alcohol",
                       "LambdaFabSelect second strip",
                       {"Phosphoric acid (as 100 %)": 0.120 * 98.00,
                        "Acetic acid (glacial)": 0.167 * 60.05,
                        "Benzyl alcohol": 22.0}),
    "rv_cex_eq": _sol("50 mM sodium acetate, 40 mM NaCl, pH 5.5",
                      "Rilveg CEX equilibration and wash", _mix(_acetate(50, 5.5), _nacl(40))),
    "rv_cex_elute": _sol("50 mM sodium acetate, 130 mM NaCl, pH 5.5",
                         "Rilveg CEX elution and virus filtration",
                         _mix(_acetate(50, 5.5), _nacl(130))),
    "rv_ufdf": _sol("20 mM histidine/histidine-HCl, 240 mM arginine-HCl, pH 6.0",
                    "Rilveg UF/DF equilibration, diafiltration and recovery",
                    {"L-histidine": 0.020 * _HIS, "L-arginine hydrochloride": 0.240 * _ARG_HCL}),
    "rv_pfr_spike": _sol("10 g/L sodium chloride", "Plug-flow reactor spike",
                         {"Sodium chloride": 10.0}),
}
RILVEG_DILUTIONS = {"naoh_25mm": (("naoh_1n", 0.025), ("wfi", 0.975))}
RILVEG_ALIASES = {"naoh_25mm": "25 mM NaOH (1 N diluted in-line with WFI)"}


def register_rilveg_solutions():
    """Add the Rilveg solutions to the shared tables (safe to call twice)."""
    SOLUTIONS.update(RILVEG_SOLUTIONS)
    DILUTIONS.update(RILVEG_DILUTIONS)
    SOLUTION_ALIASES.update(RILVEG_ALIASES)


# Every value below is from the Rilveg 2936 map unless the comment says otherwise.
RILVEG_OVERRIDES = {
    "run_days": 14.0,
    "bioreactor": {"V_L": 6000.0},
    "upstream": {"VVD": 1.2, "cell_retention": True, "bleed_frac": 0.02,   # bleed not stated
                 "sieve0": 0.98, "sieve_loss_per_day": 0.006,              # passage not stated
                 "tff_area_m2": 7.2, "tff_capacity_Lm2": 2000.0,
                 "tff_swap_flush_L": 240.0},
    "capture": {"trains": 2, "sump_max_L": 2000.0,
                "load_flow_CVh": 5.0,            # about 100 cm/h at a 20 cm bed
                "postload_CV": 3.0, "wash1_CV": 0.0, "wash2_CV": 0.0,
                "elute_CV": 3.5, "eluate_CV": 3.5, "elute_pH": 3.6,
                "strip_CV": 3.0, "cip_CV": 3.0, "cip_contact_h": 0.5,
                "reequil_CV": 11.0, "capture_yield": 0.80, "load_util": 0.90,
                "resin_cycles_max": 50,
                "sol_postload": "tris_nacl_74", "sol_elute": "lfs_elute",
                "sol_strip": "lfs_strip2", "sol_cip": "naoh_25mm",
                "sol_equil": "tris_74"},
    "vi": {"det_target_pct": 0.015, "det_conc_pct": 10.0,
           "det_min_contact_min": 4.0, "det_line_hold_min": 4.0,
           "eluate_surge_max_L": 1000.0},
    "depth": {"capacity_Lm2": 700.0},
    "aex": {"columns": 2, "load_gL": 259.0, "chase_CV": 3.0, "reequil_CV": 3.0,
            "strip_CV": 8.0, "cip_CV": 3.0, "cip_contact_h": 0.25},
    "cex": {"columns": 2, "DBC_gL": 33.0, "load_util": 0.90,   # load_util not stated
            "wash_CV": 3.0, "reequil_CV": 3.0, "elute_CV": 4.0, "eluate_CV": 1.5,
            "strip_CV": 3.0, "cip_CV": 3.0, "cip_contact_h": 0.25,
            "sol_eq": "rv_cex_eq", "sol_elute": "rv_cex_elute"},
    "polish": {"aex_yield": 0.97, "cex_yield": 0.92, "prefilter_per_day": 2.0},
    "vrf": {"sol_equil": "rv_cex_elute", "integrity_pre_h": 2.0,
            "integrity_post_h": 2.0, "max_hold_h": 24.0},
    "sptff": {"uf1_gL": 70.0, "DV": 6.0, "target_gL": 200.0, "excipient_frac": 0.0,
              "sol_df": "rv_ufdf", "sol_recover": "rv_ufdf"},
    "sizing": {"vf_Lm2": 700.0, "ufdf_load_gm2": 1400.0, "tff_flux_LMH": 70.0},
}


def apply_rilveg(p, verbose=True):
    """Set the Rilveg recipe on an ICBParams. Warns about values the parameter
    set does not have, so a typo is never dropped silently."""
    register_rilveg_solutions()
    missing = []
    for sec, kv in RILVEG_OVERRIDES.items():
        if not isinstance(kv, dict):
            if hasattr(p, sec):
                setattr(p, sec, kv)
            else:
                missing.append(sec)
            continue
        obj = getattr(p, sec, None)
        if obj is None:
            missing.append(sec)
            continue
        for k, v in kv.items():
            if hasattr(obj, k):
                setattr(obj, k, v)
            else:
                missing.append(f"{sec}.{k}")
    if missing and verbose:
        print("icb_rilveg: not in the parameter set, skipped: " + ", ".join(missing))
    return p


register_rilveg_solutions()      # importing the module makes the solution keys available


if __name__ == "__main__":
    from icb_process import ICBParams
    apply_rilveg(ICBParams())        # prints any key the parameter set does not have
    print(f"Rilveg recipe set: {len(RILVEG_SOLUTIONS)} solutions, "
          f"{sum(len(v) if isinstance(v, dict) else 1 for v in RILVEG_OVERRIDES.values())} values")
