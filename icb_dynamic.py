"""
icb_dynamic.py — Hour-by-hour simulation of the IDB line
NavNevAI_mAbWasteWater · CarbonTesseract 2026

Timeline: N-1 phase first (permeate to waste), then the N run. "day" in the
output is campaign day; "n_day" is days since N inoculation (negative in N-1).

  N-1            media in = permeate out to waste; seed grows, sets N seed VCD
  N + TFF        profile mode: titre from the molecule's permeate curve,
                 feed = outflow, a share leaves the TFF as cell concentrate;
                 mechanistic mode: legacy ATF perfusion model
  Harvest sump   collects permeate; PAT titre x volume = grams in the sump
  Capture        batch: one cycle each time the sump holds one load;
                 pcc: 3-column PCC from a level-controlled surge
  VI             PI acid dose; out-of-spec or short-hold material diverted
  Polishing      CEX, AEX, virus filter, SPTFF; buffer per gram
  Buffer         tracked per unit operation every hour
  NevAI          optional, advisory: 24 h titre forecast -> grams -> buffer

Disturbances (N day), for example
  {"type": "qp_step", "day": 10, "factor": 0.85}    scales titre (profile) or qP
  {"type": "feed_glucose_step", "day": 8, "factor": 0.7}   legacy model only
  {"type": "vcd_sensor_bias", "day": 9, "bias": 0.10}      legacy model only

Quick check:  python3 icb_dynamic.py
"""
from __future__ import annotations
import json
import math
from collections import deque
import numpy as np
from icb_process import (ICBParams, WASTE_STREAMS, BUFFER_OPS, process_design, capture_cycle,
                         buffer_L_per_g, titre_at, profile_points)

SERIES_KEYS = ["h", "day", "n_day", "phase", "X_n1", "Xv", "Xv_meas", "viability", "S",
               "titre_harvest", "titre_meas", "sieve", "VVD", "Fb_Lh", "Fh_Lh",
               "sump_L", "sump_g", "Fload_Lh", "cycles", "col_state",
               "eluate_surge_L", "F_vi_Lh", "vi_pH", "acid_frac",
               "harvest_g_h", "ds_g_cum", "buffer_L_h", "ww_total_Lh",
               "titre_fc_24h", "switch_fc_h", "ww_fc_24h_L", "g_fc_24h", "buffer_fc_24h_L",
               "g_next24", "buffer_next24_L"]


def _clamp(x, lo, hi):
    return lo if x < lo else hi if x > hi else x


def _rnd(v):
    return round(v, 4) if isinstance(v, float) and math.isfinite(v) else v


class PI:
    """Positional PI controller with clamped output and anti-windup."""
    def __init__(self, kc, ti_h, lo, hi, bias=0.0):
        self.kc, self.ti, self.lo, self.hi, self.bias = kc, max(ti_h, 1e-6), lo, hi, bias
        self.integral = 0.0

    def step(self, error, dt=1.0):
        trial = self.integral + error * dt
        u = self.bias + self.kc * (error + trial / self.ti)
        if self.lo <= u <= self.hi:
            self.integral = trial
            return u
        return _clamp(self.bias + self.kc * (error + self.integral / self.ti), self.lo, self.hi)


class ICBSimulator:
    def __init__(self, params=None, forecaster=None, seed=0, forecast_every_h=6):
        self.p = params if params is not None else ICBParams()
        self.fc = forecaster if (forecaster is not None and getattr(forecaster, "ready", False)) else None
        self.rng = np.random.default_rng(seed)
        self.fc_every = max(1, int(forecast_every_h))

    def _disturb(self, n_day):
        qp, sf, bias = 1.0, 1.0, 0.0
        for d in self.p.disturbances or []:
            if n_day >= float(d.get("day", 1e9)):
                kind = d.get("type")
                if kind == "qp_step":
                    qp *= float(d.get("factor", 1.0))
                elif kind == "feed_glucose_step":
                    sf *= float(d.get("factor", 1.0))
                elif kind == "vcd_sensor_bias":
                    bias += float(d.get("bias", 0.0))
        return qp, sf, bias

    def run(self, keep_history=False, downsample_h=1):
        p = self.p
        up, n1, br, cp = p.upstream, p.n1, p.bioreactor, p.capture
        sg, pc, vi, po, sp, pat, ct = p.surge, p.pcc, p.vi, p.polish, p.sptff, p.pat, p.control
        rng = self.rng
        profile = up.mode != "mechanistic"
        batch = cp.mode != "pcc"
        V = br.V_L
        H1 = int(round(n1.days * 24)) if n1.enabled else 0
        H = H1 + int(round(p.run_days * 24))
        n_sub, dt = 4, 0.25
        design = process_design(p)
        n1s = design["n1"]
        cyc = capture_cycle(p)
        per_g = buffer_L_per_g(p)
        pts = profile_points(p)
        target = cyc["target_g"]

        # N-1 seed
        X1 = n1.X0
        # N reactor (cell balances used by the legacy model only)
        Xv, Xd, S, P = n1s["n_seed_vcd"], 0.0, br.S0_gL, 0.0
        sieve, dev = br.sieving0, 0.0
        sig = up.variation_sigma if profile else br.qp_ou_sigma
        tau = max(up.variation_tau_h if profile else br.qp_ou_tau_h, 1e-6)
        swap_left, atf_swaps = 0, 0
        bleed_pi = PI(ct.bleed_kc, ct.bleed_ti_h, 0.0, br.bleed_max_VVD * V / 24.0,
                      bias=0.0 if profile else (design.get("bleed_Lh") or 0.0))
        bleed_on = False
        # harvest sump
        Vs_sp = sg.V_max_L * sg.level_sp
        sump_max = cp.sump_max_L if batch else sg.V_max_L
        Vs, Ms = (0.0, 0.0) if batch else (Vs_sp, 0.0)
        perm_L = 0.0
        # capture
        col, load_left, rec_left, cyc_g = 0, 0.0, 0.0, 0.0   # col: 0 idle, 1 loading, 2 recovering
        Q_dbc = pc.DBC_gL * pc.CV_L
        Q_sat = pc.sat_factor * Q_dbc
        rec_timers, m_true, m_est, pending = [], 0.0, 0.0, False
        cycles, wait_h, cycle_hours, first_cycle = 0, 0, [], None
        n_cyc_h = 0
        # eluate surge + VI
        Ve = Me = 0.0
        acid_pi = PI(ct.acid_kc, ct.acid_ti_h, 0.0, vi.acid_max, bias=vi.acid_bias)
        acid = vi.acid_bias
        tube = deque([[-1e9, vi.tube_V_L, 0.0, True]])   # [entry_h, L, g, conforming]
        # polishing / drug substance
        vf_used, vf_changes = 0.0, 0
        ds_g = ds_L = produced_g = media_L = 0.0
        lots = {}
        loss = {k: 0.0 for k in ("tff_retentate", "cell_bleed", "ramp_harvest", "sump_overflow",
                                 "capture", "breakthrough", "vi_divert", "polishing", "sptff")}
        ww_sum = {k: 0.0 for k in WASTE_STREAMS}
        buf = {k: 0.0 for k in BUFFER_OPS}
        events, alarms, alarm_hours = [], {}, []
        history, fc_log = [], []
        ww = {k: 0.0 for k in WASTE_STREAMS}

        def event(h, code, msg):
            events.append({"h": h, "day": round(h / 24, 2), "type": "event", "code": code, "msg": msg})

        def alarm(h, code, msg):
            alarms[code] = alarms.get(code, 0) + 1
            alarm_hours.append(h)
            if alarms[code] <= 3:
                events.append({"h": h, "day": round(h / 24, 2), "type": "alarm", "code": code, "msg": msg})

        def elute(h, bound):
            nonlocal Ve, Me, cycles, n_cyc_h
            Ve += cyc["eluate_L"]
            Me += bound * cyc["yield"]
            loss["capture"] += bound * (1 - cyc["yield"])
            ww["capture_buffers"] += cyc["other_buffer_waste_L"]
            ww["capture_strip_cip"] += cyc["strip_cip_L"]
            buf["capture"] += cyc["buffer_L"]
            cycles += 1
            n_cyc_h += 1
            cycle_hours.append(h)
            if cycles == cyc["resin_cycles_max"]:
                event(h, "resin_life", "Protein A resin reached its cycle limit; plan replacement")

        for h in range(H):
            t_day = h / 24.0
            in_n = h >= H1
            n_h = h - H1
            n_day = n_h / 24.0
            ww = {k: 0.0 for k in WASTE_STREAMS}
            buf0 = dict(buf)
            Fin = Fb = Fh = harv_g = 0.0
            titre_h = titre_m = Xv_m = viab = None
            n_cyc_h = 0
            Fload_vol = 0.0

            if not in_n:
                # ── N-1: media in = permeate out to waste ──
                F1 = n1.VVD * n1.V_L / 24.0
                X1 = min(X1 * math.exp(n1.mu_h), n1.X_max)
                ww["n1_permeate"] += F1
                media_L += F1
            else:
                if n_h == 0:
                    event(h, "inoculate", f"N inoculated at {n1s['n_seed_vcd']:.2f} ×10⁶ cells/mL"
                          + (f" from N-1 at {X1:.1f} ×10⁶/mL" if n1.enabled else ""))
                qp_f, s_f, vcd_bias = self._disturb(n_day)
                dev += -dev / tau + sig * math.sqrt(2.0 / tau) * rng.normal()

                if profile:
                    # ── IDB: feed = outflow (cells + broth + product) -> TFF ──
                    Fin = up.VVD * V / 24.0
                    Fb = Fin * up.tff_retentate_frac          # cell concentrate
                    Fh = Fin - Fb                             # permeate
                    c = titre_at(p, n_day + 0.5 / 24.0, pts) * qp_f * math.exp(dev)
                    titre_h = c
                    harv_g = Fh * c
                    produced_g += Fin * c
                    media_L += Fin
                    ww["tff_retentate"] += Fb
                    loss["tff_retentate"] += Fb * c
                else:
                    # ── legacy ATF perfusion model ──
                    qP = br.qP_pg_cell_h * qp_f * math.exp(dev)
                    Xv_m = max(0.0, Xv * (1 + vcd_bias + pat.vcd_noise * rng.normal()))
                    if swap_left == 0 and sieve < br.sieve_replace_below:
                        swap_left, atf_swaps = br.atf_swap_h, atf_swaps + 1
                        ww["atf_flush"] += br.atf_swap_flush_L
                        event(h, "atf_swap", f"ATF filter replaced (sieving fell to {sieve:.2f})")
                        sieve = br.sieving0
                    swapping = swap_left > 0
                    if swapping:
                        swap_left -= 1
                    else:
                        Fin = _clamp(ct.CSPR_pL_cell_day * Xv_m * V * 1e-3 / 24.0,
                                     br.VVD_min * V / 24.0, br.VVD_max * V / 24.0)
                        if not bleed_on and Xv_m >= ct.bleed_enable_frac * br.X_setpoint:
                            bleed_on = True
                            event(h, "bleed_on", "VCD reached setpoint band; bleed control active")
                        Fb = min(bleed_pi.step(Xv_m - br.X_setpoint), Fin) if bleed_on else 0.0
                    Fh = Fin - Fb
                    media_L += Fin
                    bleed_g = 0.0
                    for _ in range(n_sub):
                        mu = br.mu_max_h * S / (br.Ks_gL + S) * max(0.0, 1.0 - Xv / br.X_max)
                        prod = qP * Xv * 1e-3
                        dXv = (mu - br.kd_h) * Xv - Fb / V * Xv
                        dXd = br.kd_h * Xv - (Fb / V + br.k_lysis_h) * Xd
                        dS = Fin / V * (br.S_feed_gL * s_f - S) - br.qS_pg_cell_h * Xv * 1e-3
                        dP = prod - (Fh * sieve + Fb) / V * P
                        harv_g += Fh * sieve * P * dt
                        bleed_g += Fb * P * dt
                        produced_g += prod * V * dt
                        Xv, Xd = max(Xv + dXv * dt, 0.0), max(Xd + dXd * dt, 0.0)
                        S, P = max(S + dS * dt, 0.0), max(P + dP * dt, 0.0)
                    if not swapping:
                        sieve = max(0.0, sieve - br.sieve_loss_per_day / 24.0)
                    ww["cell_bleed"] += Fb
                    loss["cell_bleed"] += bleed_g
                    titre_h = P * sieve
                    viab = Xv / (Xv + Xd) if (Xv + Xd) > 0 else 0.0
                    if bleed_on and viab < ct.viability_alarm:
                        alarm(h, "viability", f"Viability {100 * viab:.0f}% below alarm limit")
                    if S < ct.glucose_alarm_gL:
                        alarm(h, "glucose", f"Glucose {S:.2f} g/L: feed limiting")

                titre_m = max(0.0, titre_h * (1 + pat.titre_noise * rng.normal()))

                # ── harvest routing ──
                capturing = n_day >= cp.capture_start_day
                if capturing:
                    Vs += Fh
                    Ms += harv_g
                    perm_L += Fh
                else:
                    ww["ramp_harvest"] += Fh
                    loss["ramp_harvest"] += harv_g

                # ── capture ──
                if batch:
                    budget = 1.0
                    while budget > 1e-9:
                        if col == 0:
                            Cs = Ms / Vs if Vs > 1e-9 else 0.0
                            Cs_m = max(0.0, Cs * (1 + pat.titre_noise * rng.normal()))
                            if not (capturing and Cs_m > 1e-9 and Cs_m * Vs >= target):
                                break
                            load_left = min(Vs, target / Cs_m)
                            cyc_g, col = 0.0, 1
                            if first_cycle is None:
                                first_cycle = {"n_day": round(n_day + (1.0 - budget) / 24.0, 2),
                                               "sump_L": round(perm_L, 1)}
                                event(h, "first_cycle", f"First Protein A cycle after {perm_L:.0f} L "
                                      f"of permeate (N day {n_day:.1f})")
                        elif col == 1:
                            Fl = cyc["load_flow_Lh"]
                            step = min(budget, load_left / Fl)
                            Cs = Ms / Vs if Vs > 1e-9 else 0.0
                            Vt = min(Fl * step, Vs, load_left)
                            mt = Cs * Vt
                            Vs, Ms = max(Vs - Vt, 0.0), max(Ms - mt, 0.0)
                            load_left -= Vt
                            cyc_g += mt
                            Fload_vol += Vt
                            ww["capture_flowthrough"] += Vt
                            budget -= step
                            if load_left <= 1e-6 or Vs <= 1e-9:
                                if cyc_g > cyc["capacity_g"]:
                                    loss["breakthrough"] += cyc_g - cyc["capacity_g"]
                                    alarm(h, "capture_overload",
                                          "Column loaded above capacity; product lost in flow-through")
                                    cyc_g = cyc["capacity_g"]
                                elute(h, cyc_g)
                                col, rec_left = 2, cyc["recovery_h"]
                        else:
                            step = min(budget, rec_left)
                            rec_left -= step
                            budget -= step
                            if rec_left <= 1e-9:
                                col = 0
                    if col != 0 and capturing and Ms >= target:
                        wait_h += 1
                else:
                    # PCC: several column switches may fall inside one hour, so the
                    # hour is split at each switch (avoids a 1 h minimum per cycle)
                    budget, waited = 1.0, False
                    Fl = _clamp(Fh + (Vs - Vs_sp) / ct.surge_tau_h, 0.0, cyc["load_flow_Lh"])
                    noise = 1 + pat.titre_noise * rng.normal()

                    def tick(dt_h):
                        return [t - dt_h for t in rec_timers if t - dt_h > 1e-9]

                    while budget > 1e-9:
                        if capturing and m_est >= target - 1e-6:
                            if len(rec_timers) <= pc.n_columns - 3:
                                carry = pc.bt_frac * max(0.0, m_true - Q_dbc)   # sits on column 2
                                if first_cycle is None:
                                    first_cycle = {"n_day": round(n_day + (1.0 - budget) / 24.0, 2),
                                                   "sump_L": round(perm_L, 1)}
                                    event(h, "first_cycle", f"First Protein A cycle after {perm_L:.0f} L "
                                          f"of permeate (N day {n_day:.1f})")
                                elute(h, m_true - carry)
                                rec_timers.append(cyc["recovery_h"])
                                m_true = m_est = carry
                                pending = False
                                continue
                            if not pending:
                                pending = True
                                alarm(h, "pcc_wait", "Lead column full but no regenerated column ready; loading paused")
                            waited = True
                            step = min(budget, min(rec_timers)) if rec_timers else budget
                            rec_timers = tick(step)
                            budget -= step
                            continue
                        Cs = Ms / Vs if Vs > 1e-9 else 0.0
                        Cs_m = max(0.0, Cs * noise)
                        if not capturing or Fl <= 1e-9 or Vs <= 1e-9 or Cs_m <= 1e-9:
                            rec_timers = tick(budget)
                            break
                        step = min(budget, (target - m_est) / (Cs_m * Fl), Vs / Fl)
                        Vt = Fl * step
                        mt = Cs * Vt
                        Vs, Ms = max(Vs - Vt, 0.0), max(Ms - mt, 0.0)
                        m_true += mt
                        m_est += Cs_m * Vt
                        Fload_vol += Vt
                        ww["capture_flowthrough"] += Vt
                        if m_true > Q_sat:
                            loss["breakthrough"] += m_true - Q_sat
                            m_true = Q_sat
                            alarm(h, "pcc_breakthrough", "Lead column saturated; product lost in flow-through")
                        rec_timers = tick(step)
                        budget -= step
                        if Vs <= 1e-9:
                            rec_timers = tick(budget)
                            break
                    if waited:
                        wait_h += 1
                    col = 2 if rec_timers else (1 if Fload_vol > 0 else 0)

                # ── sump limits ──
                if Vs > sump_max:
                    over = Vs - sump_max
                    c = Ms / Vs
                    ww["sump_overflow"] += over
                    loss["sump_overflow"] += over * c
                    Vs -= over
                    Ms -= over * c
                    alarm(h, "sump_overflow", "Harvest sump overflowing; capture cannot keep up")
                if batch and capturing and Fh > 0 and Vs / Fh > cp.max_sump_hold_h:
                    alarm(h, "sump_hold", f"Mean product age in the sump about {Vs / Fh:.0f} h, "
                          f"above the {cp.max_sump_hold_h:.0f} h hold limit")

            # ── eluate surge -> continuous low-pH VI (plug flow) ──
            if Ve > vi.eluate_surge_max_L:
                alarm(h, "eluate_surge", "Eluate surge above capacity")
            F_vi_max = vi.tube_V_L / vi.t_hold_min_h / (1 + acid)
            F_vi = _clamp((Ve - vi.heel_L) / ct.eluate_tau_h, 0.0, F_vi_max)
            V_draw = min(F_vi, Ve)
            Ce = Me / Ve if Ve > 1e-9 else 0.0
            m_draw = Ce * V_draw
            Ve -= V_draw
            Me -= m_draw
            pH_out, v_out, m_out = None, 0.0, 0.0
            if V_draw > 1e-9:
                buf["vi"] += V_draw * acid
                pH_in = vi.eluate_pH_mean + vi.eluate_pH_sd * rng.normal()
                pH_out = pH_in - acid * vi.ph_drop_per_acid_frac
                pH_meas = pH_out + pat.ph_noise * rng.normal()
                ok = pH_out <= vi.pH_max
                if not ok:
                    alarm(h, "vi_ph", f"VI pH {pH_out:.2f} above limit; material will be diverted")
                vol_in = V_draw * (1 + acid)
                tube.append([h, vol_in, m_draw, ok])
                acid = acid_pi.step(pH_meas - vi.pH_target)
                need = vol_in
                while need > 1e-9 and tube:
                    pkt = tube[0]
                    take = min(need, pkt[1])
                    m_take = pkt[2] * take / pkt[1]
                    pkt[1] -= take
                    pkt[2] -= m_take
                    need -= take
                    if m_take <= 1e-12:
                        ww["vi_flush"] += take
                    elif (not pkt[3]) or (h - pkt[0]) < vi.t_hold_min_h - 1e-9:
                        ww["vi_divert"] += take
                        loss["vi_divert"] += m_take
                        if pkt[3]:
                            alarm(h, "vi_hold", "VI residence below minimum; material diverted")
                    else:
                        v_out += take
                        m_out += m_take
                    if pkt[1] <= 1e-9:
                        tube.popleft()

            # ── polishing -> VF -> SPTFF/inline DF -> formulation ──
            m_ds_h = 0.0
            if m_out > 0:
                buf["vi"] += v_out * vi.base_frac
                v_neut = v_out * (1 + vi.base_frac)
                m1 = m_out * vi.vi_yield
                ww["cex_flowthrough"] += v_neut
                ww["cex_buffers"] += m1 * po.cex_buffer_L_per_g
                buf["cex"] += m1 * po.cex_buffer_L_per_g
                m2 = m1 * po.cex_yield
                v2 = m2 / po.cex_pool_gL
                ww["aex_buffers"] += m2 * po.aex_buffer_L_per_g
                buf["aex"] += m2 * po.aex_buffer_L_per_g
                m3 = m2 * po.aex_yield
                v3 = v2 * (1 + po.aex_flush_frac)
                vf_used += v3
                if vf_used >= po.vf_capacity_L:
                    vf_used -= po.vf_capacity_L
                    vf_changes += 1
                    ww["vf_flush"] += po.vf_flush_L
                    buf["vf"] += po.vf_flush_L
                    event(h, "vf_change", "Virus filter changed (capacity reached)")
                m4 = m3 * po.vf_yield
                m_ds_h = m4 * sp.sptff_yield
                v_ret = m_ds_h / sp.target_gL
                ww["sptff_permeate"] += max(0.0, v3 + sp.DV * v_ret - v_ret)
                buf["sptff"] += v_ret * (sp.DV + sp.excipient_frac)
                v_ds = v_ret * (1 + sp.excipient_frac)
                loss["polishing"] += m_out - m4
                loss["sptff"] += m4 - m_ds_h
                ds_g += m_ds_h
                ds_L += v_ds
                lot_id = int(max(n_h, 0) // p.lot_window_h)
                lot = lots.setdefault(lot_id, {"lot": lot_id + 1, "g": 0.0, "L": 0.0})
                lot["g"] += m_ds_h
                lot["L"] += v_ds

            # ── record ──
            ww_h = sum(ww.values())
            for k, v in ww.items():
                ww_sum[k] += v
            mech = (not profile) and in_n
            rec = {
                "h": h, "day": round(t_day, 3), "n_day": round(n_day, 3), "phase": int(in_n),
                "X_n1": None if in_n else X1,
                "Xv": Xv if mech else None, "Xv_meas": Xv_m, "viability": viab,
                "S": S if mech else None, "sieve": sieve if mech else None,
                "titre_harvest": titre_h, "titre_meas": titre_m,
                "Fin_Lh": Fin, "Fb_Lh": Fb, "Fh_Lh": Fh, "VVD": Fin * 24.0 / V,
                "bleed_frac": Fb / Fin if Fin > 0 else 0.0, "bleed_on": int(bleed_on),
                "sump_L": Vs, "sump_g": Ms, "Fload_Lh": Fload_vol, "cycles": cycles,
                "capture_switch": n_cyc_h, "col_state": col,
                "eluate_surge_L": Ve, "F_vi_Lh": F_vi, "vi_pH": pH_out, "acid_frac": acid,
                "harvest_g_h": harv_g, "ds_g_h": m_ds_h, "ds_g_cum": ds_g,
                "buffer_L_h": sum(buf.values()) - sum(buf0.values()),
                "ww_total_Lh": ww_h,
                "ww_cells_Lh": ww["cell_bleed"] + ww["tff_retentate"] + ww["atf_flush"],
                "titre_fc_24h": None, "switch_fc_h": None, "ww_fc_24h_L": None,
                "g_fc_24h": None, "buffer_fc_24h_L": None,
            }
            for k in BUFFER_OPS:
                rec["buf_" + k] = buf[k] - buf0[k]
            for k, v in ww.items():
                rec["ww_" + k] = v
            history.append(rec)

            # ── NevAI forecasts (advisory) ──
            if (self.fc is not None and in_n and h % self.fc_every == 0
                    and len(history) >= self.fc.seq_len):
                f = self.fc.forecast_all(history)
                tf, wf = f.get("titre"), f.get("ww")
                if tf is not None:
                    rec["titre_fc_24h"] = float(tf[-1])
                    g_fc = float(np.sum(tf)) * Fh
                    rec["g_fc_24h"] = g_fc
                    rec["buffer_fc_24h_L"] = g_fc * per_g["total"]
                    have = Ms if batch else m_est
                    remaining, cum = max(0.0, target - have), 0.0
                    for k, c in enumerate(tf):
                        cum += Fh * float(c)
                        if cum >= remaining:
                            rec["switch_fc_h"] = k + 1
                            fc_log.append((h, h + k + 1))
                            break
                    if titre_m and titre_m > 0.05 and tf[-1] < (1 - ct.titre_drop_alarm) * titre_m:
                        alarm(h, "fc_titre_drop",
                              f"NevAI forecasts permeate titre falling to {float(tf[-1]):.2f} g/L within {len(tf)} h")
                if wf is not None:
                    rec["ww_fc_24h_L"] = float(np.sum(wf))

        # actual next-24 h buffer and grams, to score the forecasts
        b = np.concatenate([[0.0], np.cumsum([r["buffer_L_h"] for r in history])])
        g = np.concatenate([[0.0], np.cumsum([r["harvest_g_h"] for r in history])])
        n = len(history)
        for i, r in enumerate(history):
            ok = r["phase"] and i + 24 <= n
            r["buffer_next24_L"] = float(b[i + 24] - b[i]) if ok else None
            r["g_next24"] = float(g[i + 24] - g[i]) if ok else None

        st = dict(ds_g=ds_g, ds_L=ds_L, produced_g=produced_g, media_L=media_L, loss=loss,
                  ww_sum=ww_sum, buf=buf, events=events, alarms=alarms, alarm_hours=alarm_hours,
                  lots=lots, cycle_hours=cycle_hours, first_cycle=first_cycle, cycles=cycles,
                  wait_h=wait_h, atf_swaps=atf_swaps, vf_changes=vf_changes, fc_log=fc_log,
                  cyc=cyc, per_g=per_g, n1s=n1s)
        return self._summarise(p, design, H, history, keep_history, downsample_h, st)

    def _summarise(self, p, design, H, history, keep_history, downsample_h, st):
        V = p.bioreactor.V_L
        ds_g, produced_g, media_L = st["ds_g"], st["produced_g"], st["media_L"]
        ww_sum, buf = st["ww_sum"], st["buf"]
        n_rows = [r for r in history if r["phase"]]
        mech = p.upstream.mode == "mechanistic"

        daily = []
        for d in range(int(math.ceil(H / 24))):
            rows = history[d * 24:(d + 1) * 24]
            if not rows:
                break
            by = {k: sum(r["ww_" + k] for r in rows) for k in WASTE_STREAMS}
            row = {"day": d + 1, "phase": "N" if any(r["phase"] for r in rows) else "N-1",
                   "total_L": round(sum(by.values()), 1),
                   "harvest_g": round(sum(r["harvest_g_h"] for r in rows), 1),
                   "cycles": sum(r["capture_switch"] for r in rows),
                   "buffer_L": round(sum(r["buffer_L_h"] for r in rows), 1),
                   "ds_g": round(sum(r["ds_g_h"] for r in rows), 1),
                   "alarms": sum(1 for a in st["alarm_hours"] if d * 24 <= a < (d + 1) * 24)}
            for k, v in by.items():
                row["L_" + k] = round(v, 1)
            for k in BUFFER_OPS:
                row["buf_" + k] = round(sum(r["buf_" + k] for r in rows), 1)
            daily.append(row)
        nd = [d for d in daily if d["phase"] == "N"] or daily
        cap_days = [d["total_L"] for d in nd] or [0.0]

        ww_total = sum(ww_sum.values())
        by_route, by_driver = {}, {}
        for k, v in ww_sum.items():
            _, route, driver, _ = WASTE_STREAMS[k]
            by_route[route] = by_route.get(route, 0.0) + v
            by_driver[driver] = by_driver.get(driver, 0.0) + v

        def pct(x):
            return round(100 * x / ww_total, 1) if ww_total > 0 else None

        buf_total = sum(buf.values())
        bdays = [d["buffer_L"] for d in nd if d["buffer_L"] > 0]
        ch = np.array(st["cycle_hours"], dtype=float)
        errs = []
        for h_made, h_pred in st["fc_log"]:
            nxt = ch[ch > h_made]
            if nxt.size:
                errs.append(abs(h_pred - float(nxt[0])))
        mape = [abs(r["buffer_fc_24h_L"] - r["buffer_next24_L"]) / r["buffer_next24_L"]
                for r in history if r.get("buffer_fc_24h_L") is not None
                and r.get("buffer_next24_L") is not None and r["buffer_next24_L"] > 1.0]
        fcyc = st["first_cycle"] or {}
        last_t = next((r["titre_harvest"] for r in reversed(history) if r["titre_harvest"] is not None), None)
        if mech:
            steady = [r for r in n_rows if r["bleed_on"]] or n_rows
            s_vcd = round(float(np.mean([r["Xv"] for r in steady])), 1) if steady else None
            s_via = round(100 * float(np.mean([r["viability"] for r in steady])), 1) if steady else None
            s_tit = round(float(np.mean([r["titre_harvest"] for r in steady])), 3) if steady else None
        else:
            s_vcd = s_via = s_tit = None

        kpi = {
            "drug_substance_kg": round(ds_g / 1000, 3),
            "ds_volume_L": round(st["ds_L"], 1),
            "product_made_kg": round(produced_g / 1000, 3),
            "overall_yield_pct": round(100 * ds_g / produced_g, 1) if produced_g > 0 else 0.0,
            "volumetric_productivity_g_L_day": round(ds_g / V / p.run_days, 3),
            "n1_seed_vcd": st["n1s"]["n_seed_vcd"],
            "n1_waste_m3": round(ww_sum["n1_permeate"] / 1000, 3),
            "media_used_m3": round(media_L / 1000, 2),
            "media_m3_per_kg_ds": round(media_L / ds_g, 2) if ds_g > 0 else None,
            "final_titre_gL": round(last_t, 3) if last_t is not None else None,
            "steady_vcd_mean": s_vcd, "steady_viability_pct": s_via, "steady_titre_gL": s_tit,
            "capture_cycles": st["cycles"],
            "first_cycle_n_day": fcyc.get("n_day"),
            "first_cycle_sump_L": fcyc.get("sump_L"),
            "capture_mean_interval_h": round(float(np.mean(np.diff(ch))), 2) if ch.size > 1 else None,
            "capture_recovery_h": round(st["cyc"]["recovery_h"], 2),
            "capture_wait_h": st["wait_h"],
            "peak_cycles_day": max((d["cycles"] for d in nd), default=0),
            "peak_g_day": max((d["harvest_g"] for d in nd), default=0.0),
            "buffer_m3": round(buf_total / 1000, 3),
            "buffer_L_per_g_ds": round(buf_total / ds_g, 3) if ds_g > 0 else None,
            "buffer_daily_peak_L": max(bdays) if bdays else 0.0,
            "buffer_peak_to_avg": round(max(bdays) / float(np.mean(bdays)), 2) if bdays else None,
            "atf_swaps": st["atf_swaps"],
            "vf_changes": st["vf_changes"],
            "wastewater_m3": round(ww_total / 1000, 2),
            "wastewater_m3_per_kg_ds": round(ww_total / ds_g, 2) if ds_g > 0 else None,
            "ww_daily_avg_m3": round(float(np.mean(cap_days)) / 1000, 3),
            "ww_daily_peak_m3": round(max(cap_days) / 1000, 3),
            "ww_peak_to_avg": round(max(cap_days) / max(float(np.mean(cap_days)), 1e-9), 2),
            "ww_hourly_peak_L": round(max(r["ww_total_Lh"] for r in history), 1) if history else 0.0,
            "ww_volume_driven_pct": pct(by_driver.get("volume", 0.0)),
            "ww_mass_driven_pct": pct(by_driver.get("mass", 0.0)),
            "ww_event_pct": pct(by_driver.get("event", 0.0)),
            "losses_g": {k: round(v, 1) for k, v in st["loss"].items()},
            "alarms": st["alarms"],
            "nevai_switch_forecast_mae_h": round(float(np.mean(errs)), 2) if errs else None,
            "nevai_buffer_forecast_mape_pct": round(100 * float(np.mean(mape)), 1) if mape else None,
        }
        step = max(1, int(downsample_h))
        out = {
            "ok": True, "params": p.to_dict(), "design": design, "kpi": kpi,
            "wastewater": {
                "by_stream_m3": {k: round(v / 1000, 3) for k, v in ww_sum.items()},
                "by_route_m3": {k: round(v / 1000, 3) for k, v in by_route.items()},
                "by_driver_m3": {k: round(v / 1000, 3) for k, v in by_driver.items()},
                "streams": {k: {"area": a, "route": r, "driver": dr, "desc": ds}
                            for k, (a, r, dr, ds) in WASTE_STREAMS.items()},
                "daily": daily,
            },
            "buffer": {"by_op_m3": {k: round(v / 1000, 4) for k, v in buf.items()},
                       "planning_L_per_g": {k: round(v, 4) for k, v in st["per_g"].items()}},
            "lots": [{"lot": L["lot"], "g": round(L["g"], 1), "L": round(L["L"], 2)}
                     for _, L in sorted(st["lots"].items())],
            "series": [{k: _rnd(r.get(k)) for k in SERIES_KEYS} for r in history[::step]],
            "events": st["events"][:300],
        }
        if keep_history:
            out["_history"] = history
        return out


if __name__ == "__main__":
    res = ICBSimulator(seed=1).run()
    d = dict(res["design"])
    d.pop("daily", None)
    print("Design:")
    print(json.dumps(d, indent=2))
    print("KPIs:")
    print(json.dumps(res["kpi"], indent=2))
