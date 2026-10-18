"""
icb_dynamic.py — Hour-by-hour simulation of the IDB line
NavNevAI_mAbWasteWater · CarbonTesseract 2026

Timeline: N-1 first (permeate to waste), then the N run, then a downstream run-out
(dsp.tail_h) with the bioreactor stopped so the last product is processed.
"day" = campaign day; "n_day" = days since N inoculation (negative in N-1).

  N-1 / N + TFF   unchanged (profile = titre curve, mechanistic = legacy ATF model)
  Harvest sump    conditioned media surge tank
  Downstream      icb_dsp.DSP: every step of the Case 5 flowsheet, in sequence
  Low-pH VI       legacy option (vi.mode = "low_ph"), between Protein A and depth filtration
  PCC             legacy option (capture.mode = "pcc"), Protein A steps lumped per cycle
  NevAI           optional, advisory: 24 h titre forecast -> grams -> solutions

Quick check:  python3 icb_dynamic.py
"""
from __future__ import annotations
import json
import math
from collections import deque
import numpy as np
from icb_process import (ICBParams, WASTE_STREAMS, BUFFER_OPS, SOLUTIONS, SOLUTION_ALIASES,
                         process_design, capture_cycle, buffer_L_per_g, titre_at, profile_points, det_frac)
from icb_dsp import DSP, Plant, Vessel, FLOW_KEYS, VESSEL_KEYS, LANE_KEYS, WFI_KEYS

_BASE_KEYS = ["h", "day", "n_day", "phase", "X_n1", "Xv", "Xv_meas", "viability", "S",
              "titre_harvest", "titre_meas", "sieve", "VVD", "Fb_Lh", "Fh_Lh",
              "sump_L", "sump_g", "Fload_Lh", "cycles", "col_state",
              "eluate_surge_L", "F_vi_Lh", "vi_pH", "acid_frac",
              "harvest_g_h", "ds_g_cum", "buffer_L_h", "ww_total_Lh",
              "titre_fc_24h", "switch_fc_h", "ww_fc_24h_L", "g_fc_24h", "buffer_fc_24h_L",
              "g_next24", "buffer_next24_L"]
SERIES_KEYS = list(dict.fromkeys(
    _BASE_KEYS + ["Fin_Lh", "capture_switch", "ds_g_h", "usp_on", "det_Lh", "ds_L_cum"]
    + list(FLOW_KEYS) + list(VESSEL_KEYS) + list(LANE_KEYS) + list(WFI_KEYS)
    + ["ww_" + k for k in WASTE_STREAMS] + ["buf_" + k for k in BUFFER_OPS]))


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
        sg, pc, vi, pat, ct = p.surge, p.pcc, p.vi, p.pat, p.control
        rng = self.rng
        profile = up.mode != "mechanistic"
        batch = cp.mode != "pcc"
        lowph = vi.mode == "low_ph"
        V = br.V_L
        H1 = int(round(n1.days * 24)) if n1.enabled else 0
        HN = int(round(p.run_days * 24))
        HT = int(round(max(p.dsp.tail_h, 0.0)))
        H = H1 + HN + HT
        n_sub, dt = 4, 0.25
        design = process_design(p)
        n1s = design["n1"]
        cyc = capture_cycle(p)
        per_g = buffer_L_per_g(p)
        dose = det_frac(p)
        pts = profile_points(p)
        target = cyc["target_g"]

        X1 = n1.X0
        Xv, Xd, S, P = n1s["n_seed_vcd"], 0.0, br.S0_gL, 0.0
        sieve, dev = br.sieving0, 0.0
        sig = up.variation_sigma if profile else br.qp_ou_sigma
        tau = max(up.variation_tau_h if profile else br.qp_ou_tau_h, 1e-6)
        swap_left, atf_swaps = 0, 0
        bleed_pi = PI(ct.bleed_kc, ct.bleed_ti_h, 0.0, br.bleed_max_VVD * V / 24.0,
                      bias=0.0 if profile else (design.get("bleed_Lh") or 0.0))
        bleed_on = False

        events, alarms, alarm_hours = [], {}, []

        def event(h, code, msg):
            events.append({"h": h, "day": round(h / 24, 2), "type": "event", "code": code, "msg": msg})

        def alarm(h, code, msg):
            alarms[code] = alarms.get(code, 0) + 1
            alarm_hours.append(h)
            if alarms[code] <= 3:
                events.append({"h": h, "day": round(h / 24, 2), "type": "alarm", "code": code, "msg": msg})

        loss = {k: 0.0 for k in ("tff_retentate", "cell_bleed", "ramp_harvest", "sump_overflow", "capture",
                                 "breakthrough", "vi", "vi_divert", "depth", "polishing", "sptff")}
        pl = Plant(p, rng, event, alarm, loss)
        led = pl.ledger
        Vs_sp = sg.V_max_L * sg.level_sp
        sump_max = cp.sump_max_L if batch else sg.V_max_L
        sump = Vessel("sump", "Conditioned media surge tank")
        if not batch:
            sump.L = Vs_sp
        dsp = DSP(pl, sump)
        perm_L = 0.0
        # PCC (legacy option)
        Q_dbc = pc.DBC_gL * pc.CV_L
        Q_sat = pc.sat_factor * Q_dbc
        rec_timers, m_true, m_est, pending = [], 0.0, 0.0, False
        pcc_charges = [("pa_wash1", pc.wash_CV / 2 * pc.CV_L, "capture_buffers"),
                       ("pa_wash2", pc.wash_CV / 2 * pc.CV_L, "capture_buffers"),
                       ("pa_strip", pc.strip_CV * pc.CV_L, "capture_strip_cip"),
                       ("naoh_05n", pc.cip_CV * pc.CV_L, "capture_strip_cip"),
                       ("tris_74", pc.reequil_CV * pc.CV_L, "capture_buffers")]
        wait_h, col = 0, 0
        # low-pH VI (legacy option)
        acid_pi = PI(ct.acid_kc, ct.acid_ti_h, 0.0, vi.acid_max, bias=vi.acid_bias)
        acid = vi.acid_bias
        tube = deque([[-1e9, vi.tube_V_L, 0.0, True]])
        vi_t0, vi_L = None, 0.0
        produced_g = media_L = 0.0
        ww_sum = {k: 0.0 for k in WASTE_STREAMS}
        buf = led.buf
        history, fc_log = [], []

        def pcc_elute(h, bound):
            dsp.pav.add(cyc["eluate_L"], bound * cyc["yield"], h + 1.0)
            loss["capture"] += bound * (1 - cyc["yield"])
            for sol, L, route in pcc_charges:
                led.use(sol, L, "capture")
                led.waste(route, L)
            el = pc.elute_CV * pc.CV_L
            led.use("pa_elute", el, "capture")
            led.waste("capture_buffers", el - cyc["eluate_L"])
            pl.cycles += 1
            pl.n_cyc_h += 1
            pl.cycle_hours.append(h)
            led.log("pa1", "Protein A (PCC)", "PCC elution", "pa_elute", "elute", h, h + 0.1, el)
            if pl.cycles == cyc["resin_cycles_max"]:
                event(h, "resin_life", "Protein A resin reached its cycle limit; plan replacement")

        for h in range(H):
            t_day = h / 24.0
            in_n = h >= H1
            n_h = h - H1
            n_day = n_h / 24.0
            usp_on = in_n and n_h < HN
            tail = in_n and not usp_on
            ww = {k: 0.0 for k in WASTE_STREAMS}
            led.ww = ww
            led.wfi_h = dict.fromkeys(BUFFER_OPS, 0.0)
            pl.flow = dict.fromkeys(FLOW_KEYS, 0.0)
            pl.h, pl.n_day, pl.n_cyc_h, pl.Fload, pl.det_L = h, n_day, 0, 0.0, 0.0
            buf0 = dict(buf)
            Fin = Fb = Fh = harv_g = 0.0
            titre_h = titre_m = Xv_m = viab = None
            Fload_vol = 0.0
            pH_out, F_vi, vi_on = None, 0.0, False
            capturing = False

            if not in_n:
                F1 = n1.VVD * n1.V_L / 24.0
                X1 = min(X1 * math.exp(n1.mu_h), n1.X_max)
                ww["n1_permeate"] += F1
                media_L += F1
            else:
                if n_h == 0:
                    event(h, "inoculate", f"N inoculated at {n1s['n_seed_vcd']:.2f} ×10⁶ cells/mL"
                          + (f" from N-1 at {X1:.1f} ×10⁶/mL" if n1.enabled else ""))
                if n_h == HN and HT > 0:
                    event(h, "run_out", "N run ended; downstream runs out the remaining product")
                if usp_on:
                    qp_f, s_f, vcd_bias = self._disturb(n_day)
                    dev += -dev / tau + sig * math.sqrt(2.0 / tau) * rng.normal()
                    if profile:
                        Fin = up.VVD * V / 24.0
                        Fb = Fin * up.tff_retentate_frac
                        Fh = Fin - Fb
                        c = titre_at(p, n_day + 0.5 / 24.0, pts) * qp_f * math.exp(dev)
                        titre_h = c
                        harv_g = Fh * c
                        produced_g += Fin * c
                        media_L += Fin
                        ww["tff_retentate"] += Fb
                        loss["tff_retentate"] += Fb * c
                    else:
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

                capturing = n_day >= cp.capture_start_day
                pl.capturing = capturing
                if usp_on:
                    if capturing:
                        sump.add(Fh, harv_g, h)
                        perm_L += Fh
                    else:
                        ww["ramp_harvest"] += Fh
                        loss["ramp_harvest"] += harv_g
                pl.perm_L = perm_L

                # ── PCC capture (legacy option) ──
                if not batch:
                    Vs, Ms = sump.L, sump.g
                    budget, waited = 1.0, False
                    Fl = (cyc["load_flow_Lh"] if tail else      # run-out: empty the surge tank
                          _clamp(Fh + (Vs - Vs_sp) / ct.surge_tau_h, 0.0, cyc["load_flow_Lh"]))
                    noise = 1 + pat.titre_noise * rng.normal()

                    def tick(dt_h):
                        return [x - dt_h for x in rec_timers if x - dt_h > 1e-9]

                    while budget > 1e-9:
                        last = tail and Vs <= 1e-9 and m_true > 1e-6     # run-out: elute the partial load
                        if capturing and (m_est >= target - 1e-6 or last):
                            if len(rec_timers) <= pc.n_columns - 3:
                                carry = 0.0 if last else pc.bt_frac * max(0.0, m_true - Q_dbc)
                                if pl.first_cycle is None:
                                    pl.first_cycle = {"n_day": round(n_day + (1.0 - budget) / 24.0, 2),
                                                      "sump_L": round(perm_L, 1)}
                                    event(h, "first_cycle", f"First Protein A cycle after {perm_L:.0f} L "
                                          f"of permeate (N day {n_day:.1f})")
                                pcc_elute(h, m_true - carry)
                                rec_timers.append(cyc["recovery_h"])
                                m_true = m_est = carry
                                pending = False
                                continue
                            if not pending and not last:
                                pending = True
                                alarm(h, "pcc_wait", "Lead column full but no regenerated column ready; loading paused")
                            waited = waited or not last
                            stp = min(budget, min(rec_timers)) if rec_timers else budget
                            rec_timers = tick(stp)
                            budget -= stp
                            continue
                        Cs = Ms / Vs if Vs > 1e-9 else 0.0
                        Cs_m = max(0.0, Cs * noise)
                        if not capturing or Fl <= 1e-9 or Vs <= 1e-9 or Cs_m <= 1e-9:
                            rec_timers = tick(budget)
                            break
                        stp = min(budget, (target - m_est) / (Cs_m * Fl), Vs / Fl)
                        Vt = Fl * stp
                        mt = Cs * Vt
                        Vs, Ms = max(Vs - Vt, 0.0), max(Ms - mt, 0.0)
                        m_true += mt
                        m_est += Cs_m * Vt
                        Fload_vol += Vt
                        det = Vt * dose
                        ww["capture_flowthrough"] += Vt + det
                        led.use("brij_10", det, "vi")
                        pl.det_L += det
                        if m_true > Q_sat:
                            loss["breakthrough"] += m_true - Q_sat
                            m_true = Q_sat
                            alarm(h, "pcc_breakthrough", "Lead column saturated; product lost in flow-through")
                        rec_timers = tick(stp)
                        budget -= stp
                        if Vs <= 1e-9:
                            rec_timers = tick(budget)
                            break
                    if waited:
                        wait_h += 1
                    col = 2 if rec_timers else (1 if Fload_vol > 0 else 0)
                    sump.L, sump.g = Vs, Ms

                # ── low-pH VI (legacy option): Protein A vessel -> VI pool ──
                cap_done = tail and sump.g <= 1e-3 and (not dsp.pa.busy() if batch else m_true <= 1e-6)
                if lowph:
                    pav = dsp.pav
                    F_vi_max = vi.tube_V_L / vi.t_hold_min_h / (1 + acid)
                    heel = 0.0 if cap_done else vi.heel_L          # run-out: drain the heel too
                    F_vi = _clamp((pav.L - heel) / ct.eluate_tau_h, 0.0, F_vi_max)
                    if cap_done and pav.L <= 2 * F_vi:
                        F_vi = pav.L
                    V_draw, m_draw = pav.take(F_vi) if F_vi > 1e-9 else (0.0, 0.0)
                    v_out = m_out = 0.0
                    if V_draw > 1e-9:
                        vi_on = True
                        led.use("vi_acid", V_draw * acid, "vi")
                        pH_in = vi.eluate_pH_mean + vi.eluate_pH_sd * rng.normal()
                        pH_out = pH_in - acid * vi.ph_drop_per_acid_frac
                        pH_meas = pH_out + pat.ph_noise * rng.normal()
                        ok = pH_out <= vi.pH_max
                        if not ok:
                            alarm(h, "vi_ph", f"VI pH {pH_out:.2f} above limit; material will be diverted")
                        vol_in = V_draw * (1 + acid)
                        tube.append([h, vol_in, m_draw, ok])
                        acid = acid_pi.step(pH_meas - vi.pH_target)
                    else:
                        vol_in = 0.0
                        if cap_done and pav.L <= 1e-4:
                            # run-out: chase the held material out of the tube with WFI
                            ready = 0.0
                            for pk in tube:
                                if pk[2] > 1e-12 and h - pk[0] < vi.t_hold_min_h - 1e-9:
                                    break
                                ready += pk[1]
                            if any(pk[2] > 1e-12 for pk in tube) and ready > 1e-9:
                                vol_in = min(ready, F_vi_max)
                                vi_on = True
                                led.use("wfi", vol_in, "vi")
                                ww["vi_flush"] += vol_in      # tube drained at the end of the campaign
                                tube.append([h, vol_in, 0.0, True])
                    if vol_in > 1e-9:
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
                    if v_out > 1e-9:
                        led.use("vi_base", v_out * vi.base_frac, "vi")
                        m1 = m_out * vi.vi_yield
                        loss["vi"] += m_out - m1
                        dsp.vip.add(v_out * (1 + vi.base_frac), m1, h + 1.0)
                    if V_draw > 1e-9 or vol_in > 1e-9:
                        if vi_t0 is None:
                            vi_t0, vi_L = float(h), 0.0
                        vi_L += V_draw or vol_in
                    elif vi_t0 is not None:
                        led.log("vi", "Viral inactivation", "Low-pH hold", "vi_acid", "hold", vi_t0, float(h), vi_L)
                        vi_t0 = None

                # ── step-by-step downstream ──
                dsp.upstream_quiet = ((batch or (sump.g <= 1e-3 and m_true <= 1e-6))
                                      and not (lowph and sum(pk[2] for pk in tube) > 1e-3))
                dsp.hour(h, tail)

                if dsp.pav.L > vi.eluate_surge_max_L:
                    alarm(h, "eluate_surge", "Protein A vessel SUM above capacity")
                if sump.L > sump_max:
                    over = sump.L - sump_max
                    c = sump.g / sump.L
                    ww["sump_overflow"] += over
                    loss["sump_overflow"] += over * c
                    sump.L -= over
                    sump.g -= over * c
                    alarm(h, "sump_overflow", "Harvest sump overflowing; capture cannot keep up")
                sump.peak_L = max(sump.peak_L, sump.L)
                if batch and capturing and Fh > 0 and sump.L / Fh > cp.max_sump_hold_h:
                    alarm(h, "sump_hold", f"Mean product age in the sump about {sump.L / Fh:.0f} h, "
                          f"above the {cp.max_sump_hold_h:.0f} h hold limit")
                if batch:
                    loading = any(tr.load is not None for tr in dsp.pa.trains)
                    if capturing and sump.g >= target and not loading:
                        wait_h += 1
                    col = 1 if loading else (2 if dsp.pa.busy() else 0)

            # ── record ──
            ww_h = sum(ww.values())
            for k, v in ww.items():
                ww_sum[k] += v
            mech = (not profile) and usp_on
            rec = {
                "h": h, "day": round(t_day, 3), "n_day": round(n_day, 3), "phase": int(in_n),
                "usp_on": int(usp_on),
                "X_n1": None if in_n else X1,
                "Xv": Xv if mech else None, "Xv_meas": Xv_m, "viability": viab,
                "S": S if mech else None, "sieve": sieve if mech else None,
                "titre_harvest": titre_h, "titre_meas": titre_m,
                "Fin_Lh": Fin, "Fb_Lh": Fb, "Fh_Lh": Fh, "VVD": Fin * 24.0 / V,
                "bleed_frac": Fb / Fin if Fin > 0 else 0.0, "bleed_on": int(bleed_on),
                "sump_L": sump.L, "sump_g": sump.g, "Fload_Lh": Fload_vol + pl.Fload, "cycles": pl.cycles,
                "capture_switch": pl.n_cyc_h, "col_state": col,
                "eluate_surge_L": dsp.pav.L, "F_vi_Lh": F_vi, "vi_pH": pH_out,
                "acid_frac": acid if lowph else 0.0, "det_Lh": pl.det_L,
                "harvest_g_h": harv_g, "ds_g_h": pl.flow["ds_g"], "ds_g_cum": pl.ds_g, "ds_L_cum": pl.ds_L,
                "buffer_L_h": sum(buf.values()) - sum(buf0.values()),
                "ww_total_Lh": ww_h,
                "ww_cells_Lh": ww["cell_bleed"] + ww["tff_retentate"] + ww["atf_flush"],
                "titre_fc_24h": None, "switch_fc_h": None, "ww_fc_24h_L": None,
                "g_fc_24h": None, "buffer_fc_24h_L": None,
            }
            rec.update(pl.flow)
            rec.update(dsp.levels())
            rec.update(dsp.lane_state(vi_on))
            for k in BUFFER_OPS:
                rec["buf_" + k] = buf[k] - buf0[k]
                rec["wfi_" + k] = led.wfi_h.get(k, 0.0)
            for k, v in ww.items():
                rec["ww_" + k] = v
            history.append(rec)

            # ── NevAI forecasts (advisory) ──
            if (self.fc is not None and usp_on and h % self.fc_every == 0
                    and len(history) >= self.fc.seq_len):
                f = self.fc.forecast_all(history)
                tf, wf = f.get("titre"), f.get("ww")
                if tf is not None:
                    rec["titre_fc_24h"] = float(tf[-1])
                    g_fc = float(np.sum(tf)) * Fh
                    rec["g_fc_24h"] = g_fc
                    rec["buffer_fc_24h_L"] = g_fc * per_g["total"] + dose * Fh * 24.0
                    have = sump.g if batch else m_est
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

        if vi_t0 is not None:
            led.log("vi", "Viral inactivation", "Low-pH hold", "vi_acid", "hold", vi_t0, float(H), vi_L)

        b = np.concatenate([[0.0], np.cumsum([r["buffer_L_h"] for r in history])])
        g = np.concatenate([[0.0], np.cumsum([r["harvest_g_h"] for r in history])])
        n = len(history)
        for i, r in enumerate(history):
            ok = r["usp_on"] and i + 24 <= n
            r["buffer_next24_L"] = float(b[i + 24] - b[i]) if ok else None
            r["g_next24"] = float(g[i + 24] - g[i]) if ok else None

        in_process_g = (dsp.in_process_g() + (m_true if not batch else 0.0) + sum(pk[2] for pk in tube)
                        + (P * V if not profile else 0.0))
        st = dict(pl=pl, dsp=dsp, produced_g=produced_g, media_L=media_L, loss=loss, ww_sum=ww_sum,
                  events=events, alarms=alarms, alarm_hours=alarm_hours, wait_h=wait_h,
                  atf_swaps=atf_swaps, fc_log=fc_log, cyc=cyc, per_g=per_g, n1s=n1s, in_process_g=in_process_g)
        return self._summarise(p, design, H, history, keep_history, downsample_h, st)

    def _summarise(self, p, design, H, history, keep_history, downsample_h, st):
        V = p.bioreactor.V_L
        pl, dsp = st["pl"], st["dsp"]
        led = pl.ledger
        ds_g, produced_g, media_L = pl.ds_g, st["produced_g"], st["media_L"]
        ww_sum, buf = st["ww_sum"], led.buf
        n_rows = [r for r in history if r["usp_on"]]
        mech = p.upstream.mode == "mechanistic"

        daily = []
        for d in range(int(math.ceil(H / 24))):
            rows = history[d * 24:(d + 1) * 24]
            if not rows:
                break
            if not any(r["phase"] for r in rows):
                ph = "N-1"
            elif any(r["usp_on"] for r in rows):
                ph = "N"
            else:
                ph = "run-out"
            by = {k: sum(r["ww_" + k] for r in rows) for k in WASTE_STREAMS}
            row = {"day": d + 1, "phase": ph, "total_L": round(sum(by.values()), 1),
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
        dsd = [d for d in daily if d["phase"] != "N-1"] or daily
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
        bdays = [d["buffer_L"] for d in dsd if d["buffer_L"] > 0]
        ch = np.array(pl.cycle_hours, dtype=float)
        errs = []
        for h_made, h_pred in st["fc_log"]:
            nxt = ch[ch > h_made]
            if nxt.size:
                errs.append(abs(h_pred - float(nxt[0])))
        mape = [abs(r["buffer_fc_24h_L"] - r["buffer_next24_L"]) / r["buffer_next24_L"]
                for r in history if r.get("buffer_fc_24h_L") is not None
                and r.get("buffer_next24_L") is not None and r["buffer_next24_L"] > 1.0]
        fcyc = pl.first_cycle or {}
        last_t = next((r["titre_harvest"] for r in reversed(history) if r["titre_harvest"] is not None), None)
        if mech:
            steady = [r for r in n_rows if r["bleed_on"]] or n_rows
            s_vcd = round(float(np.mean([r["Xv"] for r in steady])), 1) if steady else None
            s_via = round(100 * float(np.mean([r["viability"] for r in steady])), 1) if steady else None
            s_tit = round(float(np.mean([r["titre_harvest"] for r in steady])), 3) if steady else None
        else:
            s_vcd = s_via = s_tit = None
        rep = dsp.report()

        kpi = {
            "drug_substance_kg": round(ds_g / 1000, 3),
            "ds_volume_L": round(pl.ds_L, 1),
            "product_made_kg": round(produced_g / 1000, 3),
            "overall_yield_pct": round(100 * ds_g / produced_g, 1) if produced_g > 0 else 0.0,
            "volumetric_productivity_g_L_day": round(ds_g / V / p.run_days, 3),
            "n1_seed_vcd": st["n1s"]["n_seed_vcd"],
            "n1_waste_m3": round(ww_sum["n1_permeate"] / 1000, 3),
            "media_used_m3": round(media_L / 1000, 2),
            "media_m3_per_kg_ds": round(media_L / ds_g, 2) if ds_g > 0 else None,
            "final_titre_gL": round(last_t, 3) if last_t is not None else None,
            "steady_vcd_mean": s_vcd, "steady_viability_pct": s_via, "steady_titre_gL": s_tit,
            "capture_cycles": pl.cycles,
            "first_cycle_n_day": fcyc.get("n_day"),
            "first_cycle_sump_L": fcyc.get("sump_L"),
            "capture_mean_interval_h": round(float(np.mean(np.diff(ch))), 2) if ch.size > 1 else None,
            "capture_recovery_h": round(st["cyc"]["recovery_h"], 2),
            "capture_wait_h": st["wait_h"],
            "peak_cycles_day": max((d["cycles"] for d in dsd), default=0),
            "peak_g_day": max((d["harvest_g"] for d in nd), default=0.0),
            "buffer_m3": round(buf_total / 1000, 3),
            "buffer_L_per_g_ds": round(buf_total / ds_g, 3) if ds_g > 0 else None,
            "buffer_daily_peak_L": max(bdays) if bdays else 0.0,
            "buffer_peak_to_avg": round(max(bdays) / float(np.mean(bdays)), 2) if bdays else None,
            "detergent_L": round(led.sol.get("brij_10", 0.0), 1),
            "wfi_m3": round(led.sol.get("wfi", 0.0) / 1000, 3),
            "depth_filters": pl.depth_filters,
            "vf_changes": pl.vf_changes,
            "lots": len(pl.lots),
            "containers": pl.containers,
            "in_process_end_g": round(st["in_process_g"], 1),
            "atf_swaps": st["atf_swaps"],
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
        labels = {k: v["label"] for k, v in SOLUTIONS.items()}
        labels.update(SOLUTION_ALIASES)
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
            "lots": [{"lot": L["lot"], "day": L["day"], "g": round(L["g"], 1), "L": round(L["L"], 2),
                      "containers": L["containers"]} for L in pl.lots],
            "solutions": rep["solutions"], "chemicals": rep["chemicals"], "vessels": rep["vessels"],
            "steps": rep["steps"], "lanes": rep["lanes"], "solution_labels": labels,
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
    print("Solutions (L):")
    for s in res["solutions"]:
        print(f"  {s['L']:>10.1f}  {s['label']}")
    print(f"Steps logged: {len(res['steps'])}")
