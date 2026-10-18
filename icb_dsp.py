"""
icb_dsp.py — step-by-step downstream (Case 5 flowsheet) for the IDB digital twin
NavNevAI_mAbWasteWater · CarbonTesseract 2026

Each unit runs recipes: ordered steps, each with a solution, a volume, a time and
a destination (product, a waste stream, or held in the equipment). Columns run
start-up steps before first use and storage steps at the end of the campaign.
Vessels sit between units; a unit starts when its trigger is met.

  Protein A  2 trains: load + in-line Brij-9 -> post-load wash -> wash 1 -> wash 2
             -> elute (peak to Protein A vessel) -> strip -> CIP -> re-equilibrate
  Depth      new filter per 3 eluates: WFI flush -> equilibrate -> filter -> chase
             -> salt, arginine, caustic, Tris flushes (order assumed)
  AEX        2 columns: in-line Tris base + NaCl -> load (flow-through kept) -> chase
             -> strip -> CIP -> re-equilibrate
  CEX        2 columns: in-line acetic acid -> load -> wash -> elute (peak kept)
             -> strip -> CIP -> re-equilibrate
  VRF        2 holders; per new filter: WFI flush, equilibrate; per batch: filter -> chase
  UF/DF      per lot: flush -> equilibrate -> transfer in -> UF1 -> DF -> UF2 -> recovery
             -> formulation spike -> transfer to BDS mixer -> flush -> CIP -> rinse -> storage
  BDS        mix -> flush 0.2 um filter -> fill cryovault containers

Each hour is split into SLICES so short steps keep their own start and end times
in the step log (result["steps"]) used by the diagram and the films.
"""
from __future__ import annotations
import math
from icb_process import BUFFER_OPS, DILUTIONS, SOLUTIONS, capture_cycle, det_frac

EPS_G, EPS_L = 1e-3, 1e-4
MAX_LOG = 20000
SLICES = 4
FLOW_KEYS = ("L_vi", "g_vi", "L_dep", "g_dep", "L_base", "L_acid",
             "L_aex_in", "g_aex_in", "L_aex", "g_aex", "L_cex_in", "g_cex_in", "L_cex", "g_cex",
             "L_vf_in", "L_vf", "g_vf", "L_uf_in", "g_uf_in", "L_df", "L_exc", "L_ret", "g_ret",
             "L_ds", "ds_g", "elu_pa1", "elu_pa2", "elu_cex1", "elu_cex2")
VESSEL_KEYS = ("pav_L", "vip_L", "dfs_L", "aexv_L", "cexv_L", "vrfv_L", "ufdf_L", "bdsm_L")
LANES = (("det", "Brij-9 dosing"), ("pa1", "Protein A col 1 & 2"), ("pa2", "Protein A col 3 & 4"),
         ("vi", "Low-pH VI hold"), ("depth", "Depth filtration"), ("base", "Tris base + NaCl in-line"),
         ("aex1", "AEX column 1"), ("aex2", "AEX column 2"), ("acid", "Acetic acid in-line"),
         ("cex1", "CEX column 1"), ("cex2", "CEX column 2"), ("vrf1", "Virus filter holder A"),
         ("vrf2", "Virus filter holder B"), ("ufdf", "UF/DF"), ("bds", "BDS mix and fill"))
LANE_KEYS = tuple("st_" + k for k, _ in LANES) + tuple("sk_" + k for k, _ in LANES)
WFI_KEYS = tuple("wfi_" + k for k in BUFFER_OPS)


class Vessel:
    def __init__(self, key, label):
        self.key, self.label = key, label
        self.L = self.g = self.peak_L = 0.0
        self.since = None                  # when it last went from empty to holding material

    def add(self, L, g, t):
        L, g = max(L, 0.0), max(g, 0.0)
        if L <= 0.0 and g <= 0.0:
            return
        if self.L <= EPS_L:
            self.since = t
        self.L += L
        self.g += g
        self.peak_L = max(self.peak_L, self.L)

    def take(self, L):
        if self.L <= 0.0:
            return 0.0, 0.0
        L = min(max(L, 0.0), self.L)
        g = self.g * L / self.L
        self.L -= L
        self.g -= g
        if self.L <= EPS_L:                # sweep the residue with the last draw
            L += self.L
            g += self.g
            self.L = self.g = 0.0
            self.since = None
        return L, max(g, 0.0)


class Ledger:
    def __init__(self):
        self.sol = {k: 0.0 for k in SOLUTIONS}
        self.buf = {k: 0.0 for k in BUFFER_OPS}
        self.ww = {}
        self.wfi_h = dict.fromkeys(BUFFER_OPS, 0.0)
        self.steps = []

    def use(self, sol, L, op):
        if not sol or L <= 0.0:
            return
        for k, f in DILUTIONS.get(sol, ((sol, 1.0),)):
            self.sol[k] = self.sol.get(k, 0.0) + L * f
            if k == "wfi":
                self.wfi_h[op] = self.wfi_h.get(op, 0.0) + L * f
        self.buf[op] = self.buf.get(op, 0.0) + L

    def waste(self, stream, L):
        if stream and L > 0.0:
            self.ww[stream] = self.ww.get(stream, 0.0) + L

    def log(self, lane, unit, step, sol, kind, t0, t1, L):
        if len(self.steps) < MAX_LOG:
            self.steps.append({"lane": lane, "unit": unit, "step": step, "sol": sol, "kind": kind,
                               "t0": round(t0, 4), "t1": round(max(t1, t0), 4), "L": round(L, 3)})


class Plant:
    """State shared by the units (set up by icb_dynamic.ICBSimulator)."""
    def __init__(self, p, rng, event, alarm, loss):
        self.p, self.rng, self.event, self.alarm, self.loss = p, rng, event, alarm, loss
        self.ledger = Ledger()
        self.flow = dict.fromkeys(FLOW_KEYS, 0.0)
        self.h, self.n_day, self.capturing, self.perm_L = 0, 0.0, False, 0.0
        self.cycles, self.n_cyc_h, self.cycle_hours, self.first_cycle = 0, 0, [], None
        self.Fload = self.det_L = 0.0
        self.depth_filters = self.vf_changes = self.containers = 0
        self.lots, self.pending_lot = [], None
        self.ds_g = self.ds_L = 0.0

    def lose(self, key, g):
        if g > 0.0:
            self.loss[key] = self.loss.get(key, 0.0) + g


def step(name, sol, L, flow_Lh, route, kind="buf", contact_h=0.0, collect_L=0.0):
    L = max(float(L), 0.0)
    h = (L / flow_Lh if flow_Lh > 0 else 0.0) + max(float(contact_h), 0.0)
    return {"name": name, "sol": sol, "L": L, "h": h, "route": route, "kind": kind,
            "collect_L": min(max(collect_L, 0.0), L)}


def load_step(name):
    return {"name": name, "sol": None, "L": 0.0, "h": 0.0, "route": None, "kind": "load", "collect_L": 0.0}


def hold_step(name, h):
    return {"name": name, "sol": None, "L": 0.0, "h": max(float(h), 0.0), "route": None,
            "kind": "hold", "collect_L": 0.0}


def _keep(steps):
    """Drop steps with no volume and no time (setting a CV to 0 removes the step)."""
    return [s for s in steps if s["kind"] in ("load", "out") or s["L"] > 1e-9 or s["h"] > 1e-9]


class Train:
    def __init__(self, lane):
        self.lane = lane
        self.job, self.i, self.t, self.t0 = [], 0, 0.0, 0.0
        self.fresh, self.started, self.stored = True, False, False
        self.load = self.snap = self.plan = None
        self.bound = self.elute_g = self.elute_loss = 0.0
        self.L_step = self.out_L = 0.0
        self.lot, self.lot_g, self.lot_L = 0, 0.0, 0.0

    @property
    def busy(self):
        return self.i < len(self.job)

    @property
    def step(self):
        return self.job[self.i] if self.i < len(self.job) else None

    def committed(self):
        """Busy and the load step is still ahead or running."""
        return self.busy and any(s["kind"] == "load" for s in self.job[self.i:])


class Unit:
    loss_key = "polishing"

    def __init__(self, plant, key, label, op, feed, out, n_trains=1, numbered=False):
        self.pl, self.key, self.label, self.op = plant, key, label, op
        self.feed, self.out = feed, out
        n = max(1, int(n_trains))
        self.trains = [Train(f"{key}{k + 1}" if numbered else key) for k in range(n)]
        self.next = 0
        self.force = self.finish = False
        self.last_out_L = 0.0

    # recipes and behaviour (overridden)
    def startup(self): return []
    def cycle(self, tr): return []
    def storage(self): return []
    def trigger(self, t): return False
    def job_start(self, tr, t): pass
    def begin_load(self, tr, t): tr.load = {"left": self.feed.L, "add": 0.0}
    def load_rate(self, tr): return 1.0
    def take_load(self, tr, L, g, t): tr.L_step += L
    def end_load(self, tr, t): pass
    def special(self, tr, s, frac, t): pass
    def step_done(self, tr, s, t): pass

    def job_done(self, tr, t):
        if tr.out_L > 1e-9:
            self.last_out_L = tr.out_L

    def collect(self, tr, L, g, t):
        if self.out is not None:
            self.out.add(L, g, t)
        tr.out_L += L

    def busy(self):
        return any(tr.busy for tr in self.trains)

    def advance(self, t0, span):
        n = len(self.trains)
        for tr in [self.trains[(self.next + k) % n] for k in range(n)]:
            self._run(tr, t0, span)

    def _begin(self, tr, job, t):
        tr.job, tr.i, tr.t, tr.fresh = job, 0, 0.0, True
        tr.load, tr.snap, tr.out_L = None, None, 0.0
        self.job_start(tr, t)

    def _start(self, tr, t):
        if self.finish:
            if tr.started and not tr.stored:
                tr.stored = True
                job = _keep(self.storage())
                if job:
                    self._begin(tr, job, t)
                    return True
            return False
        if any(o is not tr and o.committed() for o in self.trains):
            return False
        if not self.trigger(t):
            return False
        job = _keep((self.startup() if not tr.started else []) + self.cycle(tr))
        tr.started = True
        self._begin(tr, job, t)
        self.next = (self.trains.index(tr) + 1) % len(self.trains)
        return True

    def _apply(self, tr, s, frac, t):
        k = s["kind"]
        if k in ("buf", "elute"):
            led = self.pl.ledger
            L = s["L"] * frac
            led.use(s["sol"], L, self.op)
            tr.L_step += L
            if k == "elute":
                cL = s["collect_L"] * frac
                self.collect(tr, cL, tr.elute_g * frac, t)
                self.pl.lose(self.loss_key, tr.elute_loss * frac)
                led.waste(s["route"], L - cL)
            elif s["route"] == "product":
                self.collect(tr, L, 0.0, t)
            else:
                led.waste(s["route"], L)
        elif k != "hold":
            self.special(tr, s, frac, t)

    def _done(self, tr, t):
        s = tr.step
        self.step_done(tr, s, t)
        self.pl.ledger.log(tr.lane, self.label, s["name"], s["sol"], s["kind"], tr.t0, t, tr.L_step)
        tr.i += 1
        tr.t, tr.fresh, tr.snap = 0.0, True, None
        if not tr.busy:
            self.job_done(tr, t)
            tr.job, tr.i = [], 0

    def _end_load(self, tr, t):
        self.end_load(tr, t)
        tr.load = None
        self._done(tr, t)

    def _run(self, tr, t0, span):
        left = span
        for _ in range(500):
            if left <= 1e-9:
                return
            t = t0 + span - left
            if not tr.busy and not self._start(tr, t):
                return
            s = tr.step
            if tr.fresh:
                tr.t0, tr.fresh, tr.L_step = t, False, 0.0
            if s["kind"] == "load":
                if tr.load is None:
                    if any(o is not tr and o.load is not None for o in self.trains):
                        return                              # one load pump per unit
                    tr.t0 = t
                    self.begin_load(tr, t)
                rate = self.load_rate(tr)
                need = tr.load["left"]
                if need <= 1e-9 or self.feed.L <= EPS_L or rate <= 0.0:
                    self._end_load(tr, t)
                    continue
                dt = min(left, need / rate)
                L, g = self.feed.take(min(rate * dt, need))
                tr.load["left"] -= L
                self.take_load(tr, L, g, t)
                used = min(dt, L / rate)
                left -= used
                if tr.load["left"] <= 1e-6 or self.feed.L <= EPS_L:
                    self._end_load(tr, t0 + span - left)
                elif used <= 1e-12:
                    return
            else:
                h = s["h"]
                if h <= 1e-12:
                    self._apply(tr, s, 1.0, t)
                    self._done(tr, t)
                    continue
                dt = min(left, h - tr.t)
                self._apply(tr, s, dt / h, t)
                tr.t += dt
                left -= dt
                if tr.t >= h - 1e-9:
                    self._done(tr, t0 + span - left)


class ProteinA(Unit):
    loss_key = "capture"

    def __init__(self, plant, sump, pav):
        p = plant.p
        super().__init__(plant, "pa", "Protein A", "capture", sump, pav, p.capture.trains, numbered=True)
        self.cp, self.cyc = p.capture, capture_cycle(p)
        self.CV, self.dose = self.cyc["CV_L"], det_frac(p)
        self.f = self.cp.recovery_flow_CVh * self.CV
        self._eps_h, self._eps = None, 0.0

    def startup(self):
        c, V, f = self.cp, self.CV, self.f
        return [step("WFI flush", "wfi", c.startup_wfi_CV * V, f, "unit_flushes"),
                step("Sanitise 0.5 N NaOH", "naoh_05n", c.startup_naoh_CV * V, f, "sanitisation",
                     contact_h=c.cip_contact_h),
                step("Equilibrate", "tris_74", c.equil_CV * V, f, "unit_flushes")]

    def cycle(self, tr):
        c, V, f = self.cp, self.CV, self.f
        return [load_step("Load + Brij-9" if self.dose > 0 else "Load"),
                step("Post-load wash", "tris_74", c.postload_CV * V, f, "capture_buffers"),
                step("Wash 1 pH 9.0", "pa_wash1", c.wash1_CV * V, f, "capture_buffers"),
                step("Wash 2 pH 6.5", "pa_wash2", c.wash2_CV * V, f, "capture_buffers"),
                step("Elute pH 4.45", "pa_elute", c.elute_CV * V, f, "capture_buffers", kind="elute",
                     collect_L=c.eluate_CV * V),
                step("Strip acetic acid", "pa_strip", c.strip_CV * V, f, "capture_strip_cip"),
                step("CIP 0.5 N NaOH", "naoh_05n", c.cip_CV * V, f, "capture_strip_cip", contact_h=c.cip_contact_h),
                step("Re-equilibrate", "tris_74", c.reequil_CV * V, f, "capture_buffers")]

    def storage(self):
        c, V, f = self.cp, self.CV, self.f
        return [step("WFI flush", "wfi", c.storage_wfi_CV * V, f, "unit_flushes"),
                step("Store benzyl alcohol", "storage", c.storage_CV * V, f, None)]

    def _noisy(self):
        """PAT titre reading of the surge tank: one measurement per hour."""
        f, pl = self.feed, self.pl
        if f.L <= EPS_L:
            return 0.0
        if self._eps_h != pl.h:
            self._eps_h, self._eps = pl.h, pl.rng.normal()
        return max(0.0, f.g / f.L * (1 + pl.p.pat.titre_noise * self._eps))

    def trigger(self, t):
        if not self.pl.capturing:
            return False
        if self.force:
            return self.feed.g > EPS_G
        c = self._noisy()
        return c > 1e-9 and c * self.feed.L >= self.cyc["target_g"]

    def begin_load(self, tr, t):
        pl, f = self.pl, self.feed
        c = 0.0 if self.force else self._noisy()
        left = min(f.L, self.cyc["target_g"] / c) if c > 1e-9 else f.L
        tr.load = {"left": left, "det": 0.0}
        tr.bound = 0.0
        if pl.first_cycle is None:
            pl.first_cycle = {"n_day": round(pl.n_day + (t - pl.h) / 24.0, 2), "sump_L": round(pl.perm_L, 1)}
            pl.event(pl.h, "first_cycle", f"First Protein A cycle after {pl.perm_L:.0f} L of permeate "
                                          f"(N day {pl.n_day:.1f})")

    def load_rate(self, tr):
        return self.cyc["load_flow_Lh"] / (1.0 + self.dose)

    def take_load(self, tr, L, g, t):
        pl, det = self.pl, L * self.dose
        tr.bound += g
        tr.load["det"] += det
        tr.L_step += L + det
        pl.ledger.use("brij_10", det, "vi")
        pl.ledger.waste("capture_flowthrough", L + det)
        pl.Fload += L
        pl.det_L += det

    def end_load(self, tr, t):
        pl, cap, y = self.pl, self.cyc["capacity_g"], self.cyc["yield"]
        if tr.bound > cap:
            pl.lose("breakthrough", tr.bound - cap)
            pl.alarm(pl.h, "capture_overload", "Column loaded above capacity; product lost in flow-through")
            tr.bound = cap
        tr.elute_g, tr.elute_loss = tr.bound * y, tr.bound * (1 - y)
        pl.cycles += 1
        pl.n_cyc_h += 1
        pl.cycle_hours.append(pl.h)
        if pl.cycles == self.cyc["resin_cycles_max"]:
            pl.event(pl.h, "resin_life", "Protein A resin reached its cycle limit; plan replacement")
        if tr.load["det"] > 0:
            pl.ledger.log("det", "Detergent", "Dose Brij-9", "brij_10", "buf", tr.t0, t, tr.load["det"])

    def collect(self, tr, L, g, t):
        super().collect(tr, L, g, t)
        k = "elu_" + tr.lane
        self.pl.flow[k] = self.pl.flow.get(k, 0.0) + L


class DepthFilter(Unit):
    loss_key = "depth"

    def __init__(self, plant, feed, out, trigger_L):
        super().__init__(plant, "depth", "Depth filtration", "depth", feed, out, 1)
        self.d, self.trigger_L = plant.p.depth, trigger_L

    def cycle(self, tr):
        d = self.d
        A, f = d.area_m2, d.flush_flux_LMH * d.area_m2
        return [step("WFI flush (new filter)", "wfi", d.wfi_flush_Lm2 * A, f, "depth_flush"),
                step("Equilibrate pH 5.0", "df_ac50", d.equil_Lm2 * A, f, "depth_flush"),
                load_step("Filter pool"),
                step("Chase pH 5.0 (kept)", "df_ac50", d.chase_Lm2 * A, f, "product"),
                step("Flush 0.5 M NaCl", "nacl_05m", d.salt_Lm2 * A, f, "depth_post"),
                step("Flush arginine", "df_arg", d.arginine_Lm2 * A, f, "depth_post"),
                step("Flush 1 N NaOH", "naoh_1n", d.naoh_Lm2 * A, f, "depth_post", contact_h=d.naoh_contact_h),
                step("Flush Tris pH 7.4", "tris_74", d.tris_Lm2 * A, f, "depth_post")]

    def trigger(self, t):
        if self.force:
            return self.feed.g > EPS_G
        return self.feed.g > EPS_G and self.feed.L >= self.trigger_L - 1e-6

    def job_start(self, tr, t):
        self.pl.depth_filters += 1

    def load_rate(self, tr):
        return self.d.flux_LMH * self.d.area_m2

    def take_load(self, tr, L, g, t):
        y = self.d.depth_yield
        tr.L_step += L
        self.pl.flow["L_vi"] += L
        self.pl.flow["g_vi"] += g
        self.collect(tr, L, g * y, t)
        self.pl.lose("depth", g * (1 - y))

    def collect(self, tr, L, g, t):
        super().collect(tr, L, g, t)
        self.pl.flow["L_dep"] += L
        self.pl.flow["g_dep"] += g


class AEX(Unit):
    def __init__(self, plant, feed, out):
        p = plant.p
        super().__init__(plant, "aex", "AEX", "aex", feed, out, p.aex.columns, numbered=True)
        self.a, self.y = p.aex, p.polish.aex_yield
        self.cap = self.a.CV_L * self.a.load_gL
        self.add = self.a.base_frac + self.a.nacl_frac
        self.f = self.a.flow_CVh * self.a.CV_L

    def startup(self):
        a, V, f = self.a, self.a.CV_L, self.f
        return [step("WFI flush", "wfi", a.startup_wfi_CV * V, f, "unit_flushes"),
                step("Sanitise 1 N NaOH", "naoh_1n", a.startup_naoh_CV * V, f, "sanitisation",
                     contact_h=a.cip_contact_h),
                step("Equilibrate", "tris_74", a.equil_CV * V, f, "unit_flushes")]

    def cycle(self, tr):
        a, V, f = self.a, self.a.CV_L, self.f
        return [load_step("Load (flow-through kept)"),
                step("Chase (kept)", "tris_74", a.chase_CV * V, f, "product"),
                step("Strip 0.5 M NaCl", "tris_nacl_74", a.strip_CV * V, f, "aex_buffers"),
                step("CIP 1 N NaOH", "naoh_1n", a.cip_CV * V, f, "aex_cip", contact_h=a.cip_contact_h),
                step("Re-equilibrate", "tris_74", a.reequil_CV * V, f, "aex_buffers")]

    def storage(self):
        a, V, f = self.a, self.a.CV_L, self.f
        return [step("WFI flush", "wfi", a.storage_wfi_CV * V, f, "unit_flushes"),
                step("Store benzyl alcohol", "storage", a.storage_CV * V, f, None)]

    def trigger(self, t):
        if self.force:
            return self.feed.g > EPS_G
        return self.feed.g >= self.a.min_load_frac * self.cap

    def begin_load(self, tr, t):
        f = self.feed
        m = min(f.g, self.cap)
        tr.load = {"left": f.L * m / f.g if f.g > 1e-12 else f.L, "add": 0.0}

    def load_rate(self, tr):
        return self.f / (1.0 + self.add)

    def take_load(self, tr, L, g, t):
        a, pl = self.a, self.pl
        pl.ledger.use("tris_base_2m", L * a.base_frac, "adjust")
        pl.ledger.use("nacl_2m", L * a.nacl_frac, "adjust")
        fl = pl.flow
        fl["L_aex_in"] += L
        fl["g_aex_in"] += g
        fl["L_base"] += L * self.add
        tr.load["add"] += L * self.add
        tr.L_step += L * (1 + self.add)
        self.collect(tr, L * (1 + self.add), g * self.y, t)
        pl.lose("polishing", g * (1 - self.y))

    def end_load(self, tr, t):
        if tr.load["add"] > 0:
            self.pl.ledger.log("base", "In-line addition", "Tris base + NaCl", "tris_base_2m", "buf",
                               tr.t0, t, tr.load["add"])

    def collect(self, tr, L, g, t):
        super().collect(tr, L, g, t)
        self.pl.flow["L_aex"] += L
        self.pl.flow["g_aex"] += g


class CEX(Unit):
    def __init__(self, plant, feed, out):
        p = plant.p
        super().__init__(plant, "cex", "CEX", "cex", feed, out, p.cex.columns, numbered=True)
        self.c, self.y = p.cex, p.polish.cex_yield
        self.cap = self.c.CV_L * self.c.DBC_gL
        self.f = self.c.flow_CVh * self.c.CV_L

    def startup(self):
        c, V, f = self.c, self.c.CV_L, self.f
        return [step("WFI flush", "wfi", c.startup_wfi_CV * V, f, "unit_flushes"),
                step("Sanitise 1 N NaOH", "naoh_1n", c.startup_naoh_CV * V, f, "sanitisation",
                     contact_h=c.cip_contact_h),
                step("Equilibrate pH 5.5", "cex_eq", c.equil_CV * V, f, "unit_flushes")]

    def cycle(self, tr):
        c, V, f = self.c, self.c.CV_L, self.f
        return [load_step("Load (bind)"),
                step("Wash pH 5.5", "cex_eq", c.wash_CV * V, f, "cex_buffers"),
                step("Elute 280 mM NaCl", "cex_elute", c.elute_CV * V, f, "cex_buffers", kind="elute",
                     collect_L=c.eluate_CV * V),
                step("Strip 0.5 M NaCl", "tris_nacl_74", c.strip_CV * V, f, "cex_buffers"),
                step("CIP 1 N NaOH", "naoh_1n", c.cip_CV * V, f, "cex_cip", contact_h=c.cip_contact_h),
                step("Re-equilibrate", "cex_eq", c.reequil_CV * V, f, "cex_buffers")]

    def storage(self):
        c, V, f = self.c, self.c.CV_L, self.f
        return [step("WFI flush", "wfi", c.storage_wfi_CV * V, f, "unit_flushes"),
                step("Store benzyl alcohol", "storage", c.storage_CV * V, f, None)]

    def trigger(self, t):
        if self.force:
            return self.feed.g > EPS_G
        return self.feed.g >= self.c.min_load_frac * self.cap

    def begin_load(self, tr, t):
        f = self.feed
        m = min(f.g, self.cap)
        tr.load = {"left": f.L * m / f.g if f.g > 1e-12 else f.L, "add": 0.0}
        tr.bound = 0.0

    def load_rate(self, tr):
        return self.f / (1.0 + self.c.acid_frac)

    def take_load(self, tr, L, g, t):
        a, pl = self.c.acid_frac, self.pl
        pl.ledger.use("acetic_1m", L * a, "adjust")
        pl.ledger.waste("cex_flowthrough", L * (1 + a))
        fl = pl.flow
        fl["L_cex_in"] += L
        fl["g_cex_in"] += g
        fl["L_acid"] += L * a
        tr.bound += g
        tr.load["add"] += L * a
        tr.L_step += L * (1 + a)

    def end_load(self, tr, t):
        tr.elute_g, tr.elute_loss = tr.bound * self.y, tr.bound * (1 - self.y)
        if tr.load["add"] > 0:
            self.pl.ledger.log("acid", "In-line addition", "Acetic acid", "acetic_1m", "buf",
                               tr.t0, t, tr.load["add"])

    def collect(self, tr, L, g, t):
        super().collect(tr, L, g, t)
        fl = self.pl.flow
        fl["L_cex"] += L
        fl["g_cex"] += g
        k = "elu_" + tr.lane
        fl[k] = fl.get(k, 0.0) + L


class VirusFilter(Unit):
    def __init__(self, plant, feed, out, prev):
        p = plant.p
        super().__init__(plant, "vrf", "Virus filtration", "vf", feed, out, 1)
        self.v, self.po, self.prev = p.vrf, p.polish, prev
        self.used, self.holder = 0.0, 1
        self.nominal_L = p.cex.CV_L * p.cex.eluate_CV
        self.trains[0].lane = "vrf1"

    def job_start(self, tr, t):
        tr.lane = f"vrf{self.holder}"

    def cycle(self, tr):
        v, po, f = self.v, self.po, self.v.flow_Lh
        pre = [] if self.used > 1e-9 else [
            step("WFI flush (new filter)", "wfi", po.vf_flush_L, f, "vf_flush"),
            step("Equilibrate", "cex_elute", v.equil_L, f, "vf_flush")]
        return pre + [load_step("Filter"), step("Chase (kept)", "cex_elute", v.chase_L, f, "product")]

    def trigger(self, t):
        if self.force:
            return self.feed.L > EPS_L
        unit = self.prev.last_out_L if self.prev.last_out_L > 1e-9 else self.nominal_L
        return self.feed.g > EPS_G and self.feed.L >= self.v.trigger_elutions * unit - 1e-6

    def load_rate(self, tr):
        return self.v.flow_Lh

    def take_load(self, tr, L, g, t):
        y = self.po.vf_yield
        self.used += L
        tr.L_step += L
        self.pl.flow["L_vf_in"] += L
        self.collect(tr, L, g * y, t)
        self.pl.lose("polishing", g * (1 - y))

    def collect(self, tr, L, g, t):
        super().collect(tr, L, g, t)
        self.pl.flow["L_vf"] += L
        self.pl.flow["g_vf"] += g

    def job_done(self, tr, t):
        super().job_done(tr, t)
        if self.used >= self.po.vf_capacity_L:
            self.used, self.holder = 0.0, (2 if self.holder == 1 else 1)
            self.pl.vf_changes += 1
            self.pl.event(self.pl.h, "vf_change", "Virus filter changed (capacity reached); next batch uses the other holder")


class UFDF(Unit):
    loss_key = "sptff"

    def __init__(self, plant, feed, out, tank):
        super().__init__(plant, "ufdf", "UF/DF", "sptff", feed, out, 1)
        self.s, self.tank, self.lot = plant.p.sptff, tank, 0

    def trigger(self, t):
        f = self.feed
        if f.g <= EPS_G:
            return False
        if self.force:
            return True
        return (f.since is not None and t - f.since >= self.pl.p.lot_window_h - 1e-9
                and f.g >= self.s.min_lot_g)

    def job_start(self, tr, t):
        self.lot += 1
        tr.lot = self.lot

    def cycle(self, tr):
        s = self.s
        V, M = self.feed.L, self.feed.g
        fs = min(max(s.excipient_frac, 0.0), 0.9)
        Vfin = M * s.sptff_yield / max(s.target_gL, 1e-9)
        spike = fs * Vfin
        V2 = Vfin - spike - s.recover_L
        if V2 < s.min_V_L:
            V2 = s.min_V_L
            Vfin = (V2 + s.recover_L) / (1.0 - fs)
            spike = fs * Vfin
        V1 = min(max(M / max(s.uf1_gL, 1e-9), V2), V)
        V2 = min(V2, V1)
        perm1, df_L, perm2 = max(0.0, V - V1), s.DV * V1, max(0.0, V1 - V2)
        Vout = V2 + s.recover_L + spike
        tr.plan = {"V": V}
        F, Fx = s.perm_flux_Lh, s.flush_Lh
        return [step("WFI flush", "wfi", s.wfi_flush_L, Fx, "ufdf_cip"),
                step("Equilibrate", "uf_df1", s.equil_L, Fx, "ufdf_cip"),
                load_step("Transfer in"),
                dict(step("UF1 concentrate", None, perm1, F, "sptff_permeate"), kind="perm"),
                dict(step("Diafilter pH 5.1", "uf_df1", df_L, F, "sptff_permeate"), kind="df"),
                dict(step("UF2 concentrate", None, perm2, F, "sptff_permeate"), kind="perm"),
                step("Recovery flush pH 5.3", "uf_df2", s.recover_L, Fx, "product"),
                step("Formulation spike", "uf_spike", spike, Fx, "product"),
                dict(step("Transfer to BDS mixer", None, Vout, s.transfer_Lh, None), kind="out"),
                step("WFI flush", "wfi", s.wfi_flush_L, Fx, "ufdf_cip"),
                step("CIP 1 N NaOH", "naoh_1n", s.cip_L, Fx, "ufdf_cip", contact_h=s.cip_contact_h),
                step("WFI rinse", "wfi", s.wfi_flush_L, Fx, "ufdf_cip"),
                step("Store 0.1 N NaOH", "naoh_01n", s.storage_L, Fx, None)]

    def begin_load(self, tr, t):
        tr.load = {"left": min(tr.plan["V"], self.feed.L), "add": 0.0}

    def load_rate(self, tr):
        return self.s.transfer_Lh

    def take_load(self, tr, L, g, t):
        self.tank.add(L, g, t)
        tr.L_step += L
        self.pl.flow["L_uf_in"] += L
        self.pl.flow["g_uf_in"] += g

    def end_load(self, tr, t):
        if self.feed.L > EPS_L:
            self.feed.since = t            # leftovers start the next lot window

    def collect(self, tr, L, g, t):
        self.tank.add(L, g, t)

    def _move(self, tr, L, g, t):
        L, g = min(L, self.tank.L), min(g, self.tank.g)
        self.tank.L -= L
        self.tank.g -= g
        y = self.s.sptff_yield
        self.out.add(L, g * y, t)
        self.pl.lose("sptff", g * (1 - y))
        self.pl.flow["L_ret"] += L
        self.pl.flow["g_ret"] += g * y
        tr.L_step += L

    def special(self, tr, s, frac, t):
        led, L, k = self.pl.ledger, s["L"] * frac, s["kind"]
        if k == "perm":
            self.tank.L = max(self.tank.L - L, 0.0)
            led.waste(s["route"], L)
            tr.L_step += L
        elif k == "df":
            led.use(s["sol"], L, self.op)
            led.waste(s["route"], L)
            tr.L_step += L
            self.pl.flow["L_df"] += L
        elif k == "out":
            if tr.snap is None:
                tr.snap = (self.tank.L, self.tank.g)
            self._move(tr, tr.snap[0] * frac, tr.snap[1] * frac, t)

    def step_done(self, tr, s, t):
        if s["kind"] == "out":
            if self.tank.L > 1e-9 or self.tank.g > 1e-9:
                self._move(tr, self.tank.L, self.tank.g, t)
            self.tank.L = self.tank.g = 0.0
            self.pl.pending_lot = tr.lot
        elif s["sol"] == "uf_spike":
            self.pl.flow["L_exc"] += s["L"]


class BDS(Unit):
    loss_key = "sptff"

    def __init__(self, plant, feed, ufdf):
        super().__init__(plant, "bds", "BDS", "bds", feed, None, 1)
        self.b, self.ufdf = plant.p.bds, ufdf

    def trigger(self, t):
        if self.feed.g <= EPS_G:
            return False
        s = self.ufdf.trains[0].step
        return not (s is not None and s["kind"] == "out")

    def job_start(self, tr, t):
        tr.lot = self.pl.pending_lot or (len(self.pl.lots) + 1)
        tr.lot_g = tr.lot_L = 0.0

    def cycle(self, tr):
        b = self.b
        return [hold_step("Mix in BDS hold mixer", b.mix_h),
                step("Flush 0.2 um filter", "wfi", b.filter_flush_L, b.fill_Lh, "bds_flush"),
                load_step("Fill cryovault containers")]

    def load_rate(self, tr):
        return self.b.fill_Lh

    def take_load(self, tr, L, g, t):
        pl = self.pl
        pl.ds_g += g
        pl.ds_L += L
        pl.flow["L_ds"] += L
        pl.flow["ds_g"] += g
        tr.lot_g += g
        tr.lot_L += L
        tr.L_step += L

    def end_load(self, tr, t):
        pl = self.pl
        n = int(math.ceil(tr.lot_L / max(self.b.container_L, 1e-9) - 1e-9)) if tr.lot_L > 1e-9 else 0
        pl.containers += n
        pl.lots.append({"lot": tr.lot, "day": round(t / 24.0, 2), "g": tr.lot_g, "L": tr.lot_L, "containers": n})
        pl.event(pl.h, "lot", f"Lot {tr.lot}: {tr.lot_g:.0f} g in {n} cryovault containers")


class DSP:
    def __init__(self, plant, sump):
        p = plant.p
        self.pl, self.sump = plant, sump
        self.lowph = p.vi.mode == "low_ph"
        self.pav = Vessel("pav", "Protein A vessel SUM")
        self.vip = Vessel("vip", "Low-pH VI pool")
        self.dfs = Vessel("dfs", "Depth filtration surge tank")
        self.aexv = Vessel("aexv", "AEX vessel SUM")
        self.cexv = Vessel("cexv", "CEX vessel SUM")
        self.vrfv = Vessel("vrfv", "VRF vessel SUM")
        self.tank = Vessel("ufdf", "UF/DF vessel SUM")
        self.bdsm = Vessel("bdsm", "BDS hold mixer")
        self.pa = ProteinA(plant, sump, self.pav) if p.capture.mode != "pcc" else None
        cyc = capture_cycle(p)
        vol = (1 + p.vi.acid_bias) * (1 + p.vi.base_frac) if self.lowph else 1.0
        self.depth = DepthFilter(plant, self.vip if self.lowph else self.pav, self.dfs,
                                 max(p.depth.trigger_elutions, 1.0) * cyc["eluate_L"] * vol)
        if p.polish.order == "cex_aex":
            self.u1 = CEX(plant, self.dfs, self.cexv)
            self.u2 = AEX(plant, self.cexv, self.aexv)
            last = self.aexv
        else:
            self.u1 = AEX(plant, self.dfs, self.aexv)
            self.u2 = CEX(plant, self.aexv, self.cexv)
            last = self.cexv
        self.vrf = VirusFilter(plant, last, self.vrfv, self.u2)
        self.ufdf = UFDF(plant, self.vrfv, self.bdsm, self.tank)
        self.bds = BDS(plant, self.bdsm, self.ufdf)
        self.post = [self.depth, self.u1, self.u2, self.vrf, self.ufdf, self.bds]
        self.upstream_quiet = True         # set by the simulator: PCC column / low-pH tube empty
        self.vessels = ([self.pav] + ([self.vip] if self.lowph else [])
                        + [self.dfs, self.aexv, self.cexv, self.vrfv, self.tank, self.bdsm])

    def units(self):
        return ([self.pa] if self.pa else []) + self.post

    def _flags(self, tail):
        """End of campaign: once nothing more can arrive, process what is left, then store."""
        if self.pa is not None:
            self.pa.force = tail
            self.pa.finish = tail and not self.pa.busy() and self.sump.g <= EPS_G
            quiet = self.pa.finish
        else:
            quiet = tail
        quiet = quiet and self.upstream_quiet
        if self.lowph:
            quiet = quiet and self.pav.g <= EPS_G
        for u in self.post:
            u.force = quiet
            u.finish = quiet and not u.busy() and u.feed.g <= EPS_G
            quiet = u.finish

    def hour(self, h, tail):
        span = 1.0 / SLICES
        for k in range(SLICES):
            t0 = h + k * span
            self._flags(tail)
            for u in self.units():
                u.advance(t0, span)

    def levels(self):
        return {"pav_L": self.pav.L, "vip_L": self.vip.L, "dfs_L": self.dfs.L, "aexv_L": self.aexv.L,
                "cexv_L": self.cexv.L, "vrfv_L": self.vrfv.L, "ufdf_L": self.tank.L, "bdsm_L": self.bdsm.L}

    def lane_state(self, vi_on):
        out = dict.fromkeys(LANE_KEYS, "")
        for u in self.units():
            for tr in u.trains:
                s = tr.step
                if s is not None and "st_" + tr.lane in out:
                    out["st_" + tr.lane] = s["name"]
                    out["sk_" + tr.lane] = s["kind"]
        if self.pa is not None and self.pa.dose > 0 and any(tr.load is not None for tr in self.pa.trains):
            out["st_det"], out["sk_det"] = "Dosing Brij-9", "buf"
        for u in (self.u1, self.u2):
            if any(tr.load is not None for tr in u.trains):
                lane, name = ("base", "Tris base + NaCl") if isinstance(u, AEX) else ("acid", "Acetic acid")
                out["st_" + lane], out["sk_" + lane] = name, "buf"
        if vi_on:
            out["st_vi"], out["sk_vi"] = "Low-pH hold", "hold"
        return out

    def in_process_g(self):
        g = self.sump.g + sum(v.g for v in self.vessels)
        for u in self.units():
            for tr in u.trains:
                if tr.busy and any(s["kind"] in ("load", "elute") for s in tr.job[tr.i:]):
                    g += tr.bound
        return g

    def report(self):
        led = self.pl.ledger
        sols, chem = [], {}
        for k, L in led.sol.items():
            if L <= 1e-9:
                continue
            meta = SOLUTIONS.get(k, {"label": k, "where": "", "g_per_L": {}})
            sols.append({"key": k, "label": meta["label"], "where": meta["where"], "L": round(L, 2)})
            for c, gpl in meta["g_per_L"].items():
                chem[c] = chem.get(c, 0.0) + L * gpl / 1000.0
        sols.sort(key=lambda r: -r["L"])
        chems = [{"name": c, "kg": round(kg, 3)} for c, kg in sorted(chem.items(), key=lambda kv: -kv[1])
                 if kg > 1e-6]
        vessels = [{"key": v.key, "label": v.label, "peak_L": round(v.peak_L, 1), "end_L": round(v.L, 2),
                    "end_g": round(v.g, 1)} for v in [self.sump] + self.vessels]
        used = {s["lane"] for s in led.steps}
        lanes = [{"key": k, "label": lab} for k, lab in LANES if k in used]
        return {"solutions": sols, "chemicals": chems, "vessels": vessels, "lanes": lanes, "steps": led.steps}
