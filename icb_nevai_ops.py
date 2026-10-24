"""
icb_nevai_ops.py — NevAI in operation: act on forecasts, watch for drift, plan buffer preparation
NavNevAI_mAbWasteWater · CarbonTesseract 2026

Used by icb_dynamic.ICBSimulator when a forecaster is loaded.

  HoldGuard     early partial Protein A load when the forecast says the surge-tank product would
                outlive its hold limit before a full load is ready (only with control.nevai_act)
  Watch         compares each measured permeate titre with the forecast issued nevai_watch_lead_h
                earlier and flags a shift in the forecast error (sensor fault, productivity drop,
                fouling) before the lab sees it; reports only. The error is a share of the forecast,
                but never of less than nevai_watch_floor_gL (default 1 g/L), so wobble at start-up
                titres of a few tenths of a g/L does not count as a shift
  prep_compare  buffer preparation driven by the NevAI 24 h forecast against a reactive rule
                (prepare what was used in the last 24 h): same prep capacity, same hourly use
  report        everything above for result["nevai"]

Nothing here changes the plant unless control.nevai_act is True.
"""
from __future__ import annotations
import numpy as np

WATCH_FLOOR_GL = 1.0      # watch error = (measured - forecast) / max(forecast, this); override with
                          # control.nevai_watch_floor_gL if that parameter is added


class HoldGuard:
    """Decides when NevAI should start a partial Protein A load to protect the hold limit."""

    def __init__(self, ct, cp, target_g):
        self.on = bool(ct.nevai_act)
        self.margin = max(float(ct.nevai_hold_margin_h), 0.0)
        self.min_g = max(float(ct.nevai_partial_min_frac), 0.0) * target_g
        self.limit = float(cp.max_sump_hold_h)

    def check(self, switch_fc_h, sump_L, sump_g_meas, Fh, loading):
        """True when a partial load should start now. sump_g_meas = grams in the tank by PAT."""
        if not self.on or loading or Fh <= 0 or sump_L <= 0:
            return False
        age = sump_L / Fh                                         # mean product age now, h
        wait = switch_fc_h if switch_fc_h is not None else 24.0   # beyond the horizon: 24 h or more
        return age + wait > self.limit - self.margin and sump_g_meas >= self.min_g


class Watch:
    """Flags a shift of measured titre away from the NevAI forecast trend (reports only)."""

    def __init__(self, ct):
        self.lead = max(1, int(round(ct.nevai_watch_lead_h)))
        self.win = max(2, int(round(ct.nevai_watch_window_h)))
        self.base = max(6, int(round(ct.nevai_watch_baseline_h)))
        self.base_min = min(24, self.base)          # baseline points needed before scoring
        self.thr = float(ct.nevai_watch_threshold)
        self.floor = max(float(getattr(ct, "nevai_watch_floor_gL", WATCH_FLOOR_GL)), 1e-6)
        self.issued, self.err = {}, {}
        self.flag_hours, self.first_flag_h, self.peak = 0, None, 0.0
        self.starts, self._last = [], -10

    def issue(self, h, tf):
        """Keep the 24 h titre forecast issued at hour h (tf[0] is hour h + 1)."""
        if tf is not None:
            self.issued[h] = np.asarray(tf, dtype=float)
        for k in [k for k in self.issued if k < h - 48]:
            del self.issued[k]

    def score(self, h, actual):
        """Shift of the recent forecast error from its own baseline (share of the forecast)."""
        if actual is None:
            return None
        ok = [i for i, a in self.issued.items() if i <= h - self.lead and 1 <= h - i <= len(a)]
        if not ok:
            return None
        i = max(ok)
        pred = float(self.issued[i][h - i - 1])
        if pred < 0.05:
            return None
        self.err[h] = (actual - pred) / max(pred, self.floor)
        oldest = h - self.win - self.base
        for k in [k for k in self.err if k <= oldest]:
            del self.err[k]
        recent = [self.err[k] for k in range(h - self.win + 1, h + 1) if k in self.err]
        base = [self.err[k] for k in range(oldest + 1, h - self.win + 1) if k in self.err]
        if len(recent) < max(2, self.win // 2) or len(base) < self.base_min:
            return None
        s = abs(float(np.mean(recent)) - float(np.median(base)))   # a constant model bias cancels
        self.peak = max(self.peak, s)
        return s

    def flag(self, h, s):
        if s is None or s <= self.thr:
            return False
        self.flag_hours += 1
        if self.first_flag_h is None:
            self.first_flag_h = h
        if h - self._last > 1:
            self.starts.append(h)
        self._last = h
        return True


def prep_compare(history, ct):
    """Hold-tank simulation of buffer preparation: NevAI forecast vs reactive rule."""
    use = [float(r.get("buffer_L_h") or 0.0) for r in history]
    fc = [r.get("buffer_fc_24h_L") for r in history]
    if not any(use) or not any(v is not None for v in fc):
        return None
    daily = [sum(use[i:i + 24]) for i in range(0, len(use), 24)]
    cap = float(ct.buffer_prep_Lh) if ct.buffer_prep_Lh > 0 else max(daily) / 20.0
    days = max(float(ct.buffer_hold_days), 0.1)
    out = {"prep_capacity_Lh": round(cap, 2)}
    for name in ("reactive", "nevai"):
        level = short = hold_peak = prep_peak = 0.0
        short_h, last_fc, last_h = 0, None, -1e9
        for h, u in enumerate(use):
            if fc[h] is not None:
                last_fc, last_h = fc[h], h
            past = sum(use[max(0, h - 24):h])
            fresh = name == "nevai" and last_fc is not None and h - last_h <= 6
            target = (last_fc if fresh else past) * days
            prep = min(max(target - level, 0.0), cap)
            level += prep
            prep_peak, hold_peak = max(prep_peak, prep), max(hold_peak, level)
            level -= u
            if level < 0:
                short, short_h, level = short - level, short_h + 1, 0.0
        out[name] = {"shortfall_L": round(short, 1), "short_hours": short_h,
                     "hold_peak_L": round(hold_peak, 1), "prep_peak_Lh": round(prep_peak, 2)}
    return out


def report(p, used, partial_loads, watch, H1, history):
    if not used:
        return None
    det = []
    for d in p.disturbances or []:
        sh = H1 + float(d.get("day", 1e9)) * 24.0
        after = [f for f in (watch.starts if watch else []) if f >= sh]
        det.append({"type": d.get("type"), "n_day": d.get("day"),
                    "detected_after_h": round(after[0] - sh, 1) if after else None})
    w = watch
    return {"acting": bool(p.control.nevai_act), "partial_loads": partial_loads,
            "watch": {"flag_hours": w.flag_hours if w else 0,
                      "first_flag_n_day": round((w.first_flag_h - H1) / 24.0, 2)
                      if w and w.first_flag_h is not None else None,
                      "peak_score": round(w.peak, 3) if w else None,
                      "threshold": w.thr if w else None},
            "detection": det, "prep": prep_compare(history, p.control)}
