"""
icb_nevai.py — NevAI LSTM forecasters for the IDB line
NavNevAI_mAbWasteWater · CarbonTesseract 2026

Reuses the LSTM core in nevai_engine.py. Two forecasters, hourly data:
  titre  permeate titre, next 24 h -> tomorrow's grams, cycles and buffer
  ww     total wastewater, next 24 h -> EQ tank / WWTP planning

Training starts on synthetic runs from icb_dynamic.py across all four
molecules (the model learns the simulator). Hourly plant data is added with
fine_tune(); rows need the keys in TITRE_KEYS and WW_KEYS.

Forecast bands are split-conformal 10-90 % intervals from held-out simulator runs
(train_forecasters(n_val=...)). Fine-tuning keeps the old band: re-run calibrate() on
plant histories once a few days of plant data are in.

Feature version 4: the titre forecaster is relative. It sees the last 48 h of titre divided by
the current titre (median of the last REF_H measurements) and predicts the next 24 h as a ratio
of that level, which is multiplied back. It therefore has no ceiling: a case that runs at
twice the titre of anything in training is forecast as well as one inside it (version 3
predicted absolute g/L and could not go above the highest titre it was trained on).
Windows whose 24 h target runs past the end of the N run are no longer used for training
(the run-out rows have no titre and taught the model a drop to zero). Older saved weights
are ignored and retrained.

What NevAI is trained on
  Simulator runs only (no plant data unless you fine-tune): each run is icb_dynamic.py with the
  parameters sent from the Inputs tab (the selected case is not applied), the titre curve cycled
  through the four Figure 2 molecules (plus a custom curve if one is set), its height multiplied
  by a random factor in TRAIN_DEFAULTS scale range (log-uniform), its speed by 0.9-1.1, the
  perfusion rate by 0.85-1.15, and one random titre step (x0.8-1.1) somewhere in the run.
  "runs" train the networks; "val_runs" (different molecules, other seeds) set the 10-90 % band.
  Both networks: LSTM, 48 h look-back, 24 h ahead, 32 hidden units.

Training settings (TRAIN_DEFAULTS, clamped to TRAIN_LIMITS)
  From the web page: send them with the params as params["nevai_train"] =
  {"runs", "val_runs", "epochs", "scale_lo", "scale_hi", "seed"}; ICBParams ignores the key.
  Progress is written to data/icb/training_progress.json and returned by status()["progress"].

Train it yourself (command line, then copy data/icb/*_weights.npz and *_meta.json to the server's
DATA_DIR/icb, or run it on the server):
  python3 icb_nevai.py train                                   8 runs x 30 epochs (thorough)
  python3 icb_nevai.py train --runs 12 --epochs 40 --scale-lo 0.5 --scale-hi 3 --params inputs.json
  python3 icb_nevai.py finetune plant.csv --epochs 5           adapt to hourly plant data
  python3 icb_nevai.py status
"""
from __future__ import annotations
import json
import math
import os
import sys
import time
import numpy as np
import nevai_engine as nev

SEQ_LEN, HORIZON, HIDDEN = 48, 24, 32
FEATURE_VERSION = 4       # relative titre forecast (no ceiling); run-out windows left out
ICB_DIR = nev.DATA_DIR / "icb"
ICB_DIR.mkdir(parents=True, exist_ok=True)

REF_H = 6                 # current titre = median of the last REF_H measurements in the window
MIN_OBS = 12              # measured titre hours a window needs before NevAI forecasts from it
SCALE_RANGE = (0.5, 2.0)  # training multiplies each run's titre curve by a factor in this range

TRAIN_DEFAULTS = {"runs": 4, "val_runs": 2, "epochs": 8, "scale_lo": SCALE_RANGE[0],
                  "scale_hi": SCALE_RANGE[1], "seed": 0}
TRAIN_LIMITS = {"runs": (2, 24), "val_runs": (0, 8), "epochs": (1, 80), "scale_lo": (0.1, 1.0),
                "scale_hi": (1.0, 5.0), "seed": (0, 10 ** 6)}
PROGRESS_PATH = ICB_DIR / "training_progress.json"


def train_settings(overrides=None):
    """Defaults, then any overrides, each clamped to TRAIN_LIMITS (bad values keep the default)."""
    out = dict(TRAIN_DEFAULTS)
    for k, v in (overrides or {}).items():
        if k not in out or v is None or v == "":
            continue
        try:
            lo, hi = TRAIN_LIMITS[k]
            v = float(v)
            v = min(max(v, lo), hi)
            out[k] = int(round(v)) if k in ("runs", "val_runs", "epochs", "seed") else v
        except (TypeError, ValueError):
            pass
    if out["scale_hi"] < out["scale_lo"]:
        out["scale_lo"], out["scale_hi"] = out["scale_hi"], out["scale_lo"]
    return out


def _progress(state, stage, done, total, t0, **extra):
    """Write the training progress file (read back by status())."""
    rec = {"state": state, "stage": stage, "done": int(done), "total": int(total),
           "pct": round(100.0 * done / total, 1) if total else None,
           "started": round(t0, 1), "updated": round(time.time(), 1),
           "elapsed_s": round(time.time() - t0, 1), **extra}
    try:
        tmp = PROGRESS_PATH.with_name(PROGRESS_PATH.name + ".tmp")
        tmp.write_text(json.dumps(rec))
        os.replace(tmp, PROGRESS_PATH)
    except OSError:
        pass
    return rec

TITRE_KEYS = ["titre_meas", "VVD", "Xv_meas", "viability", "sieve", "S", "n_day", "phase"]
WW_KEYS = ["ww_total_Lh", "ww_cells_Lh", "ww_capture_flowthrough", "capture_switch",
           "ww_n1_permeate", "ww_sptff_permeate", "VVD", "n_day"]


def _g(r, k):
    v = r.get(k, 0.0)
    return float(v) if v is not None else 0.0


def _nd(r):
    return min(max(_g(r, "n_day") / 30.0, 0.0), 1.0)


def titre_level(win, key="titre_meas"):
    """Current titre of a window: median of its last REF_H measurements (None if too few)."""
    vals = [float(r[key]) for r in win if r.get(key) is not None]
    if len(vals) < MIN_OBS:
        return None
    return max(float(np.median(vals[-REF_H:])), 0.02)


def titre_features(hist, ref=1.0):
    """Titre as a ratio of the current level, plus a flag for hours with a measurement."""
    return np.array([[_g(r, "titre_meas") / ref, 1.0 if r.get("titre_meas") is not None else 0.0,
                      _g(r, "VVD") / 3.0, _g(r, "Xv_meas") / 100.0, _g(r, "viability"), _g(r, "sieve"),
                      _g(r, "S") / 10.0, _nd(r), _g(r, "phase")] for r in hist], dtype=np.float32)


def ww_features(hist, scale=1.0):
    rows = []
    for r in hist:
        tot = _g(r, "ww_total_Lh")
        rows.append([tot / scale, _g(r, "ww_cells_Lh") / (tot + 1e-6),
                     _g(r, "ww_capture_flowthrough") / (tot + 1e-6), _g(r, "capture_switch") / 3.0,
                     _g(r, "ww_n1_permeate") / (tot + 1e-6), _g(r, "ww_sptff_permeate") / (tot + 1e-6),
                     _g(r, "VVD") / 3.0, _nd(r)])
    return np.array(rows, dtype=np.float32)


class Forecaster:
    """relative=True: inputs and targets are divided by the window's current level (titre_level),
    so the forecast follows any titre height. relative=False: one fixed scale (wastewater)."""

    def __init__(self, name, target_key, feature_fn, relative=False):
        self.name, self.target_key, self.feature_fn = name, target_key, feature_fn
        self.relative = relative
        self.offset = 1.0 if relative else 0.0   # relative: learn the change from now (0 = persistence)
        self.w, self.scale, self.meta = None, 1.0, {}

    @property
    def _wpath(self):
        return ICB_DIR / f"{self.name}_weights.npz"

    @property
    def _mpath(self):
        return ICB_DIR / f"{self.name}_meta.json"

    def _norm(self, win):
        """Divisor for one window: its current level (relative) or the trained scale."""
        return titre_level(win, self.target_key) if self.relative else self.scale

    def _windows(self, hist, stride):
        X, Y = [], []
        for i in range(0, len(hist) - SEQ_LEN - HORIZON + 1, max(1, stride)):
            win, tgt = hist[i:i + SEQ_LEN], hist[i + SEQ_LEN:i + SEQ_LEN + HORIZON]
            if self.relative and any(r.get(self.target_key) is None for r in tgt):
                continue                       # target runs past the N run: no titre to learn from
            norm = self._norm(win)
            if not norm:
                continue
            X.append(self.feature_fn(win, norm))
            Y.append(np.array([_g(r, self.target_key) for r in tgt], dtype=np.float32) / norm)
        return X, Y

    def _level(self, pred):
        """Network output -> forecast in units of the window divisor (never below zero)."""
        return np.clip(np.asarray(pred, dtype=float) + self.offset, 0, None)

    def _stride(self, stride):
        """The relative model uses every window a third as far apart: more examples, same runs."""
        return max(1, int(stride) // 3) if self.relative else stride

    def fit(self, histories, epochs=8, stride=12, on_epoch=None):
        self.scale = 1.0 if self.relative else (
            max(max(_g(r, self.target_key) for r in h) for h in histories if h) or 1.0)
        X, Y = [], []
        for h in histories:
            a, b = self._windows(h, self._stride(stride))
            X += a
            Y += b
        if not X:
            raise ValueError(f"Runs too short: need more than {SEQ_LEN + HORIZON} hourly rows")
        X, Y = np.array(X, dtype=np.float32), np.array(Y, dtype=np.float32) - self.offset
        self.w = nev.NevAIWeights(X.shape[2], HIDDEN, HORIZON)
        losses = nev.train_batched(self.w, X, Y, epochs=epochs, on_epoch=on_epoch)
        self.meta = {"trained": True, "feature_version": FEATURE_VERSION, "windows": int(len(X)),
                     "epochs": int(epochs), "relative": self.relative, "offset": self.offset,
                     "final_loss": round(float(losses[-1]), 5) if losses else None,
                     "seq_len_h": SEQ_LEN, "horizon_h": HORIZON, "hidden": HIDDEN,
                     "last_trained": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
        self.save()
        return self.meta

    def fine_tune(self, hist, epochs=5, stride=1):
        if self.w is None:
            raise RuntimeError("Train or load the model first")
        X, Y = self._windows(hist, stride)
        if not X:
            return {"fine_tuned": False, "reason": f"need more than {SEQ_LEN + HORIZON} hourly rows"}
        losses = nev.train_batched(self.w, np.array(X, dtype=np.float32),
                                   np.array(Y, dtype=np.float32) - self.offset, epochs=epochs, lr=nev.LR * 0.3)
        self.meta["plant_windows"] = self.meta.get("plant_windows", 0) + len(X)
        self.save()
        return {"fine_tuned": True, "windows": len(X), "loss": round(float(losses[-1]), 5)}

    def forecast(self, hist):
        """Next HORIZON hours, or None (not trained, too little history, or too few measurements)."""
        if self.w is None or len(hist) < SEQ_LEN:
            return None
        win = hist[-SEQ_LEN:]
        norm = self._norm(win)
        if not norm:
            return None
        preds, _, _ = nev.nevai_forward(self.w, self.feature_fn(win, norm))
        return self._level(preds) * norm

    def calibrate(self, histories, stride=6):
        """Split-conformal 10-90 % band from runs the model did not train on
        (in units of the window divisor: a share of the current titre for the relative model)."""
        if self.w is None:
            return None
        res = []
        for hist in histories:
            X, Y = self._windows(hist, stride)
            for x, y in zip(X, Y):
                p, _, _ = nev.nevai_forward(self.w, x)
                res.append(np.asarray(y, dtype=float) - self._level(p))
        if len(res) < 10:
            return None
        R = np.array(res)
        self.meta["band_q10"] = np.percentile(R, 10, axis=0).round(5).tolist()
        self.meta["band_q90"] = np.percentile(R, 90, axis=0).round(5).tolist()
        self.meta["band_windows"] = int(len(R))
        self.save()
        return {"windows": int(len(R)), "band": "10-90 %, split-conformal"}

    def forecast_band(self, hist, n_samples=30):
        if self.w is None or len(hist) < SEQ_LEN:
            return None
        win = hist[-SEQ_LEN:]
        norm = self._norm(win)
        if not norm:
            return None
        q10, q90 = self.meta.get("band_q10"), self.meta.get("band_q90")
        if q10 and q90 and len(q10) == HORIZON and len(q90) == HORIZON:
            mean = self.forecast(hist)
            return {"hours_ahead": list(range(1, HORIZON + 1)), "mean": mean.tolist(),
                    "p10": np.clip(mean + np.array(q10) * norm, 0, None).tolist(),
                    "p90": np.clip(mean + np.array(q90) * norm, 0, None).tolist(),
                    "note": "split-conformal band from held-out simulator runs"}
        x = self.feature_fn(win, norm)
        samples = []
        for _ in range(n_samples):                         # same weight perturbation as nev.predict
            w2 = nev.NevAIWeights(self.w.input_dim, self.w.hidden_dim, self.w.horizon)
            w2.Wx = self.w.Wx + nev.rng.standard_normal(self.w.Wx.shape) * 0.02 * np.std(self.w.Wx)
            w2.Wh = self.w.Wh + nev.rng.standard_normal(self.w.Wh.shape) * 0.02 * np.std(self.w.Wh)
            w2.b, w2.Wo, w2.bo = self.w.b.copy(), self.w.Wo.copy(), self.w.bo.copy()
            p, _, _ = nev.nevai_forward(w2, x)
            samples.append(self._level(p) * norm)
        S = np.array(samples)
        return {"hours_ahead": list(range(1, HORIZON + 1)), "mean": S.mean(axis=0).tolist(),
                "p10": np.percentile(S, 10, axis=0).tolist(), "p90": np.percentile(S, 90, axis=0).tolist(),
                "note": "band from weight perturbation, not a calibrated interval (retrain to calibrate)"}

    def save(self):
        """Write to temporary files, then rename: readers never see half-written weights."""
        tmp_w = self._wpath.with_name(self._wpath.stem + ".tmp.npz")
        np.savez_compressed(tmp_w, Wx=self.w.Wx, Wh=self.w.Wh, b=self.w.b,
                            Wo=self.w.Wo, bo=self.w.bo,
                            dims=np.array([self.w.input_dim, self.w.hidden_dim, self.w.horizon]))
        os.replace(tmp_w, self._wpath)
        tmp_m = self._mpath.with_name(self._mpath.name + ".tmp")
        tmp_m.write_text(json.dumps({**self.meta, "scale": float(self.scale)}, indent=2))
        os.replace(tmp_m, self._mpath)

    def load(self):
        if not (self._wpath.exists() and self._mpath.exists()):
            return False
        try:
            meta = json.loads(self._mpath.read_text())
            if meta.get("feature_version") != FEATURE_VERSION:
                return False
            with np.load(self._wpath) as d:
                dims = d["dims"]
                w = nev.NevAIWeights(int(dims[0]), int(dims[1]), int(dims[2]))
                w.Wx, w.Wh, w.b = np.array(d["Wx"]), np.array(d["Wh"]), np.array(d["b"])
                w.Wo, w.bo = np.array(d["Wo"]), np.array(d["bo"])
            self.w = w
            self.meta = meta
            self.scale = float(meta.get("scale", 1.0))
            return True
        except Exception:
            return False


class ICBForecasters:
    seq_len, horizon = SEQ_LEN, HORIZON

    def __init__(self):
        self.titre = Forecaster("titre", "titre_meas", titre_features, relative=True)
        self.ww = Forecaster("ww", "ww_total_Lh", ww_features)

    @property
    def ready(self):
        return self.titre.w is not None and self.ww.w is not None

    def load(self):
        a = self.titre.load()
        b = self.ww.load()
        return a and b

    def forecast_all(self, hist):
        return {"titre": self.titre.forecast(hist), "ww": self.ww.forecast(hist)}


def synthetic_histories(base_params=None, n_runs=4, seed=0, jitter=0.15, offset=0,
                        scale_range=SCALE_RANGE, on_run=None):
    """Varied simulator runs. IDB mode cycles through the molecules and varies
    curve height (scale_range, default ×0.5 to ×2, log-uniform, so cases far above the Inputs
    titre are covered),
    curve speed and feed rate; the legacy model varies growth, productivity, perfusion and
    fouling. Each run gets one random titre step.
    offset shifts the molecule cycle, so validation runs use different molecules."""
    from icb_process import ICBParams, TITRE_PROFILES
    from icb_dynamic import ICBSimulator
    rng = np.random.default_rng(seed)
    lo, hi = math.log(scale_range[0]), math.log(scale_range[1])
    hists = []
    for k in range(n_runs):
        p = ICBParams.from_dict(base_params)
        if p.upstream.mode == "mechanistic":
            b = p.bioreactor
            b.qP_pg_cell_h *= float(rng.uniform(1 - jitter, 1 + jitter))
            b.mu_max_h *= float(rng.uniform(1 - jitter, 1 + jitter))
            b.sieve_loss_per_day *= float(rng.uniform(0.5, 1.5))
            p.control.CSPR_pL_cell_day *= float(rng.uniform(1 - jitter, 1 + jitter))
        else:
            names = list(TITRE_PROFILES) + (["custom"] if p.upstream.custom_points else [])
            p.upstream.molecule = names[(k + offset) % len(names)]
            p.upstream.titre_scale *= float(math.exp(rng.uniform(lo, hi)))
            p.upstream.time_stretch *= float(rng.uniform(0.9, 1.1))
            p.upstream.VVD *= float(rng.uniform(1 - jitter, 1 + jitter))
        p.disturbances = list(p.disturbances) + [{
            "type": "qp_step", "day": float(rng.uniform(3, max(4, p.run_days - 2))),
            "factor": float(rng.uniform(0.8, 1.1))}]
        hists.append(ICBSimulator(p, forecaster=None, seed=seed + k).run(keep_history=True)["_history"])
        if on_run:
            on_run(k + 1)
    return hists


def train_forecasters(base_params=None, n_runs=None, epochs=None, stride=12, seed=None, n_val=None,
                      settings=None):
    """Train both forecasters, then calibrate their 10-90 % bands on held-out runs.
    Settings: TRAIN_DEFAULTS, then the arguments given here, then base_params["nevai_train"]
    (what the web page sends). Progress goes to PROGRESS_PATH at every run and epoch."""
    t0 = time.time()
    bp = dict(base_params or {})
    given = {"runs": n_runs, "epochs": epochs, "seed": seed, "val_runs": n_val, **(settings or {})}
    cfg = train_settings({**{k: v for k, v in given.items() if v is not None}, **(bp.pop("nevai_train", None) or {})})
    R, V, E = cfg["runs"], cfg["val_runs"], cfg["epochs"]
    total = R + V + 2 * E
    rng_scale = (cfg["scale_lo"], cfg["scale_hi"])
    done = [0]

    def step(stage, **extra):
        done[0] += 1
        _progress("running", stage, done[0], total, t0, settings=cfg, **extra)

    _progress("running", f"Simulating training run 1 of {R}", 0, total, t0, settings=cfg)
    try:
        hists = synthetic_histories(bp, R, cfg["seed"], scale_range=rng_scale,
                                    on_run=lambda k: step(f"Simulated training run {k} of {R}"))
        fs = ICBForecasters()
        out = {"titre": fs.titre.fit(hists, E, stride, on_epoch=lambda e, L: step(
                   f"Titre network: epoch {e} of {E}", loss=round(L, 5))),
               "ww": fs.ww.fit(hists, E, stride, on_epoch=lambda e, L: step(
                   f"Wastewater network: epoch {e} of {E}", loss=round(L, 5)))}
        if V > 0:
            val = synthetic_histories(bp, V, cfg["seed"] + 1000, offset=R, scale_range=rng_scale,
                                      on_run=lambda k: step(f"Simulated validation run {k} of {V}"))
            out["titre_band"] = fs.titre.calibrate(val)
            out["ww_band"] = fs.ww.calibrate(val)
        up = (bp.get("upstream") or {})
        trained_on = {"source": "simulator runs (Inputs tab parameters, no case applied)",
                      "runs": R, "val_runs": V, "epochs": E, "seed": cfg["seed"],
                      "titre_scale_range": list(rng_scale), "upstream_mode": up.get("mode", "profile"),
                      "run_days": bp.get("run_days"), "seconds": round(time.time() - t0, 1)}
        for f in (fs.titre, fs.ww):
            f.meta["trained_on"] = trained_on
            f.save()
        out.update({"runs": R, "settings": cfg, "seconds": round(time.time() - t0, 1)})
        _progress("done", "Finished", total, total, t0, settings=cfg)
        return out
    except Exception as e:
        _progress("failed", f"Failed: {e}", done[0], total, t0, settings=cfg, error=str(e))
        raise


def read_progress():
    try:
        return json.loads(PROGRESS_PATH.read_text())
    except (OSError, ValueError):
        return None


def status():
    fs = ICBForecasters()
    try:
        fs.load()
    except Exception:
        pass
    return {"ready": fs.ready, "titre": fs.titre.meta, "ww": fs.ww.meta, "progress": read_progress(),
            "train_defaults": TRAIN_DEFAULTS, "train_limits": TRAIN_LIMITS,
            "server_time": round(time.time(), 1)}


def weights_report():
    """Which weight files the server actually has (shown by /api/icb/nevai/status)."""
    out = {"dir": str(ICB_DIR.resolve()), "expected_version": FEATURE_VERSION}
    for name in ("titre", "ww"):
        w, m = ICB_DIR / f"{name}_weights.npz", ICB_DIR / f"{name}_meta.json"
        ver = None
        if m.exists():
            try:
                ver = json.loads(m.read_text()).get("feature_version")
            except Exception:
                ver = "unreadable"
        out[name] = {"weights_file": w.exists(), "meta_file": m.exists(), "feature_version": ver}
    return out


def read_plant_csv(path):
    """Hourly plant data: a header row with any of TITRE_KEYS / WW_KEYS; blank cells = no reading."""
    import csv
    rows = []
    with open(path, newline="", encoding="utf-8-sig") as f:
        for r in csv.DictReader(f):
            row = {}
            for k, v in r.items():
                if k is None:
                    continue
                v = (v or "").strip()
                try:
                    row[k.strip()] = float(v) if v else None
                except ValueError:
                    row[k.strip()] = v
            rows.append(row)
    return rows


def _cli():
    import argparse
    ap = argparse.ArgumentParser(description="Train, fine-tune or inspect the NevAI forecasters.")
    sub = ap.add_subparsers(dest="cmd")
    t = sub.add_parser("train", help="train on simulator runs")
    t.add_argument("--runs", type=int, default=8)
    t.add_argument("--val-runs", type=int, default=2)
    t.add_argument("--epochs", type=int, default=30)
    t.add_argument("--scale-lo", type=float, default=SCALE_RANGE[0])
    t.add_argument("--scale-hi", type=float, default=SCALE_RANGE[1])
    t.add_argument("--seed", type=int, default=0)
    t.add_argument("--params", help="JSON file of parameters (as the Inputs tab sends them)")
    f = sub.add_parser("finetune", help="adapt the trained networks to hourly plant data (CSV)")
    f.add_argument("csv")
    f.add_argument("--epochs", type=int, default=5)
    sub.add_parser("status", help="show what the saved networks were trained on")
    a = ap.parse_args()
    if a.cmd == "train":
        base = json.loads(open(a.params, encoding="utf-8").read()) if a.params else None
        print("Training NevAI forecasters (a few minutes)...")
        out = train_forecasters(base, settings={"runs": a.runs, "val_runs": a.val_runs, "epochs": a.epochs,
                                                "scale_lo": a.scale_lo, "scale_hi": a.scale_hi, "seed": a.seed})
        print(json.dumps(out, indent=2))
        print(f"Weights written to {ICB_DIR.resolve()}")
    elif a.cmd == "finetune":
        fs = ICBForecasters()
        if not fs.load():
            sys.exit("No trained forecasters yet. Run:  python3 icb_nevai.py train")
        rows = read_plant_csv(a.csv)
        print(f"{len(rows)} hourly rows from {a.csv}")
        for name, f_, keys in (("titre", fs.titre, TITRE_KEYS), ("ww", fs.ww, WW_KEYS)):
            if not rows or f_.target_key not in rows[0]:
                print(f"  {name}: skipped (no '{f_.target_key}' column)")
                continue
            missing = [k for k in keys if k not in rows[0]]
            print(f"  {name}: {f_.fine_tune(rows, epochs=a.epochs)}"
                  + (f"  (columns missing, treated as 0: {', '.join(missing)})" if missing else ""))
        print("The 10-90 % band is unchanged; recalibrate once more plant data is in.")
    elif a.cmd == "status":
        print(json.dumps({**status(), "files": weights_report()}, indent=2))
    else:
        from icb_dynamic import ICBSimulator
        fs = ICBForecasters()
        if fs.load():
            res = ICBSimulator(forecaster=fs, seed=7).run()
            print(json.dumps(res["kpi"], indent=2))
        else:
            print("No trained forecasters yet. Run:  python3 icb_nevai.py train")


if __name__ == "__main__":
    _cli()
