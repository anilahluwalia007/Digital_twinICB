"""
icb_nevai.py — NevAI LSTM forecasters for the IDB line
NavNevAI_mAbWasteWater · CarbonTesseract 2026

Reuses the LSTM core in nevai_engine.py. Two forecasters, hourly data:
  titre  permeate titre, next 24 h -> tomorrow's grams, cycles and buffer
  ww     total wastewater, next 24 h -> EQ tank / WWTP planning

Training starts on synthetic runs from icb_dynamic.py across all four
molecules (the model learns the simulator). Hourly plant data is added with
fine_tune(); rows need the keys in TITRE_KEYS and WW_KEYS.

Feature version 2 (IDB): older saved weights are ignored and retrained.

Train:  python3 icb_nevai.py train
"""
from __future__ import annotations
import json
import sys
import time
import numpy as np
import nevai_engine as nev

SEQ_LEN, HORIZON, HIDDEN = 48, 24, 32
FEATURE_VERSION = 2
ICB_DIR = nev.DATA_DIR / "icb"
ICB_DIR.mkdir(parents=True, exist_ok=True)

TITRE_KEYS = ["titre_meas", "VVD", "Xv_meas", "viability", "sieve", "S", "n_day", "phase"]
WW_KEYS = ["ww_total_Lh", "ww_cells_Lh", "ww_capture_flowthrough", "capture_switch",
           "ww_n1_permeate", "ww_sptff_permeate", "VVD", "n_day"]


def _g(r, k):
    v = r.get(k, 0.0)
    return float(v) if v is not None else 0.0


def _nd(r):
    return min(max(_g(r, "n_day") / 30.0, 0.0), 1.0)


def titre_features(hist, scale=1.0):
    return np.array([[_g(r, "titre_meas") / 5.0, _g(r, "VVD") / 3.0, _g(r, "Xv_meas") / 100.0,
                      _g(r, "viability"), _g(r, "sieve"), _g(r, "S") / 10.0, _nd(r),
                      _g(r, "phase")] for r in hist], dtype=np.float32)


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
    def __init__(self, name, target_key, feature_fn):
        self.name, self.target_key, self.feature_fn = name, target_key, feature_fn
        self.w, self.scale, self.meta = None, 1.0, {}

    @property
    def _wpath(self):
        return ICB_DIR / f"{self.name}_weights.npz"

    @property
    def _mpath(self):
        return ICB_DIR / f"{self.name}_meta.json"

    def _windows(self, hist, stride):
        F = self.feature_fn(hist, self.scale)
        y = np.array([_g(r, self.target_key) for r in hist], dtype=np.float32) / self.scale
        X, Y = [], []
        for i in range(0, len(hist) - SEQ_LEN - HORIZON + 1, max(1, stride)):
            X.append(F[i:i + SEQ_LEN])
            Y.append(y[i + SEQ_LEN:i + SEQ_LEN + HORIZON])
        return X, Y

    def fit(self, histories, epochs=8, stride=12):
        self.scale = max(max(_g(r, self.target_key) for r in h) for h in histories if h) or 1.0
        X, Y = [], []
        for h in histories:
            a, b = self._windows(h, stride)
            X += a
            Y += b
        if not X:
            raise ValueError(f"Runs too short: need more than {SEQ_LEN + HORIZON} hourly rows")
        X, Y = np.array(X, dtype=np.float32), np.array(Y, dtype=np.float32)
        self.w = nev.NevAIWeights(X.shape[2], HIDDEN, HORIZON)
        losses = nev.train(self.w, X, Y, epochs=epochs)
        self.meta = {"trained": True, "feature_version": FEATURE_VERSION, "windows": int(len(X)),
                     "epochs": int(epochs),
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
        losses = nev.train(self.w, np.array(X, dtype=np.float32), np.array(Y, dtype=np.float32),
                           epochs=epochs, lr=nev.LR * 0.3)
        self.meta["plant_windows"] = self.meta.get("plant_windows", 0) + len(X)
        self.save()
        return {"fine_tuned": True, "windows": len(X), "loss": round(float(losses[-1]), 5)}

    def forecast(self, hist):
        if self.w is None or len(hist) < SEQ_LEN:
            return None
        preds, _, _ = nev.nevai_forward(self.w, self.feature_fn(hist[-SEQ_LEN:], self.scale))
        return np.clip(preds, 0, None) * self.scale

    def forecast_band(self, hist, n_samples=30):
        if self.w is None or len(hist) < SEQ_LEN:
            return None
        out = nev.predict(self.w, self.feature_fn(hist[-SEQ_LEN:], self.scale), self.scale, n_samples)
        return {"hours_ahead": out["days_ahead"], "mean": out["forecast_gpd"],
                "p10": out["ci_lo_gpd"], "p90": out["ci_hi_gpd"],
                "note": "band from weight perturbation, not a calibrated interval"}

    def save(self):
        np.savez_compressed(self._wpath, Wx=self.w.Wx, Wh=self.w.Wh, b=self.w.b,
                            Wo=self.w.Wo, bo=self.w.bo,
                            dims=np.array([self.w.input_dim, self.w.hidden_dim, self.w.horizon]))
        self._mpath.write_text(json.dumps({**self.meta, "scale": float(self.scale)}, indent=2))

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
        self.titre = Forecaster("titre", "titre_meas", titre_features)
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


def synthetic_histories(base_params=None, n_runs=4, seed=0, jitter=0.15):
    """Varied simulator runs. IDB mode cycles through the molecules and varies
    curve height, curve speed and feed rate; the legacy model varies growth,
    productivity, perfusion and fouling. Each run gets one random titre step."""
    from icb_process import ICBParams, TITRE_PROFILES
    from icb_dynamic import ICBSimulator
    rng = np.random.default_rng(seed)
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
            p.upstream.molecule = names[k % len(names)]
            p.upstream.titre_scale *= float(rng.uniform(1 - jitter, 1 + jitter))
            p.upstream.time_stretch *= float(rng.uniform(0.9, 1.1))
            p.upstream.VVD *= float(rng.uniform(1 - jitter, 1 + jitter))
        p.disturbances = list(p.disturbances) + [{
            "type": "qp_step", "day": float(rng.uniform(3, max(4, p.run_days - 2))),
            "factor": float(rng.uniform(0.8, 1.1))}]
        hists.append(ICBSimulator(p, forecaster=None, seed=seed + k).run(keep_history=True)["_history"])
    return hists


def train_forecasters(base_params=None, n_runs=4, epochs=8, stride=12, seed=0):
    t0 = time.time()
    hists = synthetic_histories(base_params, n_runs, seed)
    fs = ICBForecasters()
    out = {"titre": fs.titre.fit(hists, epochs, stride), "ww": fs.ww.fit(hists, epochs, stride)}
    out.update({"runs": n_runs, "seconds": round(time.time() - t0, 1)})
    return out


def status():
    fs = ICBForecasters()
    fs.load()
    return {"ready": fs.ready, "titre": fs.titre.meta, "ww": fs.ww.meta}


if __name__ == "__main__":
    from icb_dynamic import ICBSimulator
    if len(sys.argv) > 1 and sys.argv[1] == "train":
        # Full local training: 8 runs (each molecule twice) x 30 epochs.
        # The web app's background training uses a lighter default (4 x 8).
        print("Training NevAI forecasters (a few minutes)...")
        print(json.dumps(train_forecasters(n_runs=8, epochs=30), indent=2))
    fs = ICBForecasters()
    if fs.load():
        res = ICBSimulator(forecaster=fs, seed=7).run()
        print(json.dumps(res["kpi"], indent=2))
    else:
        print("No trained forecasters yet. Run:  python3 icb_nevai.py train")
