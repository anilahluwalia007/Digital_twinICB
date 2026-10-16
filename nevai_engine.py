"""
nevai_engine.py  —  Digital_Twin_mAbWasteWater  © CarbonTesseract 2026
─────────────────────────────────────────────────────────────────────
Pure-NumPy NevAI (LSTM) for wastewater flow prediction.
No TensorFlow / PyTorch dependency; runs on Render free tier.

Compatible with Python 3.9+ (macOS system Python). Same API as the original:
NevAIWeights, nevai_forward, nevai_backward, train, predict, anomaly_score,
multi_train_correlation, ingest_scada, run_synthetic_training, run_prediction.
"""
from __future__ import annotations

import json
import math
import os
import pathlib
import time
from typing import Optional

import numpy as np

# ── Hyper-parameters ────────────────────────────────────────────────
INPUT_DIM   = 8
HIDDEN_DIM  = 64
HORIZON     = 7          # days ahead to predict (batch model)
SEQ_LEN     = 21         # look-back window (3 weeks)
LR          = 3e-3
CLIP_NORM   = 1.0
EPOCHS      = 80
BATCH_SZ    = 32
SEED        = 42

DATA_DIR    = pathlib.Path(os.environ.get("DATA_DIR", "data"))
DATA_DIR.mkdir(parents=True, exist_ok=True)
WEIGHTS_PATH   = DATA_DIR / "nevai_weights.npz"
OBS_PATH       = DATA_DIR / "nevai_observations.json"
META_PATH      = DATA_DIR / "nevai_meta.json"

rng = np.random.default_rng(SEED)


# ────────────────────────────────────────────────────────────────────
# 1.  LSTM cell — forward + backward
# ────────────────────────────────────────────────────────────────────
def sigmoid(x):
    return np.where(x >= 0,
                    1 / (1 + np.exp(-np.clip(x, -500, 500))),
                    np.exp(np.clip(x, -500, 0)) /
                    (1 + np.exp(np.clip(x, -500, 0))))


def tanh_(x):
    return np.tanh(np.clip(x, -50, 50))


def dtanh(y):
    return 1.0 - y ** 2


def dsig(y):
    return y * (1.0 - y)


class NevAIWeights:
    """All trainable parameters in one object for easy save/load."""

    def __init__(self, input_dim=INPUT_DIM, hidden_dim=HIDDEN_DIM, horizon=HORIZON):
        d, h = input_dim, hidden_dim
        scale = math.sqrt(2.0 / (d + h))
        # Gates: input, forget, output, cell
        self.Wx = rng.standard_normal((4 * h, d)) * scale
        self.Wh = rng.standard_normal((4 * h, h)) * scale
        self.b = np.zeros(4 * h)
        self.b[h:2 * h] = 1.0  # forget-gate bias initialised to 1
        # Output head: hidden -> horizon
        self.Wo = rng.standard_normal((horizon, h)) * scale * 0.1
        self.bo = np.zeros(horizon)
        self.input_dim = input_dim
        self.hidden_dim = hidden_dim
        self.horizon = horizon

    def init_adam(self):
        self.m = {k: np.zeros_like(v) for k, v in self._params().items()}
        self.v = {k: np.zeros_like(v) for k, v in self._params().items()}
        self.t = 0

    def _params(self):
        return dict(Wx=self.Wx, Wh=self.Wh, b=self.b, Wo=self.Wo, bo=self.bo)

    def adam_step(self, grads, lr=LR, b1=0.9, b2=0.999, eps=1e-8):
        self.t += 1
        for k, g in grads.items():
            self.m[k] = b1 * self.m[k] + (1 - b1) * g
            self.v[k] = b2 * self.v[k] + (1 - b2) * g * g
            mh = self.m[k] / (1 - b1 ** self.t)
            vh = self.v[k] / (1 - b2 ** self.t)
            getattr(self, k).__isub__(lr * mh / (np.sqrt(vh) + eps))

    def save(self):
        np.savez_compressed(
            WEIGHTS_PATH,
            Wx=self.Wx, Wh=self.Wh, b=self.b, Wo=self.Wo, bo=self.bo,
            meta=np.array([self.input_dim, self.hidden_dim, self.horizon]))

    @classmethod
    def load(cls):
        d = np.load(WEIGHTS_PATH, allow_pickle=False)
        meta = d['meta']
        w = cls(int(meta[0]), int(meta[1]), int(meta[2]))
        w.Wx = d['Wx']; w.Wh = d['Wh']; w.b = d['b']
        w.Wo = d['Wo']; w.bo = d['bo']
        return w


# ────────────────────────────────────────────────────────────────────
# 2.  Forward pass (one sequence, returns cache for BPTT)
# ────────────────────────────────────────────────────────────────────
def nevai_forward(w, xs):
    """xs: (T, input_dim). Returns preds (horizon,), last hidden state, cache."""
    h = w.hidden_dim
    ht = np.zeros(h)
    ct = np.zeros(h)
    cache = []
    for t in range(len(xs)):
        x = xs[t]
        z = w.Wx @ x + w.Wh @ ht + w.b
        i_g = sigmoid(z[:h])
        f_g = sigmoid(z[h:2 * h])
        o_g = sigmoid(z[2 * h:3 * h])
        g_ = tanh_(z[3 * h:])
        ct_ = f_g * ct + i_g * g_
        ht_ = o_g * tanh_(ct_)
        cache.append((x, ht, ct, i_g, f_g, o_g, g_, ct_, ht_))
        ht, ct = ht_, ct_
    preds = w.Wo @ ht + w.bo
    return preds, ht, cache


# ────────────────────────────────────────────────────────────────────
# 3.  BPTT + gradient clipping
# ────────────────────────────────────────────────────────────────────
def nevai_backward(w, xs, targets, cache, preds):
    h = w.hidden_dim
    loss = 0.5 * np.sum((preds - targets) ** 2)
    d_pred = preds - targets
    dWo = np.outer(d_pred, cache[-1][-1])
    dbo = d_pred.copy()
    dht = w.Wo.T @ d_pred
    dWx = np.zeros_like(w.Wx)
    dWh = np.zeros_like(w.Wh)
    db = np.zeros_like(w.b)
    dct = np.zeros(h)
    for t in reversed(range(len(xs))):
        x, ht_prev, ct_prev, i_g, f_g, o_g, g_, ct_, ht_ = cache[t]
        d_ht = dht
        d_ot = d_ht * tanh_(ct_)
        d_ct = d_ht * o_g * dtanh(tanh_(ct_)) + dct
        d_ft = d_ct * ct_prev
        d_it = d_ct * g_
        d_g = d_ct * i_g
        d_ct_prev = d_ct * f_g
        dz = np.concatenate([
            d_it * dsig(i_g),
            d_ft * dsig(f_g),
            d_ot * dsig(o_g),
            d_g * dtanh(g_)
        ])
        dWx += np.outer(dz, x)
        dWh += np.outer(dz, ht_prev)
        db += dz
        dht = w.Wh.T @ dz
        dct = d_ct_prev
    grads = dict(Wx=dWx, Wh=dWh, b=db, Wo=dWo, bo=dbo)
    total_norm = math.sqrt(sum(float(np.sum(g ** 2)) for g in grads.values()))
    if total_norm > CLIP_NORM:
        scale = CLIP_NORM / (total_norm + 1e-8)
        grads = {k: v * scale for k, v in grads.items()}
    return loss, grads


# ────────────────────────────────────────────────────────────────────
# 4.  Data preparation from a batch digital-twin schedule
# ────────────────────────────────────────────────────────────────────
def schedule_to_features(schedule, num_trains, cfg):
    """Convert a build_schedule() output list into (N, INPUT_DIM) array."""
    if not schedule:
        return np.zeros((0, INPUT_DIM))
    peak = max(d['tot'] for d in schedule) or 1.0
    tc = cfg.get('total_batch_cycle', 17.5)
    rows = []
    for i, d in enumerate(schedule):
        norm_flow = d['tot'] / peak
        usp_frac = d['usp'] / (d['tot'] + 1e-6)
        dsp_frac = d['dsp'] / (d['tot'] + 1e-6)
        hv_flag = 1.0 if d['hv'] > 0 else 0.0
        co_flag = float(d.get('co', 0))
        r_labels = d.get('r', []) or []
        phase = 1.0 if any(isinstance(l, str) and l == 'CIP' for l in r_labels) else 0.0
        cycle_pos = (i % max(int(tc), 1)) / max(int(tc), 1)
        nt_norm = min(num_trains / 6.0, 1.0)
        rows.append([norm_flow, usp_frac, dsp_frac, hv_flag, phase, co_flag, cycle_pos, nt_norm])
    return np.array(rows, dtype=np.float32)


def build_sequences(features, targets_raw, seq_len=SEQ_LEN, horizon=HORIZON):
    """Slide a window over features to build (X, y) pairs."""
    X, y = [], []
    peak = targets_raw.max() if targets_raw.max() > 0 else 1.0
    targets = targets_raw / peak
    for i in range(len(features) - seq_len - horizon + 1):
        X.append(features[i:i + seq_len])
        y.append(targets[i + seq_len: i + seq_len + horizon])
    if not X:
        return np.zeros((0, seq_len, features.shape[1] if features.ndim == 2 else INPUT_DIM)), np.zeros((0, horizon))
    return np.array(X, dtype=np.float32), np.array(y, dtype=np.float32)


# ────────────────────────────────────────────────────────────────────
# 5.  Training loop
# ────────────────────────────────────────────────────────────────────
def train(w, X, y, epochs=EPOCHS, batch_sz=BATCH_SZ, lr=LR, verbose=False):
    if len(X) == 0:
        return []
    w.init_adam()
    n = len(X)
    losses = []
    for ep in range(epochs):
        idx = np.random.permutation(n)
        ep_loss = 0.0
        for b in range(0, n, batch_sz):
            batch_idx = idx[b:b + batch_sz]
            bl = 0.0
            grad_acc = {k: np.zeros_like(getattr(w, k)) for k in ('Wx', 'Wh', 'b', 'Wo', 'bo')}
            for i in batch_idx:
                preds, _, cache = nevai_forward(w, X[i])
                loss, grads = nevai_backward(w, X[i], y[i], cache, preds)
                bl += loss
                for k in grad_acc:
                    grad_acc[k] += grads[k] / len(batch_idx)
            w.adam_step(grad_acc, lr=lr)
            ep_loss += bl / len(batch_idx)
        losses.append(float(ep_loss))
        if verbose and ep % 10 == 0:
            print(f"Epoch {ep:3d}  loss={ep_loss:.4f}")
    return losses


# ────────────────────────────────────────────────────────────────────
# 6.  Prediction + uncertainty (weight perturbation)
# ────────────────────────────────────────────────────────────────────
def predict(w, features_window, peak_gpd, n_samples=30):
    """features_window: (SEQ_LEN, input_dim). Returns point forecast + 10–90 % band."""
    samples = []
    noise_scale = 0.02
    for _ in range(n_samples):
        w2 = NevAIWeights(w.input_dim, w.hidden_dim, w.horizon)
        w2.Wx = w.Wx + rng.standard_normal(w.Wx.shape) * noise_scale * np.std(w.Wx)
        w2.Wh = w.Wh + rng.standard_normal(w.Wh.shape) * noise_scale * np.std(w.Wh)
        w2.b = w.b.copy()
        w2.Wo = w.Wo.copy()
        w2.bo = w.bo.copy()
        preds, _, _ = nevai_forward(w2, features_window)
        samples.append(np.clip(preds, 0, None))
    samples = np.array(samples)
    mean = samples.mean(axis=0)
    lo = np.percentile(samples, 10, axis=0)
    hi = np.percentile(samples, 90, axis=0)
    return {
        "days_ahead": list(range(1, w.horizon + 1)),
        "forecast_gpd": (mean * peak_gpd).clip(0).tolist(),
        "ci_lo_gpd": (lo * peak_gpd).clip(0).tolist(),
        "ci_hi_gpd": (hi * peak_gpd).clip(0).tolist(),
    }


# ────────────────────────────────────────────────────────────────────
# 7.  Anomaly detection
# ────────────────────────────────────────────────────────────────────
def anomaly_score(w, features, peaks):
    results = []
    for i in range(SEQ_LEN, len(features)):
        window = features[i - SEQ_LEN:i]
        preds, _, _ = nevai_forward(w, window)
        actual_next = features[i, 0]
        pred_next = float(preds[0])
        err = abs(actual_next - pred_next)
        score = float(np.clip(err / 0.3, 0, 1))
        results.append({"day": i + 1, "score": round(score, 3), "flag": score > 0.5})
    return results


# ────────────────────────────────────────────────────────────────────
# 8.  Multi-train correlation matrix
# ────────────────────────────────────────────────────────────────────
def multi_train_correlation(all_schedules, all_cfgs):
    n = len(all_schedules)
    if n == 0:
        return {}
    min_len = min(len(s) for s in all_schedules)
    series = np.array([[d['tot'] for d in sch[:min_len]] for sch in all_schedules], dtype=float)
    corr = np.corrcoef(series)
    max_lag = 14
    lag_matrix = np.zeros((n, n), dtype=int)
    lag_strength = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            best_r, best_lag = -1, 0
            a = series[i] - series[i].mean()
            b = series[j] - series[j].mean()
            na = np.linalg.norm(a); nb = np.linalg.norm(b)
            if na < 1e-6 or nb < 1e-6:
                continue
            for lag in range(-max_lag, max_lag + 1):
                if lag > 0:
                    r = np.dot(a[lag:], b[:-lag]) / (na * nb)
                elif lag < 0:
                    r = np.dot(a[:lag], b[-lag:]) / (na * nb)
                else:
                    r = np.dot(a, b) / (na * nb)
                if abs(r) > abs(best_r):
                    best_r, best_lag = r, lag
            lag_matrix[i, j] = int(best_lag)
            lag_strength[i, j] = round(float(best_r), 3)
    combined_daily = series.sum(axis=0)
    harvest_flags = np.array([[1 if d['hv'] > 0 else 0 for d in sch[:min_len]] for sch in all_schedules])
    combined_harvests = harvest_flags.sum(axis=0)
    overlap_days = int((combined_harvests >= 2).sum())
    peak_day_idx = int(np.argmax(combined_daily))
    peak_breakdown = {f"train_{i + 1}": round(float(series[i, peak_day_idx]), 1) for i in range(n)}
    p2a = combined_daily.max() / (combined_daily.mean() + 1e-6)
    eq_risk = round(float(np.clip((p2a - 1) / 4, 0, 1)), 3)
    normed = series / (np.linalg.norm(series, axis=1, keepdims=True) + 1e-8)
    cos_sim = (normed @ normed.T).tolist()
    harvest_density = [round(float(harvest_flags[i].sum()) / max(min_len, 1), 4) for i in range(n)]
    return {
        "num_trains": n,
        "sim_days": min_len,
        "corr_matrix": [[round(float(corr[i, j]), 3) for j in range(n)] for i in range(n)] if n > 1 else [[1.0]],
        "cos_sim_matrix": [[round(float(cos_sim[i][j]), 3) for j in range(n)] for i in range(n)],
        "lag_days_matrix": lag_matrix.tolist(),
        "lag_strength_matrix": [[round(float(lag_strength[i, j]), 3) for j in range(n)] for i in range(n)],
        "overlap_harvest_days": overlap_days,
        "peak_day": peak_day_idx + 1,
        "peak_gpd_combined": round(float(combined_daily.max()), 1),
        "peak_breakdown": peak_breakdown,
        "eq_risk_index": eq_risk,
        "harvest_density": harvest_density,
        "combined_flow_series": [round(float(v), 1) for v in combined_daily.tolist()],
    }


# ────────────────────────────────────────────────────────────────────
# 9.  Incremental learning from SCADA observations
# ────────────────────────────────────────────────────────────────────
def ingest_scada(rows, w, peak_gpd):
    obs = []
    if OBS_PATH.exists():
        try:
            obs = json.loads(OBS_PATH.read_text())
        except Exception:
            obs = []
    valid = []
    for r in rows:
        try:
            valid.append({
                "ts": str(r.get("ts", "")),
                "tot_gpd": float(r.get("tot_gpd", 0)),
                "usp_gpd": float(r.get("usp_gpd", 0)),
                "dsp_gpd": float(r.get("dsp_gpd", 0)),
                "harvest_flag": int(bool(r.get("harvest_flag", 0))),
                "cip_phase": int(r.get("cip_phase", 0)),
                "train_id": int(r.get("train_id", 1)),
            })
        except (ValueError, TypeError):
            continue
    obs.extend(valid)
    obs = obs[-730:]
    OBS_PATH.write_text(json.dumps(obs, indent=None))
    if len(obs) < SEQ_LEN + HORIZON + 5:
        return {"accepted": len(valid), "fine_tuned": False,
                "reason": f"Not enough observations yet (need {SEQ_LEN + HORIZON + 5}, have {len(obs)})"}
    feats = []
    for i, o in enumerate(obs):
        tot = o['tot_gpd']
        feats.append([
            tot / (peak_gpd + 1e-6),
            o['usp_gpd'] / (tot + 1e-6),
            o['dsp_gpd'] / (tot + 1e-6),
            float(o['harvest_flag']),
            float(min(o['cip_phase'], 2) / 2.0),
            0.0,
            (i % 17) / 17.0,
            0.5,
        ])
    feats_arr = np.array(feats, dtype=np.float32)
    targets_raw = np.array([o['tot_gpd'] for o in obs], dtype=np.float32)
    X, y = build_sequences(feats_arr, targets_raw)
    if len(X) == 0:
        return {"accepted": len(valid), "fine_tuned": False, "reason": "Sequences too short after windowing"}
    losses = train(w, X, y, epochs=5, lr=LR * 0.3)
    w.save()
    return {"accepted": len(valid), "fine_tuned": True, "obs_total": len(obs),
            "ft_loss_last": round(losses[-1], 5) if losses else None}


# ────────────────────────────────────────────────────────────────────
# 10. Model manager (batch model)
# ────────────────────────────────────────────────────────────────────
_MODEL: Optional[NevAIWeights] = None
_META: dict = {}


def get_model():
    global _MODEL
    if _MODEL is None:
        if WEIGHTS_PATH.exists():
            try:
                _MODEL = NevAIWeights.load()
                return _MODEL
            except Exception:
                pass
        _MODEL = NevAIWeights()
        _MODEL.init_adam()
    return _MODEL


def get_meta():
    if META_PATH.exists():
        try:
            return json.loads(META_PATH.read_text())
        except Exception:
            pass
    return {"trained": False, "epochs": 0, "samples": 0, "last_trained": None, "scada_obs": 0}


def save_meta(meta):
    META_PATH.write_text(json.dumps(meta, indent=2))


# ────────────────────────────────────────────────────────────────────
# 11. High-level API for the batch wastewater app
# ────────────────────────────────────────────────────────────────────
def run_synthetic_training(all_schedules, all_cfgs, force=False):
    meta = get_meta()
    if meta.get("trained") and not force:
        return meta
    w = get_model()
    all_X, all_y = [], []
    for sch, cfg in zip(all_schedules, all_cfgs):
        if len(sch) < SEQ_LEN + HORIZON + 2:
            continue
        num_trains = int(cfg.get("num_trains", 1))
        feats = schedule_to_features(sch, num_trains, cfg)
        targets_raw = np.array([d['tot'] for d in sch], dtype=np.float32)
        X, y = build_sequences(feats, targets_raw)
        if len(X):
            all_X.append(X)
            all_y.append(y)
    if not all_X:
        return {"trained": False, "error": "No usable schedule data"}
    X_all = np.concatenate(all_X, axis=0)
    y_all = np.concatenate(all_y, axis=0)
    idx = np.random.permutation(len(X_all))
    X_all, y_all = X_all[idx], y_all[idx]
    losses = train(w, X_all, y_all, epochs=EPOCHS, verbose=False)
    w.save()
    meta = {
        "trained": True,
        "epochs": EPOCHS,
        "samples": int(len(X_all)),
        "last_trained": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "final_loss": round(losses[-1], 5) if losses else None,
        "scada_obs": meta.get("scada_obs", 0),
    }
    save_meta(meta)
    return meta


def run_prediction(schedule, cfg, num_trains):
    w = get_model()
    feats = schedule_to_features(schedule, num_trains, cfg)
    if len(feats) < SEQ_LEN:
        return {"error": f"Need at least {SEQ_LEN} days of schedule"}
    window = feats[-SEQ_LEN:]
    peak_gpd = max(d['tot'] for d in schedule) or 1.0
    return predict(w, window, peak_gpd)


def run_anomaly(schedule, cfg, num_trains):
    w = get_model()
    feats = schedule_to_features(schedule, num_trains, cfg)
    if len(feats) < SEQ_LEN + 2:
        return []
    peaks = np.array([d['tot'] for d in schedule], dtype=float)
    return anomaly_score(w, feats, peaks)
