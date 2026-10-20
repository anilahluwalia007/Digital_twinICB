"""
icb_routes.py — Flask routes for the IDB digital twin
NavNevAI_mAbWasteWater · CarbonTesseract 2026

  GET  /api/icb/defaults        default parameters, design check, titre curves
  GET  /api/icb/cases           the five basis-of-design cases and their default settings
  POST /api/icb/cases/preview   noise-free plan for every case {"params", "case": {"settings", "values"}}
  POST /api/icb/cases/compare   simulate all five cases, design envelope (Cases 1–3), Case 5 vs Case 1
  POST /api/icb/design          design check for {"params": {...}}
  POST /api/icb/portfolio       buffer / sump sizing across all titre curves
  POST /api/icb/sizing          equipment sizing, worst case across all titre curves
  POST /api/icb/simulate        dynamic run {"params", "seed", "downsample_h", "use_forecaster",
                                "case": {"id", "settings", "values"} | null}
  POST /api/icb/forecast        run, then NevAI 24 h bands for titre and wastewater
  GET  /api/icb/nevai/status    model state, weights on disk, training state, why_off
  POST /api/icb/nevai/train     start background training (returns immediately)
  POST /api/icb/nevai/ingest    hourly plant rows -> fine-tune both forecasters

A case is applied on top of the posted parameters (icb_cases.apply_case); the run returns
result["case"] = {"info": what the case changed, "annual": per-year totals, occupancy and checks}.

Training runs in a background thread. A lock file stops several gunicorn
workers from training at once; a lock left by a worker that died is ignored.
Weights are saved atomically, so simulations keep using the current weights
while a retraining runs. When NevAI is not used, the response says why
(forecaster_reason). With AUTO_TRAIN=1 (default) training starts automatically
when no current-version model exists.
"""
import copy
import os
import socket
import threading
import time
from flask import jsonify, request
from icb_process import (ICBParams, WASTE_STREAMS, BUFFER_OPS, TITRE_PROFILES,
                         process_design, portfolio_sizing)
from icb_dynamic import ICBSimulator
from icb_sizing import size_portfolio
from icb_cases import (apply_case, annual, case_api_info, preview_cases, compare_cases, rebase_suggestions)
import icb_nevai

_TRAIN = {"running": False, "started": None, "finished": None, "last": None, "error": None}
_TRAIN_LOCK = threading.Lock()
_LOCKFILE = icb_nevai.ICB_DIR / ".training.lock"
_LOCK_STALE_S = 7200
_HOST = socket.gethostname()


def _now():
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def _pid_alive(pid):
    try:
        os.kill(pid, 0)
        return True
    except ProcessLookupError:
        return False
    except Exception:
        return True          # exists but not ours, or cannot tell: assume alive


def _lock_owner():
    try:
        host, _, pid = _LOCKFILE.read_text().strip().rpartition(":")
        return host, int(pid)
    except (OSError, ValueError):
        return None, None


def _lock_held():
    try:
        if not _LOCKFILE.exists():
            return False
        if time.time() - _LOCKFILE.stat().st_mtime >= _LOCK_STALE_S:
            return False
        host, pid = _lock_owner()
        if pid and host in ("", _HOST) and not _pid_alive(pid):
            return False     # left behind by a worker that died mid-training
        return True
    except OSError:
        return False


def _lock_claim():
    try:
        if _LOCKFILE.exists() and not _lock_held():
            try:
                _LOCKFILE.unlink()
            except OSError:
                pass
        fd = os.open(str(_LOCKFILE), os.O_CREAT | os.O_EXCL | os.O_WRONLY)
        os.write(fd, f"{_HOST}:{os.getpid()}".encode())
        os.close(fd)
        return True
    except FileExistsError:
        return False
    except OSError:
        return True   # no shared file system: fall back to the in-process lock


def _lock_release():
    try:
        _LOCKFILE.unlink()
    except OSError:
        pass


def _start_training(kwargs):
    with _TRAIN_LOCK:
        if _TRAIN["running"] or not _lock_claim():
            return False
        _TRAIN.update(running=True, started=_now(), finished=None, error=None)

    def job():
        try:
            _TRAIN["last"] = icb_nevai.train_forecasters(**kwargs)
        except Exception as e:
            _TRAIN["error"] = str(e)
        finally:
            _lock_release()
            _TRAIN["running"] = False
            _TRAIN["finished"] = _now()

    threading.Thread(target=job, daemon=True).start()
    return True


def _why_off():
    """Plain-language reason why NevAI cannot be used right now."""
    rep = icb_nevai.weights_report()
    missing = [n for n in ("titre", "ww") if not (rep[n]["weights_file"] and rep[n]["meta_file"])]
    old = [n for n in ("titre", "ww") if n not in missing
           and rep[n]["feature_version"] != icb_nevai.FEATURE_VERSION]
    if _TRAIN["running"] or _lock_held():
        return "NevAI is training on the server; forecasts switch on when it finishes"
    if missing:
        return ("no trained weights in data/icb/ on the server (missing: " + ", ".join(missing) +
                "); check data/icb/ is in the GitHub repo and not excluded by .gitignore, or click Train NevAI")
    if old:
        return "the weights in data/icb/ are from an older version; click Train NevAI"
    if _TRAIN.get("error"):
        return "the last training failed: " + _TRAIN["error"]
    return "the weights in data/icb/ could not be loaded"


def _forecaster_or_reason(use=True):
    if not use:
        return None, "NevAI forecasts were unticked for this run"
    fs = icb_nevai.ICBForecasters()
    try:
        if fs.load():
            return fs, None
    except Exception as e:
        return None, f"the weights could not be loaded ({e})"
    return None, _why_off()


def register_icb_routes(app, login_required=None):
    def _guard(fn):
        return login_required(fn) if callable(login_required) else fn

    def _params(body):
        p = ICBParams.from_dict(body.get("params"))
        p.run_days = max(3.0, min(float(p.run_days), 120.0))
        return p

    def _case_params(body):
        """Inputs + the selected case (if any): (case-applied params, Inputs params, case info)."""
        p = _params(body)
        base = copy.deepcopy(p)
        info = apply_case(p, body.get("case"))
        return p, base, info

    @app.route("/api/icb/defaults", methods=["GET"])
    @_guard
    def icb_defaults():
        p = ICBParams()
        return jsonify({"ok": True, "params": p.to_dict(), "design": process_design(p),
                        "profiles": {k: {"label": v["label"], "points": v["points"]}
                                     for k, v in TITRE_PROFILES.items()},
                        "buffer_ops": list(BUFFER_OPS),
                        "waste_streams": {k: {"area": a, "route": r, "driver": dr, "desc": d}
                                          for k, (a, r, dr, d) in WASTE_STREAMS.items()}})

    @app.route("/api/icb/cases", methods=["GET"])
    @_guard
    def icb_cases():
        return jsonify({"ok": True, **case_api_info()})

    @app.route("/api/icb/cases/preview", methods=["POST"])
    @_guard
    def icb_cases_preview():
        try:
            body = request.get_json(force=True) or {}
            return jsonify({"ok": True, **preview_cases(_params(body), body.get("case"))})
        except Exception as e:
            app.logger.exception("icb case preview error")
            return jsonify({"ok": False, "error": str(e)}), 400

    @app.route("/api/icb/cases/compare", methods=["POST"])
    @_guard
    def icb_cases_compare():
        try:
            body = request.get_json(force=True) or {}
            return jsonify(compare_cases(_params(body), body.get("case"), int(body.get("seed", 0)),
                                         body.get("ids")))
        except Exception as e:
            app.logger.exception("icb case compare error")
            return jsonify({"ok": False, "error": str(e)}), 400

    @app.route("/api/icb/design", methods=["POST"])
    @_guard
    def icb_design():
        body = request.get_json(force=True) or {}
        p, _, _ = _case_params(body)
        return jsonify({"ok": True, "design": process_design(p)})

    @app.route("/api/icb/portfolio", methods=["POST"])
    @_guard
    def icb_portfolio():
        try:
            body = request.get_json(force=True) or {}
            return jsonify({"ok": True, **portfolio_sizing(_params(body))})
        except Exception as e:
            app.logger.exception("icb portfolio error")
            return jsonify({"ok": False, "error": str(e)}), 400

    @app.route("/api/icb/sizing", methods=["POST"])
    @_guard
    def icb_sizing():
        try:
            body = request.get_json(force=True) or {}
            return jsonify({"ok": True, **size_portfolio(_params(body), int(body.get("seed", 0)))})
        except Exception as e:
            app.logger.exception("icb sizing error")
            return jsonify({"ok": False, "error": str(e)}), 400

    @app.route("/api/icb/simulate", methods=["POST"])
    @_guard
    def icb_simulate():
        try:
            body = request.get_json(force=True) or {}
            p, base, info = _case_params(body)
            fc, why = _forecaster_or_reason(bool(body.get("use_forecaster", True)))
            res = ICBSimulator(p, forecaster=fc, seed=int(body.get("seed", 0))).run(
                downsample_h=int(body.get("downsample_h", 1)))
            res["forecaster_used"] = fc is not None
            if fc is None:
                res["forecaster_reason"] = why
            res["case"] = {"info": info, "annual": annual(p, res, info)}
            if info:
                rebase_suggestions((res.get("sizing") or {}).get("suggested_params"), base, p)
            return jsonify(res)
        except Exception as e:
            app.logger.exception("icb simulate error")
            return jsonify({"ok": False, "error": str(e)}), 400

    @app.route("/api/icb/forecast", methods=["POST"])
    @_guard
    def icb_forecast():
        try:
            body = request.get_json(force=True) or {}
            fs, why = _forecaster_or_reason(True)
            if fs is None:
                return jsonify({"ok": False, "error": "NevAI not available: " + why}), 400
            p, _, _ = _case_params(body)
            res = ICBSimulator(p, forecaster=fs, seed=int(body.get("seed", 0))).run(keep_history=True)
            hist = res.pop("_history")
            n = int(body.get("n_samples", 30))
            return jsonify({
                "ok": True,
                "titre": fs.titre.forecast_band(hist, n),
                "wastewater": fs.ww.forecast_band(hist, n),
                "recent": [{"h": r["h"], "titre_meas": round(r["titre_meas"] or 0.0, 4),
                            "ww_total_Lh": round(r["ww_total_Lh"], 2)} for r in hist[-72:]],
                "kpi": res["kpi"],
            })
        except Exception as e:
            app.logger.exception("icb forecast error")
            return jsonify({"ok": False, "error": str(e)}), 400

    @app.route("/api/icb/nevai/status", methods=["GET"])
    @_guard
    def icb_nevai_status():
        st = icb_nevai.status()
        training = dict(_TRAIN)
        training["running"] = bool(training["running"] or _lock_held())
        out = {"ok": True, **st, "training": training, "weights": icb_nevai.weights_report()}
        if not st.get("ready"):
            out["why_off"] = _why_off()
        return jsonify(out)

    @app.route("/api/icb/nevai/train", methods=["POST"])
    @_guard
    def icb_nevai_train():
        body = request.get_json(force=True) or {}
        kwargs = {
            "base_params": body.get("params"),
            "n_runs": max(2, min(int(body.get("n_runs", 4)), 8)),
            "epochs": max(3, min(int(body.get("epochs", 8)), 40)),
            "stride": max(2, int(body.get("stride", 12))),
        }
        started = _start_training(kwargs)
        return jsonify({"ok": True, "started": started, "training": dict(_TRAIN)})

    @app.route("/api/icb/nevai/ingest", methods=["POST"])
    @_guard
    def icb_nevai_ingest():
        try:
            body = request.get_json(force=True) or {}
            rows = body.get("rows") or []
            fs, why = _forecaster_or_reason(True)
            if fs is None:
                return jsonify({"ok": False, "error": "Train the forecasters before fine-tuning: " + why}), 400
            return jsonify({"ok": True,
                            "titre": fs.titre.fine_tune(rows),
                            "wastewater": fs.ww.fine_tune(rows),
                            "required_keys": sorted(set(icb_nevai.TITRE_KEYS + icb_nevai.WW_KEYS))})
        except Exception as e:
            app.logger.exception("icb nevai ingest error")
            return jsonify({"ok": False, "error": str(e)}), 400

    if _LOCKFILE.exists() and not _lock_held():          # clear a lock left by a dead worker
        _lock_release()
    if os.environ.get("AUTO_TRAIN", "1") == "1" and not icb_nevai.ICBForecasters().load():
        _start_training({})
