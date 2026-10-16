"""
icb_routes.py — Flask routes for the IDB digital twin
NavNevAI_mAbWasteWater · CarbonTesseract 2026

  GET  /api/icb/defaults        default parameters, design check, titre curves
  POST /api/icb/design          design check for {"params": {...}}
  POST /api/icb/portfolio       buffer / sump sizing across all molecules
  POST /api/icb/simulate        dynamic run {"params", "seed", "downsample_h", "use_forecaster"}
  POST /api/icb/forecast        run, then NevAI 24 h bands for titre and wastewater
  GET  /api/icb/nevai/status    model state + background training state
  POST /api/icb/nevai/train     start background training (returns immediately)
  POST /api/icb/nevai/ingest    hourly plant rows -> fine-tune both forecasters

Training runs in a background thread. A lock file stops several gunicorn
workers from training at once. With AUTO_TRAIN=1 (default) training starts
automatically when no current-version model exists.
"""
import os
import threading
import time
from flask import jsonify, request
from icb_process import (ICBParams, WASTE_STREAMS, BUFFER_OPS, TITRE_PROFILES,
                         process_design, portfolio_sizing)
from icb_dynamic import ICBSimulator
import icb_nevai

_TRAIN = {"running": False, "started": None, "finished": None, "last": None, "error": None}
_TRAIN_LOCK = threading.Lock()
_LOCKFILE = icb_nevai.ICB_DIR / ".training.lock"
_LOCK_STALE_S = 7200


def _now():
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def _lock_held():
    try:
        return _LOCKFILE.exists() and time.time() - _LOCKFILE.stat().st_mtime < _LOCK_STALE_S
    except OSError:
        return False


def _lock_claim():
    try:
        if _LOCKFILE.exists() and not _lock_held():
            _LOCKFILE.unlink()
        fd = os.open(str(_LOCKFILE), os.O_CREAT | os.O_EXCL | os.O_WRONLY)
        os.write(fd, str(os.getpid()).encode())
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


def register_icb_routes(app, login_required=None):
    def _guard(fn):
        return login_required(fn) if callable(login_required) else fn

    def _params(body):
        p = ICBParams.from_dict(body.get("params"))
        p.run_days = max(3.0, min(float(p.run_days), 120.0))
        return p

    def _forecaster(body):
        if not body.get("use_forecaster", True) or _TRAIN["running"] or _lock_held():
            return None
        fs = icb_nevai.ICBForecasters()
        return fs if fs.load() else None

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

    @app.route("/api/icb/design", methods=["POST"])
    @_guard
    def icb_design():
        body = request.get_json(force=True) or {}
        return jsonify({"ok": True, "design": process_design(_params(body))})

    @app.route("/api/icb/portfolio", methods=["POST"])
    @_guard
    def icb_portfolio():
        try:
            body = request.get_json(force=True) or {}
            return jsonify({"ok": True, **portfolio_sizing(_params(body))})
        except Exception as e:
            app.logger.exception("icb portfolio error")
            return jsonify({"ok": False, "error": str(e)}), 400

    @app.route("/api/icb/simulate", methods=["POST"])
    @_guard
    def icb_simulate():
        try:
            body = request.get_json(force=True) or {}
            fc = _forecaster(body)
            res = ICBSimulator(_params(body), forecaster=fc, seed=int(body.get("seed", 0))).run(
                downsample_h=int(body.get("downsample_h", 1)))
            res["forecaster_used"] = fc is not None
            return jsonify(res)
        except Exception as e:
            app.logger.exception("icb simulate error")
            return jsonify({"ok": False, "error": str(e)}), 400

    @app.route("/api/icb/forecast", methods=["POST"])
    @_guard
    def icb_forecast():
        try:
            body = request.get_json(force=True) or {}
            fs = _forecaster({"use_forecaster": True})
            if fs is None:
                return jsonify({"ok": False, "error": "NevAI forecasters not trained yet"}), 400
            res = ICBSimulator(_params(body), forecaster=fs, seed=int(body.get("seed", 0))).run(keep_history=True)
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
        if _TRAIN["running"] or _lock_held():
            st = {"ready": False, "titre": {}, "ww": {}}
        else:
            st = icb_nevai.status()
        training = dict(_TRAIN)
        training["running"] = training["running"] or _lock_held()
        return jsonify({"ok": True, **st, "training": training})

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
            fs = _forecaster({"use_forecaster": True})
            if fs is None:
                return jsonify({"ok": False, "error": "Train the forecasters before fine-tuning"}), 400
            return jsonify({"ok": True,
                            "titre": fs.titre.fine_tune(rows),
                            "wastewater": fs.ww.fine_tune(rows),
                            "required_keys": sorted(set(icb_nevai.TITRE_KEYS + icb_nevai.WW_KEYS))})
        except Exception as e:
            app.logger.exception("icb nevai ingest error")
            return jsonify({"ok": False, "error": str(e)}), 400

    if os.environ.get("AUTO_TRAIN", "1") == "1" and not icb_nevai.ICBForecasters().load():
        _start_training({})
