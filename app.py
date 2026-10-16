"""
app.py — web app for the Integrated Dynamic BioProcessing (IDB) mAb digital twin
NavNevAI_mAbWasteWater · CarbonTesseract 2026

Local (Mac):  python3 app.py   ->  open http://127.0.0.1:8000
Render:       gunicorn app:app (see render.yaml)

Login
  Set ADMIN_EMAIL and ADMIN_PW (Render: Environment tab). The app then shows a
  login page; a signed session cookie keeps you logged in for SESSION_HOURS
  (default 12). The dashboard has a Log out button.
  SECRET_KEY signs the cookie; render.yaml generates one. Without it a random key
  is used, so everyone is logged out whenever the server restarts.
  Locally, with neither ADMIN_EMAIL nor ADMIN_PW set, the app runs without login.
  On Render it refuses to open until both are set.
"""
import hmac
import os
import secrets
import threading
import time
from datetime import timedelta
from urllib.parse import urlparse

from flask import (Flask, jsonify, redirect, render_template_string, request,
                   session, url_for)
from werkzeug.middleware.proxy_fix import ProxyFix
from icb_routes import register_icb_routes

app = Flask(__name__, static_folder="static", static_url_path="/static")
app.wsgi_app = ProxyFix(app.wsgi_app, x_for=1, x_proto=1)   # Render terminates HTTPS in front of the app

ON_RENDER = bool(os.environ.get("RENDER"))
app.config.update(
    SECRET_KEY=os.environ.get("SECRET_KEY") or secrets.token_hex(32),
    SESSION_COOKIE_HTTPONLY=True,
    SESSION_COOKIE_SAMESITE="Lax",            # blocks cross-site POSTs carrying the cookie
    SESSION_COOKIE_SECURE=ON_RENDER,          # HTTPS-only cookie on Render; plain http works locally
    PERMANENT_SESSION_LIFETIME=timedelta(hours=float(os.environ.get("SESSION_HOURS", "12"))),
)

PUBLIC_PATHS = {"/login", "/healthz"}
MAX_FAILS, LOCK_S = 5, 15 * 60                # 5 wrong passwords per IP -> 15 min lock-out
_fails = {}
_fails_lock = threading.Lock()


def _creds():
    return (os.environ.get("ADMIN_EMAIL", "").strip().lower(), os.environ.get("ADMIN_PW", ""))


def _auth_enabled():
    email, pw = _creds()
    return bool(email and pw) or ON_RENDER    # on Render, missing credentials lock the app


def _signed_in():
    """A session counts only for the current ADMIN_EMAIL (changing it signs everyone out)."""
    user = session.get("user")
    return bool(user) and user == _creds()[0]


def _safe_next(target):
    """Only allow redirects back into this app."""
    if not target:
        return "/"
    u = urlparse(target)
    if u.scheme or u.netloc or not target.startswith("/") or target.startswith("//"):
        return "/"
    return target


def _locked(ip):
    with _fails_lock:
        n, since = _fails.get(ip, (0, 0.0))
        if n >= MAX_FAILS and time.time() - since < LOCK_S:
            return int(LOCK_S - (time.time() - since))
        if n >= MAX_FAILS:
            _fails.pop(ip, None)
    return 0


def _record_fail(ip):
    with _fails_lock:
        n, since = _fails.get(ip, (0, time.time()))
        _fails[ip] = (n + 1, since if n else time.time())


LOGIN_HTML = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Sign in — IDB digital twin</title>
<style>
:root{--bg:#f4f6f9;--card:#fff;--bd:#e3e7ee;--t1:#16202c;--t2:#5b6675;--ac:#0e7c7b;--er:#b42318}
*{box-sizing:border-box}
body{margin:0;min-height:100vh;display:flex;align-items:center;justify-content:center;
     font:14px/1.45 -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,sans-serif;background:var(--bg);color:var(--t1)}
.card{width:min(380px,calc(100% - 32px));background:var(--card);border:1px solid var(--bd);border-radius:12px;padding:26px 24px}
h1{font-size:18px;margin:0 0 4px}
p{margin:0 0 18px;color:var(--t2);font-size:13px}
label{display:block;font-size:12px;color:var(--t2);margin:12px 0 4px}
input{width:100%;padding:9px 10px;border:1px solid var(--bd);border-radius:8px;font-size:14px}
input:focus{outline:2px solid var(--ac);outline-offset:-1px}
button{width:100%;margin-top:18px;background:var(--ac);color:#fff;border:0;border-radius:8px;padding:10px;font-weight:600;font-size:14px;cursor:pointer}
.err{color:var(--er);font-size:13px;margin-top:12px}
</style>
</head>
<body>
<form class="card" method="post" action="/login">
  <h1>IDB mAb digital twin</h1>
  <p>Sign in to continue.</p>
  <input type="hidden" name="csrf" value="{{ csrf }}">
  <input type="hidden" name="next" value="{{ next }}">
  <label for="email">Email</label>
  <input id="email" name="email" type="email" autocomplete="username" value="{{ email }}" required autofocus>
  <label for="pw">Password</label>
  <input id="pw" name="password" type="password" autocomplete="current-password" required>
  <button type="submit" {{ 'disabled' if disabled else '' }}>Sign in</button>
  {% if error %}<div class="err">{{ error }}</div>{% endif %}
</form>
</body>
</html>"""


def _login_page(error="", email="", nxt="/", status=200, disabled=False):
    session.setdefault("csrf", secrets.token_urlsafe(24))
    return render_template_string(LOGIN_HTML, error=error, email=email, next=nxt,
                                  csrf=session["csrf"], disabled=disabled), status


@app.before_request
def _require_login():
    if not _auth_enabled() or request.path in PUBLIC_PATHS or _signed_in():
        return None
    if request.path.startswith("/api/"):
        return jsonify({"ok": False, "error": "Not signed in", "login": "/login"}), 401
    return redirect(url_for("login", next=request.full_path.rstrip("?")))


@app.route("/login", methods=["GET", "POST"])
def login():
    nxt = _safe_next(request.values.get("next"))
    if not _auth_enabled() or _signed_in():
        return redirect(nxt)
    admin_email, admin_pw = _creds()
    if not (admin_email and admin_pw):
        return _login_page("Login is not configured: set ADMIN_EMAIL and ADMIN_PW in the "
                           "Render Environment tab, then redeploy.", status=503, disabled=True)
    if request.method == "GET":
        return _login_page(nxt=nxt)

    ip = request.remote_addr or "?"
    email = request.form.get("email", "").strip().lower()
    wait = _locked(ip)
    if wait:
        return _login_page(f"Too many failed attempts. Try again in {wait // 60 + 1} min.",
                           email, nxt, 429)
    if not hmac.compare_digest(request.form.get("csrf", ""), session.get("csrf", "")):
        return _login_page("Your sign-in form expired. Please try again.", email, nxt, 400)
    ok_email = hmac.compare_digest(email.encode(), admin_email.encode())
    ok_pw = hmac.compare_digest(request.form.get("password", "").encode(), admin_pw.encode())
    if not (ok_email and ok_pw):
        _record_fail(ip)
        time.sleep(0.5)
        return _login_page("Email or password is incorrect.", email, nxt, 401)

    with _fails_lock:
        _fails.pop(ip, None)
    session.clear()                           # new session id on sign-in
    session.permanent = True
    session["user"] = admin_email
    return redirect(nxt)


@app.route("/logout", methods=["POST"])
def logout():
    session.clear()
    return redirect(url_for("login"))


@app.route("/api/me")
def me():
    return jsonify({"ok": True, "auth": _auth_enabled(), "user": session.get("user")})


@app.route("/")
def index():
    return app.send_static_file("index.html")


@app.route("/healthz")
def healthz():
    return jsonify({"ok": True})


register_icb_routes(app)

if __name__ == "__main__":
    port = int(os.environ.get("PORT", "8000"))   # macOS uses 5000 for AirPlay
    app.run(host="127.0.0.1", port=port, debug=False)
