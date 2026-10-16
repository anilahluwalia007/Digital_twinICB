# Integrated Dynamic BioProcessing (IDB) — mAb digital twin

N-1 perfusion seed (media in = cell-free permeate out to waste) → N reactor
(feed in = cells + broth + product out) → TFF cell removal (cell concentrate to
bio-kill) → harvest sump → standard (batch) Protein A → low-pH viral inactivation
→ CEX/AEX polishing → virus filter → SPTFF with inline DF → drug substance.

Each molecule's permeate titre curve (Figure 2) drives the process. Titre changes
every hour, so buffer and wastewater are planned in litres per gram recovered, not
per day. The twin includes:

- four digitised titre curves (purple, green, blue, red) plus a custom curve
- one Protein A cycle each time the sump holds one column load (3-column PCC is an option)
- buffer per unit operation, per day and in L/g
- wastewater by stream, route and driver (volume, mass, event)
- NevAI (LSTM) 24 h forecasts of titre, grams, buffer and wastewater
- portfolio sizing across all molecules (which molecule sets buffer supply, which sets the sump)
- the old ATF perfusion model (cells returned to the reactor), kept for comparison

Default parameters are illustrative, and the titre curves were read off a chart
image. Replace both with site data before relying on absolute numbers.

Trained NevAI weights (feature version 2) are included in `data/icb/`, so the app
shows forecasts immediately, locally and on Render. Weights from the older ATF/PCC
version are ignored automatically and retrained.

## Run on a Mac

    python3 -m venv .venv
    source .venv/bin/activate
    pip install -r requirements.txt
    python3 icb_dynamic.py          # quick check: prints design + KPIs
    python3 app.py                  # open http://127.0.0.1:8000

Optional: `python3 icb_nevai.py train` retrains NevAI (8 runs × 30 epochs, about
2 minutes) and overwrites the weights in `data/icb/`.

Port 8000 is used because macOS reserves 5000 for AirPlay.

## Using your own titre data

- Quick test: pick **Custom curve** in the Molecule menu, then paste
  `[[N day, g/L], ...]` into the box that opens.
- Permanent: replace the points in `TITRE_PROFILES` in `icb_process.py`, then
  retrain NevAI.

The titre is interpolated linearly between points and held at the last value after
the curve ends.

## Put it on GitHub

Create an empty repository on github.com (no README, no .gitignore), for example
`icb-mab-twin`. Then use either option.

### Option A — no terminal (drag and drop)

1. Unzip `icb-mab-twin.zip`.
2. In Finder, open the `icb-mab-twin` folder and press **Cmd + Shift + .** so the
   hidden files (`.gitignore`, `.python-version`) show.
3. On the new repository page, click **uploading an existing file**
   (or **Add file → Upload files**).
4. Select everything inside the folder (Cmd + A) and drag it onto the page.
   The `static` and `data` folders must be included. If folders don't upload,
   use Chrome.
5. Click **Commit changes**.

Check the repository shows `app.py`, `render.yaml`, `requirements.txt`,
`static/index.html`, `static/chart.umd.js` and `data/icb/` with four files.

### Option B — terminal

    cd icb-mab-twin
    git init
    git add .
    git commit -m "IDB dynamic twin"
    git branch -M main
    git remote add origin https://github.com/YOUR-USER/icb-mab-twin.git
    git push -u origin main

## Deploy on Render

1. render.com → **New → Blueprint** → connect GitHub if asked → pick the repository.
   Render reads `render.yaml`.
2. Enter your sign-in details when asked: `ADMIN_EMAIL` (your email) and `ADMIN_PW`
   (a strong password). `SECRET_KEY` is generated for you. Render only asks during
   this first Blueprint setup; to change them later, edit them in the service's
   **Environment** tab and save (the service restarts).
3. Click **Apply** / **Deploy**, wait for the build, then open the service URL and
   sign in. The **Log out** button is at the top right of the dashboard.

Already deployed with the old `APP_PASSWORD` pop-up? In the service's **Environment**
tab, delete `APP_PASSWORD`, add `ADMIN_EMAIL` and `ADMIN_PW`, and add `SECRET_KEY`
with any long random string (or click **Generate**). Then push the new files.

Sign-in details:
- You stay signed in for 12 hours (`SESSION_HOURS`). Changing `ADMIN_EMAIL` or
  `SECRET_KEY` signs everyone out.
- Five wrong passwords from one address lock sign-in from it for 15 minutes.
- On Render the app stays locked until both `ADMIN_EMAIL` and `ADMIN_PW` are set.
  Locally, leaving both unset runs the app without a login; set them in the
  terminal to test it, e.g. `ADMIN_EMAIL=me@x.com ADMIN_PW=test python3 app.py`.

Every later push to `main` redeploys automatically.

Notes:
- The free instance sleeps when idle, so the first visit after a break takes up to
  a minute to wake.
- The free instance's disk resets on restart. Weights committed in `data/icb/`
  are always there; anything trained or fine-tuned on Render is lost on restart.
- If `data/icb/` is missing or holds old-version weights, NevAI retrains in the
  background (`AUTO_TRAIN=1`), which is slow on the free CPU. Simulations still run
  meanwhile, just without forecasts. A lock file stops several workers training at once.

## API

    GET  /api/icb/defaults        default parameters, design check, titre curves
    POST /api/icb/design          design check for {"params": {...}}
    POST /api/icb/portfolio       buffer / sump sizing across all molecules
    POST /api/icb/simulate        dynamic run {"params", "seed", "downsample_h", "use_forecaster"}
    POST /api/icb/forecast        run, then NevAI 24 h bands for titre and wastewater
    GET  /api/icb/nevai/status    model state + background training state
    POST /api/icb/nevai/train     start background training
    POST /api/icb/nevai/ingest    hourly plant rows -> fine-tune both forecasters

## Files

- `icb_process.py` — parameters, titre curves, waste streams, design check, portfolio sizing
- `icb_dynamic.py` — hour-by-hour simulation: N-1, N + TFF, sump, capture, VI, polishing, buffer
- `icb_nevai.py` — NevAI titre and wastewater forecasters (hourly, 24 h ahead)
- `nevai_engine.py` — NevAI LSTM core (Python 3.9+ compatible)
- `icb_routes.py` — API routes; background training
- `app.py` — web app entry point; login page, sessions and logout
- `static/index.html` — dashboard
- `static/chart.umd.js` — Chart.js 4.4.1 (MIT), bundled so no CDN is needed
- `data/icb/` — trained NevAI weights
- `render.yaml`, `requirements.txt`, `.python-version` — Render deployment
