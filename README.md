# Integrated Dynamic BioProcessing (IDB) — mAb digital twin

N-1 perfusion seed (media in = cell-free permeate out to waste) → N reactor
(feed in = cells + broth + product out) → TFF cell removal (cell concentrate to
bio-kill) → conditioned media surge tank → Brij-9 in-line → Protein A (columns
1&2, 3&4) → depth filtration → Tris base/NaCl in-line → AEX → acetic acid in-line
→ CEX → virus filtration (2 holders) → UF/DF → formulation spike → BDS hold →
0.2 µm → cryovault. Every downstream step runs in sequence (Case 5 flowsheet).

Each molecule's permeate titre curve (Figure 2) drives the process. Titre changes
every hour, so buffer and wastewater are planned in litres per gram recovered, not
per day. The twin includes:

- four digitised titre curves (purple, green, blue, red) plus a custom curve
- one Protein A cycle each time the sump holds one column load (3-column PCC is an option)
- every downstream step of every unit, with each named solution in litres (WFI
  included) and the raw chemicals in kg
- buffer per unit operation, per day and in L/g
- wastewater by stream, route and driver (volume, mass, event)
- NevAI (LSTM) 24 h forecasts of titre, grams, buffer and wastewater
- portfolio sizing across all molecules (which molecule sets buffer supply, which sets the sump)
- a live process flow diagram that follows the site flowsheet (streams 1–63) and
  plays the last run hour by hour, with a stream table and SVG download
- a start-up film and a step-sequence film, both recordable as MP4 / WebM
- a one-click Print / Save as PDF report (landscape: diagram, KPIs, charts, tables,
  all parameters)
- the old ATF perfusion model (cells returned to the reactor), kept for comparison

Default parameters are illustrative, and the titre curves were read off a chart
image. Replace both with site data before relying on absolute numbers.

Trained NevAI weights (feature version 3) are included in `data/icb/`, so the app
shows forecasts immediately, locally and on Render. Weights from older versions
are ignored automatically and retrained.

## Run on a Mac

    python3 -m venv .venv
    source .venv/bin/activate
    pip install -r requirements.txt
    python3 icb_dynamic.py          # quick check: prints design, KPIs and solutions
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

## Process flow diagram and PDF report

After each run the diagram under the KPIs plays the simulation hour by hour.
Flowing lines turn teal and animate, the surge tank and Protein A SUM levels move,
and columns colour up while loading or eluting. Hover a line or stream number for
L/h, g/h and g/L; the stream table below the diagram shows the same values. Use
the slider to pick an hour, then **Print / Save as PDF** and choose "Save as PDF"
in the print dialog.

Notes:
- All 63 streams carry values. Each column, filter and UF/DF unit shows the step
  it is running; a grey dashed line has no flow in the hour shown.
- Polishing order: `polish.order` defaults to `aex_cex` (AEX then CEX, as on the
  site PFD). Set it to `cex_aex` in the Parameters panel for the legacy order.

## Downstream, step by step (Case 5 flowsheet)

`icb_dsp.py` runs every downstream unit as a recipe of steps, each with a named
solution, volume, time and destination:

- Brij-9 dosed in-line into each Protein A load (vi.mode = "detergent")
- Protein A, 2 trains (columns 1&2, 3&4): post-load wash, wash 1, wash 2, elute,
  strip, CIP, re-equilibrate; start-up sanitisation and benzyl alcohol storage
- Depth filtration: new filter every 3 eluates; flush, equilibrate, filter, chase,
  post-use flushes (solution order assumed)
- Tris base + NaCl in-line, AEX (2 columns); acetic acid in-line, CEX (2 columns)
- Virus filtration, 2 holders, alternating at each filter change
- UF/DF per lot (lot_window_h): UF, DF, UF, recovery, formulation spike, CIP, storage
- BDS hold mixer, 0.2 µm fill filter, cryovault containers
- dsp.tail_h: the downstream keeps running after the N run to process the rest

Each run returns `solutions` (litres per tote, WFI included), `chemicals` (kg),
`vessels` (peak volumes), `lots`, and `steps` (every step with start/end hour).
The diagram shows the current step on each unit; the start-up film shows the
steps running now; the step sequence film plays them as a timeline and records
MP4/WebM. Volumes and times are placeholders: replace them with the process
description. Options kept: vi.mode = "low_ph", polish.order = "cex_aex",
capture.mode = "pcc" (Protein A steps lumped per cycle in PCC mode).

Files added: `icb_dsp.py`, `static/seqfilm.js`.

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
   The `static`, `templates` and `data` folders must be included. If folders don't upload,
   use Chrome.
5. Click **Commit changes**.

Check the repository shows `app.py`, `render.yaml`, `requirements.txt`,
`templates/index.html`, `static/chart.umd.js`, `static/pfd.js`, `static/pfdvideo.js`,
`static/seqfilm.js`, `icb_dsp.py` and `data/icb/`
with four files. `index.html` must be in `templates/` (not `static/`): that is
where `app.py` loads it from.

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

- `icb_process.py` — parameters, titre curves, solutions (with compositions), waste streams,
  design check, portfolio sizing
- `icb_dsp.py` — step-by-step downstream (Case 5 flowsheet): every unit as a recipe of steps
- `icb_dynamic.py` — hour-by-hour simulation: N-1, N + TFF, surge tank, downstream, run-out;
  also records per-stream values and the step log for the diagram and films
- `icb_nevai.py` — NevAI titre and wastewater forecasters (hourly, 24 h ahead)
- `nevai_engine.py` — NevAI LSTM core (Python 3.9+ compatible)
- `icb_routes.py` — API routes; background training
- `app.py` — web app entry point; login page, sessions and logout
- `templates/index.html` — dashboard, print layout and PDF report
- `static/pfd.js` — live process flow diagram (no external libraries)
- `static/pfdvideo.js` — start-up film with the steps running now (MP4 / WebM)
- `static/seqfilm.js` — step-sequence film: one lane per unit, coloured by solution
- `static/chart.umd.js` — Chart.js 4.4.1 (MIT), bundled so no CDN is needed
- `data/icb/` — trained NevAI weights
- `render.yaml`, `requirements.txt`, `.python-version` — Render deployment
