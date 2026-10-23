# Integrated Dynamic BioProcessing (IDB) — mAb digital twin

N-1 perfusion seed (media in = cell-free permeate out to waste) → N reactor
(feed in = cells + broth + product out) → TFF cell removal (cell concentrate to
bio-kill) → conditioned media surge tank → Brij-9 in-line → Protein A (columns
1&2, 3&4) → depth filtration → Tris base/NaCl in-line → AEX → acetic acid in-line
→ CEX → virus filtration (2 holders) → UF/DF → formulation spike → BDS hold →
0.2 µm → cryovault. Every downstream step runs in sequence.

## Basis-of-design cases

| # | Case | Avg titre (g/L/day) | Batches/yr | Yield | Capture load (g/L) |
|---|------|------|----|-----|----|
| 1 | Design Case: High Yield/Volume | 3.0 | 20 | 71 % | 80 |
| 2 | Design Case: Turndown Volume | 3.0 | 3 | – | – |
| 3 | Design Case: Low Yield/Volume | 0.4 | 3 | 25 % | 35 |
| 4 | Rating Case: Rilveg | 0.6 | 20 | – | – |
| 5 | Scenario Case: High-High Yield | 4.8 | 20 | 71 % | – |

Cases 1–3 bracket the operating range; the **design case** is their envelope (each
item sized for the case that needs the most). Case 4 rates Rilveg for cost of goods.
Case 5 is simulated and compared with Case 1 for the incremental-capital decision
(Note 1).

How a case is applied (`icb_cases.py`, on top of the Inputs, at run time):

- **Titre**: the chosen titre-curve shape is scaled (`upstream.titre_scale`) so its
  run average equals the case value. Default basis: g made per L of N working volume
  per day (= mean permeate titre × VVD); the alternative basis is mean permeate titre
  in g/L. Shape "auto" uses the Figure 2 curve needing the least scaling.
- **Yield**: harvest-to-DS train yield (Protein A → fill); every step yield is scaled
  by the same factor.
- **Capture load**: g per L resin; `DBC_gL = load / load_util`, column volume kept.
- **Batches/yr**: per-year totals = one batch × batches; reactor and downstream lines
  from days per batch + turnaround (default 3 d) against operating days (default 330).
- **Blank cells**: Cases 2 and 5 inherit yield and load from Case 1; Case 4 keeps the
  Inputs values. All values and assumptions can be edited on the Cases tab.

## Dashboard tabs

- **Case bar**: Case 1 … Case 5, Design case (envelope of Cases 1–3), No case.
- **Cases**: editable case table, how cases are applied, noise-free preview of every
  case, and "Run all 5 cases": side-by-side results, per-year totals, line counts,
  the design envelope sizing (with "Set by") and the Case 5 vs Case 1 capacity table
  (ratio and six-tenths scaling factor per item).
- **Inputs**: run set-up and all parameters. Orange-edged fields are set by the active
  case at run time; the table of what the case changed is shown here.
- **HMB**: product (mAb) balance with closure, water balance, totals per stream
  (streams 1–63) per batch and per year, and the stream table at one hour.
- **Design**: design check and daily plan, equipment sizing (apply and re-run), sizing
  across titre curves.
- **Outputs**: case checks, per-batch and per-year table, occupancy, KPIs, charts,
  solutions, raw materials, vessels, lots, wastewater, losses, events.
- **Schematics & diagrams**: live PFD, start-up film, step-sequence film.
- **Print / Save as PDF** prints every tab except Inputs, plus all parameters.

Default parameters are illustrative, and the titre curves were read off a chart
image. Replace both with site data before relying on absolute numbers.

## Run on a Mac

    python3 -m venv .venv
    source .venv/bin/activate
    pip install -r requirements.txt
    python3 icb_dynamic.py          # quick check: prints design, KPIs and solutions
    python3 app.py                  # open http://127.0.0.1:8000

Optional: `python3 icb_nevai.py train` retrains NevAI (8 runs × 30 epochs, about
2 minutes) and overwrites the weights in `data/icb/`. NevAI is trained on the base
curves; forecasts for strongly scaled cases (e.g. Case 5) are extrapolations.
Training also calibrates split-conformal 10–90 % forecast bands on two held-out runs.

## NevAI in operation

NevAI is advisory by default. Set `control.nevai_act` to let it start an early partial
Protein A load when the forecast says the surge-tank product would pass
`max_sump_hold_h` (minus `nevai_hold_margin_h`) before a full load is ready; the
smallest partial load is `nevai_partial_min_frac` of a full one.

The watch compares each measured permeate titre with the forecast issued
`nevai_watch_lead_h` earlier and raises the `nevai_watch` alarm when the forecast
error shifts by more than `nevai_watch_threshold`. Calibrate the threshold from a run
without disturbances: set it to about 1.5 × the reported `nevai_watch_peak_score`.
The flag also shows in the PFD readout and in the start-up film.

The NevAI card on the Outputs tab compares buffer preparation driven by the 24 h
forecast with a reactive rule (prepare what was used in the last 24 h), at the same
preparation capacity (`buffer_prep_Lh`, 0 = peak day ÷ 20 h). **Measure NevAI benefit**
runs the selected case advisory and acting with the same seed and tables the difference.
Retrain once after this update to get calibrated bands; existing weights keep working.

Port 8000 is used because macOS reserves 5000 for AirPlay.

## Using your own titre data

- Quick test: pick **Custom curve** (Inputs tab, or as the case shape on the Cases
  tab), then paste `[[N day, g/L], ...]` into the box on the Inputs tab.
- Permanent: replace the points in `TITRE_PROFILES` in `icb_process.py`, then retrain NevAI.

## Put it on GitHub

Create an empty repository, then drag-and-drop the unzipped folder (show hidden files
with Cmd + Shift + .) or:

    cd icb-mab-twin
    git init
    git add .
    git commit -m "IDB dynamic twin with Cases 1-5"
    git branch -M main
    git remote add origin https://github.com/YOUR-USER/icb-mab-twin.git
    git push -u origin main

Check the repository shows `app.py`, `icb_cases.py`, `icb_routes.py`, `icb_dsp.py`, `icb_nevai_ops.py`,
`icb_sizing.py`, `render.yaml`, `requirements.txt`, `templates/index.html`,
`static/chart.umd.js`, `static/pfd.js`, `static/pfdvideo.js`, `static/seqfilm.js`
and `data/icb/` with four files.

## Deploy on Render

1. render.com → **New → Blueprint** → pick the repository (Render reads `render.yaml`).
2. Set `ADMIN_EMAIL` and `ADMIN_PW`; `SECRET_KEY` is generated.
3. Deploy, open the service URL and sign in.

"Run all 5 cases" runs five simulations in one request: if the start command uses
gunicorn's default 30 s timeout, use e.g. `gunicorn app:app --timeout 180`.
The free instance sleeps when idle and its disk resets on restart (weights committed
in `data/icb/` are always there).

## API

    GET  /api/icb/defaults        default parameters, design check, titre curves
    GET  /api/icb/cases           Cases 1-5 and default case settings
    POST /api/icb/cases/preview   noise-free plan of every case {"params", "case": {"settings", "values"}}
    POST /api/icb/cases/compare   simulate all cases; design envelope; Case 5 vs Case 1
    POST /api/icb/design          design check
    POST /api/icb/portfolio       buffer / sump sizing across titre curves
    POST /api/icb/sizing          equipment sizing, worst case across titre curves
    POST /api/icb/simulate        {"params", "seed", "downsample_h", "use_forecaster", "case": {"id", "settings", "values"}}
    POST /api/icb/forecast        run, then NevAI 24 h bands
    GET  /api/icb/nevai/status    model and training state
    POST /api/icb/nevai/train     start background training
    POST /api/icb/nevai/ingest    hourly plant rows -> fine-tune
    POST /api/icb/nevai/benefit   same case advisory vs acting (control.nevai_act), same seed

A run with a case returns `case.info` (what the case changed) and `case.annual`
(per-year totals, occupancy, checks against the case values).

## Files

- `icb_cases.py` — Cases 1–5: apply a case, per-year totals, occupancy, preview, all-case
  comparison, design envelope, Case 5 vs Case 1
- `icb_process.py` — parameters, titre curves, solutions, waste streams, design check, portfolio sizing
- `icb_dsp.py` — step-by-step downstream: every unit as a recipe of steps
- `icb_sizing.py` — equipment sizing from the hourly run
- `icb_dynamic.py` — hour-by-hour simulation
- `icb_nevai.py`, `nevai_engine.py` — NevAI forecasters and calibrated bands
- `icb_nevai_ops.py` — NevAI in operation: hold guard, drift watch, buffer-preparation comparison
- `icb_routes.py` — API routes; background training
- `app.py` — web app entry point; login and sessions
- `templates/index.html` — tabbed dashboard and PDF report
- `static/pfd.js`, `static/pfdvideo.js`, `static/seqfilm.js`, `static/chart.umd.js`
- `data/icb/` — trained NevAI weights
- `render.yaml`, `requirements.txt`, `.python-version` — Render deployment
