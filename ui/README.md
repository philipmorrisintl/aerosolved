# AeroSolved Interactive GUI

A desktop front-end for the AeroSolved OpenFOAM library, designed to run on a
Linux machine with an installed OpenFOAM-2412 (or newer). The UI is built on
**PySide6** (the official Python binding for Qt 6) and wraps the existing
`Allrun` → `m4` → `wmake`/`run` workflow without modifying the C++ solver code.

## What the UI does

| Capability | Detail |
|------------|--------|
| **Case discovery** | Automatically enumerates every `cases/*` that has an `Allrun` — no manual config needed. |
| **Knob editing** | Every `-DVAR...=$SHELL_VAR` pair the case's `Allrun` feeds to `m4` becomes a text box. Defaults are pre-filled from the Allrun's shell assignments. |
| **Model selection** | The positional args (`mesh`, `sectional`/`moment`) are exposed as dropdowns with the case's default values. |
| **Run** | One click runs the (patched) `Allrun` in an isolated run directory, streams the full log with phase highlighting (mesh → fields → solve → postProcess). |
| **Dashboard** | After the run, the UI parses `postProcessing/` into a dashboard: (a) a plain-text summary of every plot + scalar; (b) PDF/PNG figures via the matplotlib `Agg` backend; (c) a JSON dashboard spec for scripting. |

## Install

From an aerosolved checkout with the `ui/` package on your branch (the one
you're looking at now):

```bash
# System Qt / Python — one of these
sudo apt-get install python3-pyqt6          # for PyQt6, OR
pip install --user pySide6                  # for PySide6, OR
conda install -c conda-forge pyside6

# Python deps (the solver's own deps live in the repo's requirements.txt)
pip install --user -r ui/requirements-ui.txt
```

## Run

```bash
# The launcher
python -m ui
# or, explicitly
python ui/frontend/app.py
```

The main window opens with every available case in a drop-down. Pick a case,
adjust the knob values in the form, hit **Run** — the run directory is
created, the log streams, and after completion the **Results** tab is filled
with the dashboard.

## File layout

```
ui/
    __init__.py              # package root, re-exports `ui.frontend.app`
    README.md                # this file
    requirements-ui.txt      # pip deps for the UI (PySide6, numpy, matplotlib)
    .gitignore

    backend/                 # pure Python — no Qt, fully unit-testable
    __init__.py
    case_config.py           # parse Allrun -> VARS shell vars, model choices, solver
    m4_expand.py             # mirror of setMacros(); runs m4 with -D flags
    runner.py                # RunManager + RunEvent stream + patch_allrun_knobs
    results.py               # per-shape extractors + generic scanner -> dashboard spec
    dashboard.py             # render_text / render_matplotlib / spec_json

    frontend/                # Qt layer (guarded import, importable headless)
    __init__.py
    app.py                   # MainWindow, KnobStore, RunWorker, DashboardWorker

    tests/                   # pytest suite (runs headless, no display)
    conftest.py
    test_case_config.py
    test_runner.py
    test_results.py
    test_dashboard.py
    test_frontend.py
```

## Design notes

- **No C++ changes.** The UI is a pure-orchestration layer around the existing
   `Allrun` shell scripts. It *patches the Allrun copy in the run directory* to
   override the shell-variable defaults, then runs it — exactly as if the user
   had edited the script.
- **Headless-testable.** `ui.backend` has zero Qt deps and can be imported and
   tested without a display. `ui.frontend.app` guards Qt with a try/except so
   non-graphical environments (CI, SSH without X, remote desktops) can still
   import the package and test the orchestration logic.
- **Case-agnostic.** `case_config.py` parses **every** case in the repo, not
   just `rain`/`CAG`. New cases that follow the AeroSolved Allrun VARS
   conventions work automatically — no UI update needed.
- **Dashboard spec is JSON**. The `DashboardSpec` dict is JSON-serialisable, so
   the same data structure drives the Qt figure viewer, a future browser
   dashboard, and the CLI reporter.
- **Idempotent run directory.** Re-running the same case with the same knobs
   reuses the existing run directory and cleans transient artefacts (no need
   to `rm -rf`); a knob change causes a fresh copy.

## Limitations

- The UI does **not** auto-generate `0/` or `Time/` fields — it relies on the
   case's `0.org` → `0` copy which `Allrun` already does.
- The "dashboard" is deliberately *not* a 3D visualisation — it shows the
   1-D sectional flux distributions and mass budgets that the existing
   `plot.py` renders. A 3D field viewer would require additional parsing of
   the `Time/<t>/W` or `U` OpenFOAM binary fields, which is out of scope here.
- `matplotlib`'s `Agg` backend is used for the figures (headless-safe); no
   Qt backend is needed so the backend can be exercised in CI without a
   display.
- The Qt binding is auto-detected (PySide6 → PyQt6 fallback → None). If neither
   is installed the UI still starts its Qt-free backend and reports a clean
   error.
