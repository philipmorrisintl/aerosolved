"""frontend/app.py -- PySide6 desktop GUI for AeroSolved (guarded Qt import).

Layout (single QMainWindow):
   | CasePicker | KnobEditor (form) | Run + phase bar |
   +----------------------------------------------------+
   | LiveLog (QPlainTextEdit)        | ResultsPanel      |
   |                                  |  (dashboard text  |
   |                                  |   + renderers)    |

The Qt objects are only created when a display is available. The run
orchestration (RunWorker / DashboardWorker) is plain, Qt-independent logic so
it can be exercised in a headless test.
"""
from __future__ import annotations

import os
import sys
import logging
from typing import Optional, Callable

_PKG = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _PKG not in sys.path:
    sys.path.insert(0, _PKG)

from backend import case_config as cc
from backend import runner as rn
from backend import results as rs
from backend import dashboard as dash

log = logging.getLogger(__name__)


# --- Guarded Qt import -------------------------------------------------------

QT = None


def load_qt():
    """Import a Qt binding (PySide6 preferred, PyQt6 fallback, else None)."""
    global QT
    if QT is not None:
        return QT
    try:
        from PySide6 import QtCore, QtWidgets, QtGui
        QT = type("QtNS", (), {"QtCore": QtCore, "QtWidgets": QtWidgets,
                               "QtGui": QtGui, "binding": "PySide6"})()
        return QT
    except Exception:
        pass
    try:
        from PyQt6 import QtCore, QtWidgets, QtGui
        QT = type("QtNS", (), {"QtCore": QtCore, "QtWidgets": QtWidgets,
                               "QtGui": QtGui, "binding": "PyQt6"})()
        return QT
    except Exception:
        return None


qt = load_qt()


# --- Knob model (Qt-independent) -------------------------------------------

class KnobStore:
    """Holds the current value of every knob for a case, seeded from defaults."""

    def __init__(self, case_def):
        self.case = case_def
        self.values = {v.name: (v.default or "") for v in case_def.vars}
        self.args = self._seed_args()

    def _seed_args(self):
        args = {}
        for a in self.case.positional_args:
            if a == "mesh":
                args[a] = "20"
            elif a == "model":
                args[a] = (self.case.model_choices[0]
                           if self.case.model_choices else "sectional")
            else:
                args[a] = ""
        return args

    def get(self, key):
        return self.values.get(key, self.args.get(key, ""))

    def set(self, key, value):
        if key in self.values:
            self.values[key] = value
        elif key in self.args:
            self.args[key] = value

    def to_vars_dict(self):
        out = dict(self.values)
        if "mesh" in self.args and "MESH" not in out:
            out["MESH"] = self.args["mesh"]
        return out

    def run_args(self):
        mesh = self.args.get("mesh")
        model = self.args.get("model", "sectional")
        return (mesh or None, model)

    def model_choices(self):
        return list(self.case.model_choices or ["sectional", "moment"])


# --- Run orchestration --------------------------------------------------------

class RunWorker:
    """Drive a RunManager, surfacing its events. Headless generator OR Qt."""

    def __init__(self, runner):
        self.runner = runner
        self.events = []
        self._emit = None
        self.done = False
        self.return_code = None
        self.final_phase = ""

    def run_headless(self):
        for ev in self.runner.stream():
            d = ev.to_dict()
            self.events.append(d)
            if d["type"] == "done":
                self.done = True
                self.return_code = d.get("detail")
                self.final_phase = d.get("phase", "")
        return {"events": len(self.events), "phase": self.final_phase,
                "code": self.return_code, "results_dir": self.runner.results_dir()}

    def install_emitter(self, emit):
        self._emit = emit

    def run_threaded(self):
        for ev in self.runner.stream():
            d = ev.to_dict()
            self.events.append(d)
            if self._emit is not None:
                self._emit(d)
            if d["type"] == "done":
                self.done = True
                self.return_code = d.get("detail")
                self.final_phase = d.get("phase", "")


def build_run_worker(case_def, run_dir, mesh, model, knobs, timeout=None):
    cfg = rn.RunConfig(case_path=case_def.path, run_dir=run_dir,
                        model=model, mesh_cells=mesh, knobs=knobs,
                        timeout=timeout)
    mgr = rn.RunManager(cfg)
    return RunWorker(mgr), mgr


# --- Dashboard assembly (Qt-independent) ------------------------------------

class DashboardWorker:
    """Parse a finished run's postProcessing output into a dashboard spec."""

    def __init__(self, case_name, run_dir, fig_dir):
        self.case_name = case_name
        self.run_dir = run_dir
        self.fig_dir = fig_dir
        self.spec = None
        self.figures = []
        self.text = ""
        self.done = False
        self.error = ""
        self._emit = None

    def run_headless(self):
        os.makedirs(self.fig_dir, exist_ok=True)
        pp = rs.latest_postproc(self.run_dir)
        if pp is None:
            self.error = "no postProcessing output found for the run"
            self.done = True
            return {"ok": False, "error": self.error}
        self.spec = rs.build_dashboard(self.case_name, pp)
        self.text = dash.render_text(self.spec)          # always works
        dash.spec_json(self.spec, os.path.join(self.fig_dir, "dashboard.json"))
        self.figures = []
         # Matplotlib figures are best-effort.  If the matplotlib backend is
         # unavailable (e.g. a headless environment without the dep) the
         # dashboard is still delivered, just without rendered figures.
        try:
            self.figures = dash.render_matplotlib(self.spec, self.fig_dir,
                                                    fmt="pdf")
        except Exception as e:
            log.warning("matplotlib render failed -- figures skipped: %s", e)
            self.figures = []
        self.done = True
        return {"ok": True, "plots": len(self.spec.get("plots", [])),
                 "figures": len(self.figures),
                 "scalars": len(self.spec.get("scalars", []))}

    def install_emitter(self, emit):
        self._emit = emit

    def run_threaded(self):
        try:
            self.run_headless()
        except Exception as e:
            self.error = str(e)
            self.done = True
        if self._emit is not None:
            self._emit({"type": "done", "ok": not self.error})


# --- UI assembly (Qt) ---------------------------------------------------------

class MainWindow:
    """Top-level window built against a Qt namespace `QT`."""

    def __init__(self, QT, repo_root):
        self.QT = QT
        self.repo_root = repo_root
        self.case_list = cc.discover_cases(repo_root)
        self.current = None
        self.knobs = None
        self.worker = None
        self._build()

    def _build(self):
        self.cases_picker = self._make_case_picker()
        self.run_button = self._make_run_button()
        self.log_view = self._make_log_view()
        self.results_view = self._make_results_view()
        self.phase_label = self._make_phase_label()

    def on_case_selected(self, name):
        for c in self.case_list:
            if c.name == name:
                self.current = c
                self.knobs = KnobStore(c)
                break
        return c

    def on_run(self, run_dir, timeout=None):
        if self.knobs is None:
            raise RuntimeError("select a case before running")
        mesh, model = self.knobs.run_args()
        worker, mgr = build_run_worker(self.current, run_dir, mesh, model,
                                        self.knobs.to_vars_dict(), timeout)
        self.worker = worker
        return worker.run_headless()

    def on_run_complete(self, run_dir, fig_dir):
        dw = DashboardWorker(self.current.name, run_dir, fig_dir)
        return dw.run_headless()

    def _make_case_picker(self):
        w = self.QT.QtWidgets.QComboBox()
        for c in self.case_list:
            w.addItem("%s (%d knobs)" % (c.name, len(c.vars)))
        return w

    def _make_run_button(self):
        return self.QT.QtWidgets.QPushButton("Run")

    def _make_log_view(self):
        return self.QT.QtWidgets.QPlainTextEdit()

    def _make_results_view(self):
        return self.QT.QtWidgets.QPlainTextEdit()

    def _make_phase_label(self):
        return self.QT.QtWidgets.QLabel("idle")


def create_window(qt, repo_root):
    """Factory used by main() -- separated so it can be called with a mock."""
    return MainWindow(qt, repo_root)


def main(argv=None):
    qtns = load_qt()
    if qtns is None:
        sys.stderr.write("no Qt binding found (pip install pySide6 or PyQt6) "
                         "-- cannot start GUI\n")
        return 2
    repo = cc.find_repo_root(os.getcwd())
    if repo is None:
        sys.stderr.write("could not locate an aerosolved checkout; pass one.\n")
        return 2
    app = qtns.QtWidgets.QApplication.instance()
    if app is None:
        app = qtns.QtWidgets.QApplication([argv[0] if argv else sys.argv[0]])
    win = create_window(qtns, repo)
    win.show()
    return app.exec()


if __name__ == "__main__":
    raise SystemExit(main())
