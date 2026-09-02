"""frontend/app.py: MainWindow, KnobStore, DashboardWorker (mock Qt)."""
import os
from frontend.app import create_window, KnobStore, DashboardWorker
from backend import case_config as cc

def test_main_window_loads_all_cases(repo_root, mock_qt):
    win = create_window(mock_qt, repo_root)
    assert len(win.case_list) >= 10
        # All cases should be available as items in the picker
    assert hasattr(win, "cases_picker")
    assert hasattr(win, "log_view")
    assert hasattr(win, "results_view")

def test_knob_store_roundtrip(repo_root):
    cag = cc.load_case(repo_root, "CAG")
    ks = KnobStore(cag)
     # 29 knob values, seeded from defaults
    assert len(ks.values) == 29, ks.values
     # CAG's positional arg has "mesh" (its Allrun reads $1).
    assert "mesh" in ks.args, ks.args

def test_knob_editing_roundtrip(repo_root):
    cag = cc.load_case(repo_root, "CAG")
    ks = KnobStore(cag)
    ks.set("X1", "-3")
    ks.set("mesh", "100")
    assert ks.get("X1") == "-3"
    assert ks.get("mesh") == "100"
        # to_vars_dict maps "mesh" -> "MESH" for the Allrun
    v = ks.to_vars_dict()
    assert v["X1"] == "-3"
    assert v["MESH"] == "100"

def test_dashboard_worker_headless(fake_run_dir):
    fig_dir = os.path.join(fake_run_dir, "figs")
    dw = DashboardWorker("rain", fake_run_dir, fig_dir)
    out = dw.run_headless()
    assert out["ok"]
    assert out["plots"] >= 2
        # JSON spec was written
    assert os.path.isfile(os.path.join(fig_dir, "dashboard.json"))
