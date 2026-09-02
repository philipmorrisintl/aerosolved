"""runner.py: patch_allrun_knobs + RunManager stream."""
from backend import runner as rn
import os
import tempfile

def test_patch_allrun_knobs_selective():
    src = ("NY1=$MESH\nNX1=5\nY3=1.5\n# comment line\n"
              "case $2 in\n  sectional) MODEL=fixedSectional ;;\n  moment) MODEL=twomoment;;\nesac\n")
    patched = rn.patch_allrun_knobs(src, {"NY1": "32", "NX1": "64", "Y3": "2.0"})
    assert "NY1=32" in patched
    assert "NX1=64" in patched
    assert "Y3=2.0" in patched
    assert "case $2 in" in patched, "over-matched the model block"
    assert "MODEL=fixedSectional" in patched,   "over-matched the model block"

def test_patch_allrun_knobs_noop_when_no_knobs():
    src = "X1=1\nX2=2\n"
    assert rn.patch_allrun_knobs(src, {}) == src

def test_runmanager_full_pipeline_stream(fake_case):
    run_dir = tempfile.mkdtemp() + "_run"
    cfg = rn.RunConfig(case_path=fake_case, run_dir=run_dir,
                          model="sectional", mesh_cells="20", knobs={"NY1":"32"})
    events = list(rn.RunManager(cfg).stream())
    done = [e for e in events if e.type == "done"][-1]
    assert done.phase == "success"
    # At least one log event should carry the "Courant Number" string.
    assert any("Courant Number" in e.message for e in events if e.type == "log")
    pp = rn.RunManager(rn.RunConfig(case_path=fake_case, run_dir=run_dir,
                                       model="sectional")).results_dir()
    assert pp is not None
    assert os.path.isdir(pp)

def test_runmanager_prepare_reuses_existing_run_dir(fake_case):
    # Second call to prepare() with the same run_dir must not raise.
    run_dir = tempfile.mkdtemp() + "_reused"
    cfg = rn.RunConfig(case_path=fake_case, run_dir=run_dir,
                          model="sectional", mesh_cells="20", knobs={"NY1":"32"})
    rn.RunManager(cfg).prepare()
    # After the first prepare, Allrun is present
    assert os.path.isfile(os.path.join(run_dir, "Allrun"))
    # A second prepare with the same run_dir should NOT crash
    rn.RunManager(rn.RunConfig(case_path=fake_case, run_dir=run_dir,
                                  model="sectional")).prepare()
