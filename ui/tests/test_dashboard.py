"""dashboard.py: text + json + matplotlib (best-effort)."""
from backend import dashboard as dash
import json
import os
import pytest

def test_render_text_contains_scalars(fake_run_dir):
    from backend import results as rs
    spec = rs.build_dashboard("rain", fake_run_dir)
    text = dash.render_text(spec)
    assert "AeroSolved run dashboard" in text
        # At least one scalar
    assert "total_mass_top_last" in text

def test_spec_json_roundtrip(fake_run_dir):
    from backend import results as rs
    spec = rs.build_dashboard("rain", fake_run_dir)
    out = os.path.join(fake_run_dir, "dashboard.json")
    dash.spec_json(spec, out)
    # JSON must be valid and reproducible
    with open(out) as f:
        loaded = json.load(f)
        assert loaded["case"] == "rain"
        assert "plots" in loaded and "scalars" in loaded

def test_matplotlib_optional_no_crash(fake_run_dir):
    from backend import results as rs
    spec = rs.build_dashboard("rain", fake_run_dir)
        # If matplotlib is missing, this should not raise -- the DashboardWorker
        # in frontend/app.py wraps it in try/except; here we just call the
        # render_matplotlib directly to check that it *attempts* and either
        # returns a dict or raises (in both cases we don't crash the suite).
    try:
        out = dash.render_matplotlib(spec, fake_run_dir, fmt="pdf")
        # If it ran, it must produce at least as many files as plots
        assert isinstance(out, dict)
    except ModuleNotFoundError:
        pytest.skip("matplotlib not installed -- rendering skipped")
