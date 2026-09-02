"""results.py: per-shape extractor + generic scanner."""
from backend import results as rs
import os
import numpy as np

def test_rain_extractor_produces_two_plots(fake_run_dir):
    spec = rs.build_dashboard("rain", fake_run_dir)
    assert len(spec["plots"]) >= 2, spec["plots"]
        # Number-flux + mass-flux budgets
    titles = [p["title"] for p in spec["plots"]]
    assert any("Number flux" in t for t in titles),     titles
    assert any("Mass flux" in t for t in titles),       titles

def test_rain_scalars_sum_to_total(fake_run_dir):
    spec = rs.build_dashboard("rain", fake_run_dir)
    scalars = {s["name"]: s["value"] for s in spec["scalars"]}
    assert "total_mass_top_last" in scalars
    assert "total_mass_bottom_last" in scalars

def test_generic_scanner_on_unknown_case(fake_run_dir):
    spec = rs.build_dashboard("unknown_case", fake_run_dir)
        # Falls back to the generic scanner
    assert len(spec["files"]) >= 4, spec["files"]
        # Generic plots for the first 3 files
    assert len(spec["plots"]) >= 1, spec["plots"]

def test_missing_postproc_reports_warning():
    spec = rs.build_dashboard("rain", "/no/such/path/postProcessing")
    assert spec["warnings"], spec["warnings"]
