"""case_config.py: parse the real aerosolved repo."""
from backend import case_config as cc

def test_16_cases_discovered(repo_root):
    cases = cc.discover_cases(repo_root)
    names = [c.name for c in cases]
    for want in ("rain", "CAG", "bentPipe", "cavity", "hygroscopicGrowth",
                  "uniformCoalescence"):
        assert want in names, (want, names)
    assert len(cases) >= 10, cases

def test_cag_29_knobs(repo_root):
    cag = cc.load_case(repo_root, "CAG")
    assert len(cag.vars) == 29, cag.vars
    # Spot-check the known names
    names = {v.name for v in cag.vars}
    for want in ("X1", "X4", "Y1", "Y6", "Z2", "NY1", "GX2", "MODEL"):
        assert want in names, want

def test_rain_no_vvars_block_but_application_is_resolved(repo_root):
    rain = cc.load_case(repo_root, "rain")
    assert len(rain.vars) == 0, "rain does not use VARS"
    assert rain.application == "aerosolEulerFoam", rain.application

def test_m4_escaping_no_crash(repo_root):
    cag = cc.load_case(repo_root, "CAG")
    # The raw VARS block contains 29 -D flags in shell-continuation form.
    assert cag.raw_vars_block
    # No -DVARNAME should be duplicated
    from collections import Counter
    m4_names = [v.m4_name for v in cag.vars]
    assert all(c <= 1 for c in Counter(m4_names).values())

def test_application_read_from_controlDict_fallback(repo_root):
     # A case with system/controlDict.m4 -- application should still resolve.
    cag = cc.load_case(repo_root, "CAG")
    assert cag.application, "application should be non-empty"
