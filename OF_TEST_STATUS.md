# OpenFOAM Build & Smoke-Test Status Report
**Date**: 2026-08-30  (first-run diagnosis appended §5)
**Repository**: `philipmorrisintl/aerosolved`
**Local checkout**: `/Users/aspen_mb/Documents/VS_CODE_PROJECTS/Aerosolved/aerosolved`
**Changeset tarball**: `aerosolved.changes.tar.gz` (37 files, 40 KB; re-extract to pick up the v2 harness)

---

## 1. What was done locally (all verified on this machine)

### Bug fixes — verified
| # | File | Fix | Verification method | Status |
|---|------|-----|--------------------|--------|
| 1 | `scripts/CAG.py:72` | `massFlowAirHot` duplicate → `massFlowAirCold` | Computed old/new totals; old=13/60, new=14/60 | ✅ DONE |
| 2 | `scripts/Cloud.py:38` | `"mixutre"` → `"mixture"` | grep confirm | ✅ DONE |
| 3 | `doc/Chap5_Nomenclature.md:29` | `"object3"` → `"object"` | grep confirm | ✅ DONE |
| 4 | `doc/Chap1_About.md` | `OpenFOAM-v1812` → `v2406/v2412` + link to new Installation chapter | grep confirm | ✅ DONE |
| 5 | `cases/bentPipe/Allrun:116` + `cases/Brownian/Allrun:28` | Add missing trailing backslash in `VARS` block | Empirically tested via `m4` — see §2 | ✅ DONE (cosmetic) |
| 6 | `scripts/thermophysicalFunctions.py:69` | air molar mass `28.810` → `28.9596` | grep confirm | ✅ DONE |
| 7 | 23 `*.py` shebangs | `#!/usr/bin/python` → `#!/usr/bin/env python3` | count check: 0 remaining | ✅ DONE |
| 8 | `scripts/AeroSolvedRunFunctions` | `set -o pipefail` + dep probe; `bash -n` OK | `bash -n` | ✅ DONE |

### Environment fixes — verified
- **`requirements.txt`** created (`numpy>=1.25`, `scipy>=1.9`, `matplotlib>=3.5`).
- Project venv at `Aerosolved/.venv-aerosolved/` (Python 3.14 + numpy 2.5.2 + scipy + matplotlib).
- All self-contained scripts **run end-to-end**: `CAG.py`, `Cloud.py`, `logNormalSectionalGrid.py` (PDF produced), `blendedCoalescence.py` (PDF produced).
- All `*.py` byte-compile clean.
- All `Allrun` + helper scripts pass `bash -n`.
- CI `ci.yml` passes YAML parse (`jobs: lint / python-scripts / build`).

### New deliverables — added
| File | Purpose |
|------|---------|
| `doc/Chap0_Installation.md` | Full installation/build/run chapter; numpy-ABI troubleshooting baked in |
| `.github/workflows/ci.yml` | 3-job CI: lint / python-scripts / build matrix (OpenFOAM 2406+2412) |
| `build_and_smoketest.sh` | One-shot self-diagnosing build+smoke-test harness (v2 — auto-detects OpenFOAM, falls back to `./Allwmake`, captures full logs) |
| `requirements.txt` | Python dependency manifest |
| `IMPROVEMENTS.txt` | 28 tiered recommendations with `file:line` |
| `OF_TEST_STATUS.md` (this file) | Handoff report + verification matrix + first-run diagnosis |

---

## 2. Critical correction: the `Allrun` backslash fix

The original IMPROVEMENTS classified the missing `\` on `-DVARMWALLBC=$MWALLBC`
(bentPipe) and `-DVARMBC=$MBC` (Brownian) as a **"critical silent
misconfiguration."** Empirical testing **disproved this**:

```
[fixed] $VARS = -DVARMINLETBC=... -DVARMWALLBC=... -DVARZWALLBC=...
[buggy] $VARS = -DVARMINLETBC=... -DVARMWALLBC=... -DVARZWALLBC=...
```
Because `VARS="..."` uses **double quotes**, a missing backslash leaves a
*liter*al newline in the string, which `m4` still parses as whitespace — so
**all three macros reach `m4` in both cases**. The fix is now documented as
`[INFO/LOW]` (cosmetic/style), not a functional bug. The only files that
actually change runtime behaviour are items 1, 2, 6, 7. `IMPROVEMENTS.txt`
has been corrected accordingly.

---

## 3. What still requires the OpenFOAM machine

| Item | Change | Verification needed |
|------|--------|---------------------|
| `scripts/AeroSolvedRunFunctions` | `set -o pipefail` added | `wmake` + `Allrun` runs without failing on unset vars |
| `cases/*/Allrun` backslash fixes | Cosmetic, but confirm smoke still passes | `build_and_smoketest.sh` |
| `generateGitInfo.sh` | Portable rewrite (no GNU `-mmin`/`touch -d`) | `make compile` on GNU *and* BSD |
| `requirements.txt` | New Python deps | `pip install -r requirements.txt` succeeds |
| `.github/workflows/ci.yml` | New 3-job CI | Push to GitHub, confirm lint + build jobs green |
| `*New.C` (13 files) | Left in place — need `wmake` to confirm removal is safe | Deferred; requires build confirmation |
| `Continuum coalescence` sub-model | Not yet implemented | Deferred |
| `Allrun` mesh math → `awk` | Not changed | Deferred |

---

## 4. How to test on the OpenFOAM machine

**Where the tarball lands: on your aerosolved *source checkout* (the thing you
`git clone`d) — NOT on the OpenFOAM install dir, and NOT into a fresh/empty
directory. It is an overlay of the changed files onto the existing tree; the
build needs the unchanged `Makefile`, `Allwmake`, `cases/*` etc. that stay.**

The simplest path is the bundled `bootstrap_openfoam.sh`, which finds the
aerosolved root, overlays the tarball, activates OpenFOAM by trying every
common method, and runs the test. Run it with **no OpenFOAM path guessed** —
it discovers the activation itself:

       bash bootstrap_openfoam.sh                    # auto-detect
       bash bootstrap_openfoam.sh /path/to/aerosolved <tarball>
        (if auto-detect misses)

If it can't activate OpenFOAM it prints discovery commands — run **one** of
these and paste the output back so the exact line can be pinned for THIS box
(the earlier `bin/foam` under `WM_PROJECT_DIR` was wrong for a Debian/Ubuntu
`openfoam2412` install, which puts `WM_PROJECT_DIR` at
`/usr/lib/openfoam/openfoam2412` with *no `bin/foam` under it*):

       find / -name bashrc -path '*openfoam*' 2>/dev/null | head
       dpkg -L openfoam2412 2>/dev/null | grep -E 'bashrc|/bin/foam$' | head
       echo "WM_PROJECT_DIR=[$WM_PROJECT_DIR] FOAM_INST_DIR=[$FOAM_INST_DIR]"

**Manual equivalent** (once `bashrc` is known):

     cd <aerosolved-checkout>
     tar -xzf aerosolved.changes.tar.gz        # overlay the 37 files
     . <path-to-OpenFOAM>/etc/bashrc            # the real env script
     pip install -r requirements.txt
     bash build_and_smoketest.sh

### v2 harness features
The updated `build_and_smoketest.sh`:
- Detects OpenFOAM at multiple locations (distro, ETSI, `$WM_PROJECT_DIR`, `$FOAM_INST_DIR`).
- Auto-sources `bin/foam` if OpenFOAM isn't already activated.
- Falls back to `./Allwmake` when `make compile` fails or has no `compile` target.
- Captures full build and per-case logs to `/tmp/aerosolved_*.log`.
- Prints a structured `FINAL SUMMARY` at the end.

### Expected output on success
```
[0/6] environment diagnostics    ... OpenFOAM found at /usr/lib/openfoam/openfoam2412
[1/6] checking prerequisites     ... wmake present, m4 present
[2/6] building                   ... strategy: make compile
            ... wmake output ...
        build OK — checking for built executables
        -rwxr-xr-x  aerosolEulerFoam
[3/6] Python scripts             ... all python scripts OK
[4/6] VARS templates             ... Allrun VARS templates OK
[5/6] smoke case                 ... cases/rain/Allrun: OK
==================================================================
  FINAL SUMMARY
  Python scripts   : PASSED
  VARS templates   : PASSED
  Smoke cases      : PASSED
================================================================
ALL SMOKE TESTS PASSED
  build log: /tmp/aerosolved_build_YYYYMMDD_HHMMSS.log
```

### On failure
The script prints `FINAL SUMMARY` with which steps failed, then prints
`tail -n 30` of the build log. **Paste the last 30 lines of
`/tmp/aerosolved_build_YYYYMMDD_HHMMSS.log`** when reporting a failure.

---

## 5. Diagnosis of 2026-08-30 first run on aarch64 box

**Observed output**:
```
[1/6] checking environment
WARNING: git not found (commit info will be empty)
    OpenFOAM version: <unset>
    WM_PROJECT_DIR   : /usr/lib/openfoam/openfoam2412
[2/6] building (make compile)
make: *** No rule to make target 'compile'.  Stop.
real    0m0.003s
ERROR: build failed
```

### Root causes identified
1. **`WM_PROJECT_DIR` is set but `FOAM_RELEASE` is unset**: The box has a
   **Debian/Ubuntu distro-package OpenFOAM 2412** at `/usr/lib/openfoam/openfoam2412/`.
   Distro packages set `WM_PROJECT_DIR` but do NOT set `FOAM_RELEASE`.
   **Fix**: source `bin/foam` explicitly: `. /usr/lib/openfoam/openfoam2412/bin/foam`.
2. **`git not found`**: `generateGitInfo.sh` was already written to handle this
   gracefully (returns empty commit/branch when git is absent). No action needed.
3. **`make: No rule to make target 'compile'`**: The tarball was almost certainly
   extracted into a **fresh/empty directory** instead of an existing `aerosolved/`
    checkout. Since the tarball is a *delta* (only 37 files), it does **NOT**
    include the top-level `Makefile`, which defines the `compile` target.
    **Fix**: `cd` *into the existing `aerosolved/`* directory, *then* extract
    the tarball there. Verify `Makefile` is present: `ls Makefile Allwmake`.
4. **aarch64**: No issue. The code is architecture-independent.
   Building natively on ARM is fine — no QEMU emulation needed.

### Steps for the next run
1. Extract the **new** `aerosolved.changes.tar.gz` (contains v2 harness):
      ```bash
     cd <aerosolved-existing-checkout>/
    tar -xzf aerosolved.changes.tar.gz      # extract here, into the repo
    ```
2. Verify the top-level build files are present:
      ```bash
    ls Makefile Allwmake scripts/AeroSolvedRunFunctions
      ```
3. Source OpenFOAM:
      ```bash
     . /usr/lib/openfoam/openfoam2412/bin/foam
      ```
4. Run the v2 harness:
      ```bash
     bash build_and_smoketest.sh
      ```
5. If it fails, paste `/tmp/aerosolved_build_*.log` tail (30 lines) back.

---

## 6. RESOLVED root cause of the `make compile` failure (from `aerosolved_build_20260901_171001.log`)

**Symptom:**
```
aerosolModel/aerosolModelGitInfo.H:3:12: error: stray '\' in program
     3 |     Info<< \"Git state at compilation:\" << nl << nl
...
/usr/bin/aarch64-linux-gnu-ld.bfd: cannot find -laerosolModels
```

**Root cause:** the `generateGitInfo.sh` that was producing
`aerosolModelGitInfo.H` emitted **over-escaped C++ string literals**
(`Info<< \\"...\\\\"` — a stray backslash before every quote). g++ reads `\"` as a stray
backslash + an unterminated string, so every header line fails to compile. Because
`aerosolModel.C` does `#include "aerosolModelGitInfo.H"` at line 21, the *first* object
of `libaerosolModels` fails, the shared library never builds, and **every downstream
solver/utility dies at link time** with `cannot find -laerosolModels`
(aerosolEulerFoam, aerosolBuoyantEulerFoam, setLogNormal, setSaturatedMixture, …).

**Why it happened:** the script wrote the header with `echo "…\\\"…"` — whether that
expands to `\"` (bad) or `"` (good) depends on the shell's `xpg_echo`/`POSIXLY_CORRECT`
mode, which differs between distros. On the aarch64 box it produced `\"`.

**The fix (shipped in this tarball):** `libraries/aerosolModels/generateGitInfo.sh` now
emits the header with `printf '%s\n' '…'`, which does **no backslash interpretation** —
the output is byte-for-byte valid C++ on every platform, independent of `xpg_echo`.
Verified locally: the generated header compiles under `clang++ -std=c++17
-Wall -Wextra` with no diagnostics (both the git-present and git-absent cases). The script
is also now **idempotent** (commit-sidecar), so it no longer thrashes on every `.C` build.

**The "PASSED" that looked real was a false positive.** With `SKIP_BUILD=1` the *second*
run found `log.blockMesh`/`log.aerosolEulerFoam` left over from the *first* (failed) run;
`Allrun`'s "already run … remove log file to re-run" guard short-circuited `exit 0`, so
the old harness printed "OK" without the solver ever running. The **v3 harness** (now
shipped) (a) deletes the case's `log.*` before each run, (b) treats an "already run … to
re-run" short-circuit as **FAILED**, not PASSED, and (c) warns when a solver is absent
from `PATH` under `SKIP_BUILD=1`. A genuine "PASSED" now requires the solvers to be
present *and* the run to actually execute.

---

## 7. Known environment caveat for the local agent session

The `python3` on `PATH` (`/usr/local/bin/python3`, Python 3.14) inherits
`PYTHONPATH=<hermes-venv>` — a different Python. Always run with the
explicit venv path, or activate it first:
```bash
/Users/aspen_mb/Documents/VS_CODE_PROJECTS/Aerosolved/.venv-aerosolved/bin/python3 <script>
```

---

## 8. Summary

| Status | Count |
|--------|-------|
| ✅ Verified locally | 8 Tier-1 items + v3 harness |
| ✅ Verified locally | All 38 files compile / byte-compile / bash -n clean |
| ✅ RESOLVED | `generateGitInfo.sh` over-escape — now `printf`-based, byte-for-byte valid C++ |
| ✅ **VERIFIED ON aarch64 OF 2412** | **Full `Allwmake` build (libaerosolModels + solvers) links; `rain` smoke case runs end-to-end** |
| ➕ Deferred (needs build-time opt-in) | `*New.C` deletion, coalescence model, awk Allrun — all documented |

**Verified on `aspen_vm` (Linux 7.0.0-30-generic aarch64, OpenFOAM 2412 distro, `/usr/lib/openfoam/openfoam2412`): 2026-09-01 18:13 PDT — full `Allwmake` build succeeded; `cases/rain` smoke case PASSED (blockMesh → aerosolEulerFoam → postProcess); all 5 compute scripts + 2 plot scripts ran on the project venv; Allrun VARS templates passed. Log: `/tmp/aerosolved_build_20260901_181329.log`. Tarball sha256 `de797aac…` (38 files, 31 modified + 7 new) contains the fix, the self-diagnosing harness, `bootstrap_openfoam.sh`, `requirements.txt`, `IMPROVEMENTS.txt`, `doc/Chap0_Installation.md`, `.github/workflows/ci.yml`, and the Tier-1 source fixes.**

Re-extracting the tarball into an **existing** `aerosolved/` checkout picks up all changes; the harness leaves no build cruft in the case directories.

