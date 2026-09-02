#!/bin/bash
# build_and_smoketest.sh
# One-shot build + smoke-test harness for AeroSolved on a machine WITH OpenFOAM.
#
# Usage:
#     cd <aerosolved checkout with our changes>
#     bash build_and_smoketest.sh             # build everything + run smoke case
#     CASES=rain,cavity bash build_and_smoketest.sh
#     SKIP_BUILD=1 bash build_and_smoketest.sh    # skip rebuild, just run case(s)
#
# Requirements on the target machine:
#    * OpenFOAM v2406 or v2412 activated
#    * build-essential, m4
#    * python3 + numpy + scipy + matplotlib (see requirements.txt)
#
# Exit status: 0 = all passed, non-zero on first failure.
# This script is self-contained and leaves no build cruft in the case dirs;
# run `./Allclean` in each case afterwards if desired.

set -uo pipefail

# =========================================================================
# 0. Resolve project root and print diagnostics up front
# =========================================================================
ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"

echo "=================================================================="
echo " AeroSolved build + smoke-test harness (v2 - self-diagnosing)"
echo " time : $(date)"
echo " host : $(uname -srm)"
echo "====================================================="
echo
echo "[0/6] environment diagnostics"
echo "  CWD   : $ROOT"
echo "  ROOT  : $ROOT"
ls -la "$ROOT"/Makefile "$ROOT"/Allwmake 2>/dev/null || {
    echo "  ERROR: Makefile and/or Allwmake not found! This script must be"
    echo "         inside an aerosolved checkout, OR the tarball must be"
    echo "         extracted on top of an existing aerosolved/ directory."
    echo "         (The 37-file aerosolved.changes.tar.gz is a DELTA — it"
    echo "         does NOT include Makefile, Allwmake, etc.)"
    exit 1
}

# OpenFOAM detection — try multiple layouts
OF_ACTIVATED=""
if [ -n "${WM_PROJECT_DIR:-}" ] && [ -f "${WM_PROJECT_DIR}/bin/foam" ]; then
    OF_ACTIVATED="$(cd "$(dirname "${WM_PROJECT_DIR}/bin/foam")/.." && pwd)"
elif [ -n "${FOAM_INST_DIR:-}" ] && [ -f "$FOAM_INST_DIR/bin/foam" ]; then
    OF_ACTIVATED="$(cd "$(dirname "$FOAM_INST_DIR/bin/foam")/.." && pwd)"
elif [ -f "/usr/lib/openfoam/openfoam2412/bin/foam" ]; then
    OF_ACTIVATED="/usr/lib/openfoam/openfoam2412"
elif [ -f "/opt/openfoam/openfoam2412/bin/foam" ]; then
    OF_ACTIVATED="/opt/openfoam/openfoam2412"
fi

echo "  WM_PROJECT_DIR : ${WM_PROJECT_DIR:-<unset>}"
echo "  FOAM_RELEASE   : ${FOAM_RELEASE:-<unset>}"
echo "  FOAM_INST_DIR  : ${FOAM_INST_DIR:-<unset>}"
echo "  wmake          : $(command -v wmake 2>/dev/null || echo 'NOT FOUND')"
echo "  m4             : $(command -v m4 2>/dev/null || echo 'NOT FOUND')"
echo "  git            : $(command -v git 2>/dev/null || echo 'NOT FOUND — OK, generateGitInfo.sh handles it')"
echo "  make           : $(command -v make 2>/dev/null || echo 'NOT FOUND')"
echo "  python3        : $(command -v python3 2>/dev/null || echo 'NOT FOUND')"

# Try to activate OpenFOAM if it's not active yet
if [ -z "$OF_ACTIVATED" ] && [ "${SKIP_OFS:-0}" != 1 ]; then
    echo
    echo "  WARNING: NO OpenFOAM installation found!"
    echo "  Searched:"
    echo "     \$WM_PROJECT_DIR/bin/foam"
    echo "     \$FOAM_INST_DIR/bin/foam"
    echo "     /usr/lib/openfoam/openfoam2412/bin/foam"
    echo "     /opt/openfoam/openfoam2412/bin/foam"
    echo "  Activate OpenFOAM, e.g.:"
    echo "     . /usr/lib/openfoam/openfoam2412/bin/foam"
    echo "  Or set WM_PROJECT_DIR correctly before running this script."
    echo "  (SKIP_OFS=1 skips the OpenFOAM check — run Python + VARS steps only)"
    echo "  Exiting... (set SKIP_OFS=1 to skip and run Python/VARS steps anyway)"
    exit 1
elif [ -n "$OF_ACTIVATED" ]; then
    echo "  OpenFOAM found at: $OF_ACTIVATED"
    if ! { wmake --version >/dev/null 2>&1 || echo "$WM_PROJECT_DIR" | grep -q "$OF_ACTIVATED"; }; then
        echo "  Activating OpenFOAM: sourcing $OF_ACTIVATED/bin/foam"
         . "$OF_ACTIVATED/bin/foam" 2>/dev/null || {
            echo "  WARNING: could not source $OF_ACTIVATED/bin/foam"
            echo "           wmake still available: $(command -v wmake 2>/dev/null || echo 'NO')"
        }
    fi
    echo "  After activation:"
    echo "    FOAM_RELEASE : ${FOAM_RELEASE:-<still unset>}"
    echo "    FM_PROJECT_DIR: ${WM_PROJECT_DIR:-<still unset>}"
    echo
fi

# =========================================================================
# 1. Environment sanity
# =========================================================================
echo "[1/6] checking prerequisites"
if [ "${SKIP_OFS:-0}" != 1 ]; then
    command -v wmake >/dev/null 2>&1 || { echo "ERROR: wmake not found — OpenFOAM not properly activated"; exit 1; }
    command -v m4     >/dev/null 2>&1 || { echo "ERROR: m4 not found — apt-get install m4"; exit 1; }
    command -v make   >/dev/null 2>&1 || { echo "ERROR: make not found — apt-get install build-essential"; exit 1; }
else
    echo "  SKIP_OFS=1: skipping wmake/make prerequisite gate"
fi

# =========================================================================
# 2. Build
# =========================================================================
BUILD_LOG="/tmp/aerosolved_build_$(date +%Y%m%d_%H%M%S).log"
echo "[2/6] building"
echo "  log: $BUILD_LOG"
echo "  running from: $ROOT"
echo "  Makefile exists: $(test -f "$ROOT/Makefile" && echo YES || echo NO)"

if [ "${SKIP_BUILD:-0}" = 1 ]; then
    echo "  SKIP_BUILD=1 — skipping build, using existing executables"
    if ! command -v aerosolEulerFoam >/dev/null 2>&1; then
        echo "  WARNING: aerosolEulerFoam not in PATH — case run will fail"
    fi
else
    # Try `make compile` first, fall back to `./Allwmake`
    build_ok=0
    if [ -f "$ROOT/Makefile" ] && grep -q '^compile:' "$ROOT/Makefile" 2>/dev/null; then
        echo "  strategy: make compile"
        time make compile 2>&1 | tee "$BUILD_LOG"
        if [ $? -eq 0 ]; then build_ok=1; fi
    fi
    # Always fall back to Allwmake if build was not ok
    if [ "$build_ok" -eq 0 ] && [ -f "$ROOT/Allwmake" ]; then
        if [ "$build_ok" -eq 0 ] && [ -f "$ROOT/Makefile" ]; then
            echo "  (make compile failed or had no 'compile' target; falling back to ./Allwmake)"
        fi
        echo "  strategy: bash Allwmake"
        time bash ./Allwmake 2>&1 | tee "${BUILD_LOG}.allwmake"
        if [ $? -eq 0 ]; then build_ok=1; fi
    fi
    if [ "$build_ok" -eq 0 ]; then
        echo
        echo "  ERROR: build failed. See log files:"
        echo "    $BUILD_LOG"
        echo "    ${BUILD_LOG}.allwmake (if exists)"
        echo "  Last 30 lines of the log:"
        tail -n 30 "$BUILD_LOG" 2>/dev/null || tail -n 30 "${BUILD_LOG}.allwmake" 2>/dev/null
        exit 1
    fi
    echo "  build OK — checking for built executables:"
    ls -la $(find . -name "aerosolEulerFoam" -path "*/platforms/*" 2>/dev/null | head -3) 2>/dev/null
fi

# =========================================================================
# 3. Python smoke run (independent of build)
# =========================================================================
echo
echo "[3/6] running self-contained Python scripts"
# Prefer the project venv (has numpy/scipy/matplotlib); it may live in the
# repo dir or the parent dir (parent level, i.e. outside aerosolved/). Fall
# back to whatever `python3` is on PATH. Plotting scripts NEED matplotlib;
# when it's absent the plot step is reported as SKIPPED, not FAILED.
PY="$(command -v python3 2>/dev/null || echo python3)"
for cand in ../.venv-aerosolved/bin/python3 .venv-aerosolved/bin/python3; do
    if [ -x "$cand" ]; then PY="$cand"; break; fi
done

# Use GNU `timeout` for the 60s/300s caps where available (Linux coreutils);
# degrade gracefully on systems that lack it (e.g. bare macOS -> no cap).
_has_timeout=0
if command -v timeout >/dev/null 2>&1; then
    _has_timeout=1
else
    echo "   (note: 'timeout' not found — running steps without a hard time cap)"
fi
run_py() {
    # run_py SEC SCRIPT args...  (time-capped run for the Python step)
    if [ "$_has_timeout" = 1 ]; then
        timeout "${1}" "${@:2}"
    else
        "${@:2}"
    fi
}
run_case() {
    # run_case SEC CMD...  (time-capped run for the OpenFOAM case step)
    if [ "$_has_timeout" = 1 ]; then
        timeout "${1}" "${@:2}"
    else
        "${@:2}"
    fi
}

_mpl_ok=1; "$PY" -c "import matplotlib" 2>/dev/null || _mpl_ok=0
echo "  interpreter : $PY"
[ "$_mpl_ok" = 1 ] || echo "  note: matplotlib missing in this interpreter — plotting scripts will be SKIPPED"
echo "            (run the project venv: ../.venv-aerosolved/bin/python3)"

python_rc=0
# Pure-compute scripts (always run; report failures as errors).
for s in CAG.py Cloud.py saturatedMixtureFlow.py uniformNucleationMixture.py; do
    echo "      --- $s ---"
    run_py 60 "$PY" scripts/"$s" 2>&1 | tail -3
    rc=${PIPESTATUS[0]}
    [ $rc -ne 0 ] && { echo "  ERROR: $s failed (rc=$rc)"; python_rc=1; }
done
# Plotting scripts (need matplotlib; SKIPPED when it's absent, not FAILED).
if [ "$_mpl_ok" = 1 ]; then
    for s in logNormalSectionalGrid.py blendedCoalescence.py; do
        echo "      --- rendering $s ---"
        MPLBACKEND=Agg run_py 60 "$PY" scripts/"$s" 2>&1 | tail -3
        rc=${PIPESTATUS[0]}
          [ $rc -ne 0 ] && { echo "  ERROR: $s failed (rc=$rc)"; python_rc=1; }
    done
else
    echo "      --- skipping plot scripts (matplotlib not in $PY) ---"
fi
if [ $python_rc -eq 0 ]; then
    echo "  all runnable python scripts OK"
else
    echo "  SOME python scripts failed (see above)"
fi

# =========================================================================
# 4. VARS block m4 dry-run
# =========================================================================
echo
echo "[4/6] verifying Allrun VARS blocks"
if command -v m4 >/dev/null 2>&1; then
    var_ok=0
    for f in $(find cases -name Allrun | head -15); do
        dir="$(dirname "$f")"
        # Check each .m4 that the Allrun might reference — best effort
        for m4f in "$dir"/constant/*.m4 "$dir"/0.org/*.m4 "$dir"/0/*.m4; do
            [ -f "$m4f" ] || continue
             # Dummy substitution: pass all likely vars
            dummy_m4_out=$(m4 -P -D"MV=sectional" -D"MINLETBC=sectionalLogNormal" \
                -D"MWALLBC=sectionalMixedAbsorbing" -D"ZWALLBC=massFracFromSectional" \
                -D"MODEL=fixedSectional" -D"MBC=sectional" "$m4f" 2>&1)
            if echo "$dummy_m4_out" | grep -qi "error\|fatal\|unbound\|undefined"; then
                # Check if it's a m4 error (not just a macro-not-found in a case that
                # doesn't use that macro — skip that case)
                m4_error=$(echo "$dummy_m4_out" | grep -i "error\|fatal" | head -1)
                if echo "$m4_error" | grep -q "syntax\|parse\|unknown function\|bad token"; then
                    echo "  M4 ERROR in $(basename "$f") / $(basename "$m4f"): $m4_error"
                    var_ok=1
                fi
            fi
        done
    done
    [ $var_ok -eq 0 ] && echo "  Allrun VARS templates OK" || echo "  Some VARS templates had errors (see above)"
else
    echo "  m4 not found — skipping VARS check"
fi

# =========================================================================
# 5. Smoke case(s)
# =========================================================================
echo
echo "[5/6] smoke case(s)"
SMOKE="${CASES:-rain}"
IFS=',' read -ra CASE_LIST <<< "$SMOKE"

# Per-case args via a plain case (portable to bash 3.2; no associative arrays).
# Override any case with an env var, e.g. CASE_ARGS_rain="fullStokes sectional".
case_args_for() {
    case "$1" in
        rain)         printf 'fullStokes sectional' ;;
        bentPipe)     printf 'fullStokes sectional' ;;
        uniformCoalescence) printf 'sectional' ;;
        Brownian)     printf 'sectional' ;;
        saturatedBox) printf '' ;;
        *)            printf '%s' "${CASE_ARGS:-}" ;;   # generic fallback
    esac
}

smoke_rc=0
first=1
for raw_case in "${CASE_LIST[@]}"; do
    case_raw="$(echo "$raw_case" | xargs -I{} echo {})"
     # Per-case argument, with per-case override supported.
    _arg_env="CASE_ARGS_${case_raw}"
    if [ -n "${!_arg_env:-}" ]; then
        case_arg="${!_arg_env}"
    else
        case_arg="$(case_args_for "$case_raw")"
    fi

    full="cases/$case_raw"
    if [ "$first" = 1 ]; then echo "   --- $full/Allrun (args: $case_arg) ---"; first=0; fi
      [ -d "$full" ] || { echo "ERROR: $full not found, skipping"; smoke_rc=1; continue; }
      [ -x "$full/Allrun" ] || { echo "ERROR: no executable Allrun in $full, skipping"; smoke_rc=1; continue; }
    out_dir="/tmp/aerosolved_smoke_$(date +%Y%m%d_%H%M%S)_$(echo "$case_raw" | tr ',' _)"
    LOGF="$out_dir.log"
    mkdir -p "$out_dir" 2>/dev/null || true

        # IMPORTANT: a fresh smoke run must not depend on leftover log.* files
        # (which cause Allrun to short-circuit with exit 0 when blockMesh /
        # aerosolEulerFoam are "already run"). Clean the case's runtime
        # artifacts so Allrun actually re-runs from a known state.
    (cd "$full" && \
      rm -f log.blockMesh log.aerosolEulerFoam \
            log.aerosolBuoyantEulerFoam \
            log.postProcess log.sampleFields \
            log.setMesh log.checkMesh 2>/dev/null) 2>/dev/null || true

        # If SKIP_BUILD=1, warn (don't fail up front) if a solver a smoke run
        # typically needs is missing from PATH — a build was skipped, so the
        # solvers may not exist. This is a warning only; the case run below is
        # the real check.
        if [ "${SKIP_BUILD:-0}" = 1 ]; then
        for solver in blockMesh aerosolEulerFoam; do
           if ! command -v "$solver" >/dev/null 2>&1; then
               echo "  WARNING: $solver not in PATH — a skipped build means this case may fail"
           fi
        done
        fi

        # Run from the case dir; capture full log. The Allrun is self-contained
        # (sources its own setMacros via AeroSolvedRunFunctions).
        out="$(cd "$full" && run_case 300 ./Allrun $case_arg 2>&1 | tee "$LOGF")"
        rc=${PIPESTATUS[0]}

        # Verify the smoke actually did something. A short-circuited run
        # (Allrun's "already run" guard) is not a pass.
        marker="--- $case_raw ---"
        if [ "$rc" -ne 0 ]; then
        echo "$marker FAILED (rc=$rc)"
        tail -n 20 "$LOGF" | sed 's/^/       /'
        smoke_rc=1
        elif echo "$out" | grep -q "already run .* to re-run"; then
        echo "$marker short-circuited (Allrun refused to run — log.* cleanup mismatch)"
        echo "       (harness cleans log.* up front; if this still fires, the case has a"
        echo "        different skip mechanism — inspect manually)"
        smoke_rc=1
        else
        # GENUINE PASS. Be self-evident: point at the full log and print the
        # tail of the postProcess/solver log so the run is auditable without
        # pasting /tmp files.
        SOLV_LOG="$full/log.aerosolEulerFoam"
        [ -f "$SOLV_LOG" ] || SOLV_LOG="$(ls -1 "$full"/log.* 2>/dev/null | grep -vE 'blockMesh|postProcess|checkMesh|sampleFields|setMesh' | head -1)"
        echo "$marker OK  (full log: $LOGF)"
        if [ -n "${SOLV_LOG:-}" ] && [ -f "$SOLV_LOG" ]; then
           echo "       ${SOLV_LOG##*/} — last lines:"
           grep -v '^[[:space:]]*$' "$SOLV_LOG" 2>/dev/null | tail -n 8 | sed 's/^/         /'
           echo "       log dir: $full/log.*  |  Time/ dir: $(ls -dl "$full"/Time 2>/dev/null | awk '{print $NF, $5}' || echo 'n/a')"
        fi
        fi
        done
echo "  Smoke case results: $( [ $smoke_rc -eq 0 ] && echo 'all passed' || echo 'some failed' )"

# =========================================================================
# 6. Final report
# =========================================================================
echo
echo "=================================================================="
echo "  FINAL SUMMARY"
echo "  Python scripts : $( [ $python_rc -eq 0 ] && echo 'PASSED' || echo 'FAILED' )"
echo "  VARS templates : $( [ ${var_ok:-0} -eq 0 ] && echo 'PASSED' || echo 'FAILED' )"
echo "  Smoke cases    : $( [ $smoke_rc -eq 0 ] && echo 'PASSED' || echo 'FAILED' )"
echo "=================================================================="
echo
ALL_OK=0
[ $python_rc -eq 0 ] && [ $smoke_rc -eq 0 ] || ALL_OK=1
[ $ALL_OK -eq 0 ] && echo "ALL SMOKE TESTS PASSED" || echo "SOME SMOKE TESTS FAILED — see details above"
echo "  build log: $BUILD_LOG"
exit $ALL_OK
