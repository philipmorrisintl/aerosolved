#!/bin/bash
# bootstrap_openfoam.sh — land the changeset and run the smoke test on any
# OpenFOAM box WITHOUT hard-coded paths. Run this from anywhere on the box.
#
# It:
#   1. finds your aerosolved source root (by marker files), or a path you pass,
#   2. overlays aerosolved.changes.tar.gz into it,
#   3. activates OpenFOAM by trying every common activation method,
#   4. runs build_and_smoketest.sh.
#
# Usage:
#   bash bootstrap_openfoam.sh                         # auto-detect everything
#   bash bootstrap_openfoam.sh /path/to/aerosolved     # if auto-detect misses
#   bash bootstrap_openfoam.sh /path/to/aerosolved TARBALL.tgz
#
# If it cannot find the OpenFOAM activation, it prints the discovery commands
# you can run to find the exact one for THIS box (paste the output back).

set -uo pipefail

# ---- 1. Locate the aerosolved root ----
# An aerosolved root has: Makefile + Allwmake + cases/rain/Allrun together.
find_root() {
    for base in "${1:-$PWD}" "$PWD" "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; do
        find "$base" -maxdepth 4 -type f -name Allwmake 2>/dev/null | while read -r aw; do
            d="$(dirname "$aw")"
            if [ -f "$d/Makefile" ] && [ -f "$d/cases/rain/Allrun" ]; then
                echo "$d"; return 0
            fi
        done
    done
}

ROOT_HINT="${1:-}"
AER_ROOT="$(find_root "$ROOT_HINT")"
if [ -z "${AER_ROOT:-}" ]; then
    echo "Could not auto-locate an aerosolved checkout."
    echo "Pass it explicitly:  bash bootstrap_openfoam.sh /path/to/aerosolved"
    echo "Find your clone with:  find ~ -maxdepth 3 -type f -name Allwmake 2>/dev/null"
    exit 1
fi
echo "[1/4] aerosolved root: $AER_ROOT"
cd "$AER_ROOT" || { echo "cannot cd to $AER_ROOT"; exit 1; }

# ---- 2. Overlay the changeset ----
TARBALL="${2:-aerosolved.changes.tar.gz}"
case "$TARBALL" in
    /*) : ;;                       # absolute
    *)  for cand in "$TARBALL" "$AER_ROOT/$TARBALL" "../$TARBALL"; do \
           [ -f "$cand" ] && TARBALL="$cand" && break; \
        done ;;
esac
if [ -f "$TARBALL" ]; then
    echo "[2/4] overlaying $TARBALL"
    tar -xzf "$TARBALL" -C "$AER_ROOT"
else
    echo "[2/4] WARN: $TARBALL not found — assuming changes already applied"
fi

# ---- 3. Activate OpenFOAM (try every common method) ----
activate_openfoam() {
    if command -v wmake >/dev/null 2>&1; then echo "OpenFOAM already active"; return 0; fi
     # Debian/Ubuntu: package may provide its own env script
    for cand in \
        /usr/share/openfoam/openfoam2412/etc/bashrc \
        /usr/lib/openfoam/openfoam2412/etc/bashrc \
        /opt/openfoam/openfoam-2412/OpenFOAM-2412/etc/bashrc \
        "$WM_PROJECT_DIR/etc/bashrc"; do
        if [ -f "$cand" ]; then
            echo "sourcing env script: $cand"
            . "$cand" 2>/dev/null && command -v wmake >/dev/null 2>&1 && return 0
        fi
    done
     # distro `foam` wrapper on PATH
    for cand in /usr/bin/openfoam /usr/bin/foam /usr/local/bin/foam; do
        [ -x "$cand" ] && { "$cand" >/dev/null 2>&1 && command -v wmake >/dev/null 2>&1 && return 0; }
    done
    return 1
}

if ! activate_openfoam; then
    echo
    echo "[3/4] OpenFOAM NOT activated automatically. Discover its activation:"
    echo "   A) which wmake RunFunctions getApplication 2>/dev/null"
    echo "   B) dpkg -L openfoam2412 2>/dev/null | grep -E 'bashrc|/foam$|etc/' | head"
    echo "   C) ls -d /usr/lib/openfoam/openfoam2412*/bin /opt/openfoam/*/OpenFOAM-*/etc 2>/dev/null"
    echo "   D) echo \$WM_PROJECT_DIR  \$FOAM_INST_DIR  \$WM_BINARY_DIR"
    echo "   E) find / -name bashrc -path '*openfoam*' 2>/dev/null | head"
    echo "Paste that output back so I can pin the exact activation line for this box."
    exit 1
fi
echo "[3/4] OpenFOAM active: wmake=$(command -v wmake), FOAM_RELEASE=${FOAM_RELEASE:-<unset>}"

# ---- 4. Run the smoke test ----
echo "[4/4] running build_and_smoketest.sh"
bash "$AER_ROOT/build_and_smoketest.sh" "$@"
