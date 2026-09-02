#!/usr/bin/env python3
"""AeroSolved GUI launcher.

    python ui/run_gui.py                    # auto-detect the repo root
    python ui/run_gui.py /path/to/aerosolved   # explicit repo root

Requires: PySide6 (or PySide2, fall back on older systems).
"""
import os, sys, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from frontend.app import main as gui_main

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="AeroSolved GUI")
    parser.add_argument("--repo", default=os.environ.get("AEROSOLVED_ROOT",""),
                        help="aerosolved repo root (auto-detect from CWD if omitted)")
    args = parser.parse_args()
    if args.repo:
        os.environ["AEROSOLVED_ROOT"] = os.path.abspath(args.repo)
    sys.exit(gui_main())
