"""Backend layer: pure Python, no Qt.

Submodules:
   case_config  — case discovery + Allrun/VARS static analysis
   m4_expand    — `m4 -DVAR=...` macro expansion runner
   runner       — Allrun orchestration (mesh → solve → postProcess) + events
   results      — per-shape post-processing parsers + dashboard build
   dashboard    — a data-only "dashboard spec" the frontend renders

All of this is importable and testable without a display.
"""

__all__ = ["case_config", "m4_expand", "runner", "results", "dashboard"]

__version__ = "0.1.0"
