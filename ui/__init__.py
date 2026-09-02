"""AeroSolved interactive GUI package (backend + frontend).

The GUI is split into two layers so the orchestration logic can be exercised
without a display:

  ui/backend   — pure Python, no Qt. Case discovery, VARS expansion, run
                 orchestration, result parsing, dashboard extraction.
                 Fully unit-testable headless.
  ui/frontend  — PySide6/Qt wiring on top of the backend. Not importable in a
                 headless environment.
"""

__all__ = ["backend", "frontend"]

__version__ = "0.1.0"
