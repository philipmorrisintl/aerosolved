"""ui.frontend -- PySide6 GUI layer (guarded; importable headless).

This layer wraps the pure-Python backend (case discovery, VARS knob editing,
run orchestration, dashboard rendering) in a Qt desktop UI. Qt is imported
lazily and guarded so the package can be imported -- and its non-Qt helpers
tested -- on a headless machine without a display server.
"""
from . import app

__all__ = ["app"]
