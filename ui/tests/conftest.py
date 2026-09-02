"""Pytest fixtures for the ui/ package tests."""
import os
import sys
import tempfile
import numpy as np
import pytest

UI_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for p in (UI_ROOT, os.path.join(UI_ROOT, "frontend")):
    sys.path.insert(0, p)
REPO_ROOT = os.path.abspath(os.path.join(UI_ROOT, ".."))


@pytest.fixture
def repo_root():
    return REPO_ROOT


@pytest.fixture
def fake_case():
    """Minimal fake case: Allrun + postProcessing output matching rain."""
    src = tempfile.mkdtemp(prefix="aerosolved_fake_case_")
    os.makedirs(os.path.join(src, "system"), exist_ok=True)
    os.makedirs(os.path.join(src, "constant"), exist_ok=True)
    open(os.path.join(src, "system", "controlDict"), "w").write("application aerosolEulerFoam;\n")
    open(os.path.join(src, "constant", "aerosolProperties"), "w").write("fixedSectionalCoeffs{ distribution{ yMin 1E-24; yMax 1E-10; N 10; } }\n")
    allrun = '#!/bin/bash\nset -e\necho "Running blockMesh"\necho "Running setFields on $PWD"\necho "Courant Number: 0.12 deltaT: 0.5"\npython3 - <<"PY"\nimport numpy as np, os\nfor sub in ("numberFlux","massFlux"):\n    os.makedirs("postProcessing/"+sub+"/0", exist_ok=True)\nt = np.linspace(0,1000,31); N = 10\nnp.savetxt("postProcessing/numberFlux/0/patch.top.dat",    np.column_stack([t, np.random.rand(31,N)]))\nnp.savetxt("postProcessing/numberFlux/0/patch.bottom.dat", np.column_stack([t, np.random.rand(31,N)*0.3]))\nnp.savetxt("postProcessing/massFlux/0/patch.top.dat",    np.column_stack([t, np.random.rand(31,N)*1e-9]))\nnp.savetxt("postProcessing/massFlux/0/patch.bottom.dat", np.column_stack([t, np.random.rand(31,N)*3e-10]))\nPY\necho "Done"\n'
    open(os.path.join(src, "Allrun"), "w").write(allrun)
    os.chmod(os.path.join(src, "Allrun"), 0o755)
    return src


@pytest.fixture
def fake_run_dir():
    d = tempfile.mkdtemp(prefix="aerosolved_run_")
    pp = os.path.join(d, "postProcessing")
    os.makedirs(os.path.join(pp, "numberFlux", "0"), exist_ok=True)
    os.makedirs(os.path.join(pp, "massFlux", "0"),       exist_ok=True)
    t = np.linspace(0, 1000, 31); N = 10
    np.savetxt(os.path.join(pp, "numberFlux", "0", "patch.top.dat"),    np.column_stack([t, np.random.rand(31, N)]))
    np.savetxt(os.path.join(pp, "numberFlux", "0", "patch.bottom.dat"), np.column_stack([t, np.random.rand(31, N) * 0.3]))
    np.savetxt(os.path.join(pp, "massFlux", "0", "patch.top.dat"),    np.column_stack([t, np.random.rand(31, N) * 1e-9]))
    np.savetxt(os.path.join(pp, "massFlux", "0", "patch.bottom.dat"), np.column_stack([t, np.random.rand(31, N) * 3e-10]))
    return d


@pytest.fixture
def mock_qt():
    class QtWidgets:
        @classmethod
        def QComboBox(cls):
            class Combo:
                def __init__(self, *a, **k):
                    self._items = []
                def addItem(self, v):
                    self._items.append(v)
                def count(self):
                    return len(self._items)
            return Combo()
        @classmethod
        def QPushButton(cls, *a, **k):
            class B:
                def click(self, *a, **k): pass
            return B()
        @classmethod
        def QPlainTextEdit(cls, *a, **k):
            class T:
                def appendPlainText(self, v): pass
            return T()
        @classmethod
        def QLabel(cls, *a, **k):
            class L:
                def setText(self, v): pass
            return L()
        @classmethod
        def QApplication(cls, *a, **k):
            class App:
                @staticmethod
                def instance(): return None
                @staticmethod
                def exec(): return 11
            return App()
    return type("qt_mock", (), {
         "QtWidgets": QtWidgets,
         "QtGui": type("QtGui", (), {"QColor": lambda *a, **k: None}),
         "QtCore": type("QtCore", (), {}),
         "binding": "mock",
     })()
