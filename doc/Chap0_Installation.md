# Installation

_Navigation_

1. [Installation](Chap0_Installation.md)
2. [About](Chap1_About.md)
3. [Model](Chap2_Model.md)
4. [Tutorial](Chap3_Tutorial.md)
5. [Cases](Chap4_Cases.md)
6. [Nomenclature](Chap5_Nomenclature.md)
7. [Classes](Chap6_Classes.md)
8. [References](Chap7_References.md)

## 1. Overview

AeroSolved is an OpenFOAM library package. It builds and runs **on top of** an
installed OpenFOAM environment; it does not bundle OpenFOAM. Installation has
four stages:

1. install **OpenFOAM** (a supported version, see §2)
2. obtain the **AeroSolved** sources
3. **build** the libraries and solvers with `make`
4. (optional) install the **Python** dependencies for post-processing

## 2. Prerequisites

### 2.1 OpenFOAM version

AeroSolved is developed and tested against **OpenFOAM-v2406** and
**OpenFOAM-v2412** (OpenFOAM.com / ETSI). Other versions may work but are
unsupported; the `Make/options` files reference library and header paths that
only exist in recent OpenFOAM. Obtain OpenFOAM from
<https://www.openfoam.com/> or, for a reproducible Linux build, the official
OpenCFD container image `opencfd/openfoam-docker`
(see <https://www.opencfd.co.uk/>).

> **Note on the About page:** the older "Dependencies" section of
> `Chap1_About.md` historically referenced `OpenFOAM-v1812`. The current,
> supported versions are **v2406 / v2412**.

### 2.2 Command-line tools

| Tool       | Needed for                                   |
|------------|----------------------------------------------|
| `bash`     | all `Allrun` scripts                         |
| `m4`       | expanding the `.m4` case-template files      |
| `git`      | optional; embeds commit info into the solver |
| `python3`  | post-processing / plot scripts               |

On Debian/Ubuntu: `sudo apt-get install build-essential gfortran m4 git`.

### 2.3 Python packages (post-processing only)

The `Allrun` scripts and `AeroSolvedRunFunctions` detect `python3` and warn
(and will hard-fail the post-processing step) if `numpy` and `matplotlib` are
missing. `scripts/functions.py` additionally requires `scipy`.

Recommended — install into a **dedicated virtual environment** so AeroSolved
does not clash with the system Python (and vice versa):

```bash
python3 -m venv .venv-aerosolved
source .venv-aerosolved/bin/activate        # Linux / macOS
pip install -r requirements.txt              # numpy, scipy, matplotlib
```

On Debian/Ubuntu without a venv: `sudo apt-get install python3-numpy
python3-matplotlib python3-scipy`.

## 3. Obtain the sources

```bash
git clone https://github.com/philipmorrisintl/aerosolved.git
cd aerosolved
```

The source tree layout is:

```
aerosolved/
├── applications/   solvers (aerosolEulerFoam, aerosolBuoyantEulerFoam) + utilities
├── libraries/      aerosolThermo, aerosolModels, customFunctions,
│                   customTurbulenceModels
├── cases/          example / tutorial cases
├── scripts/        setup & post-processing Python helpers
├── doc/            documentation
├── Allwmake        build everything
├── Allwclean       clean build artefacts
├── Makefile        make / make clean / make doc
└── requirements.txt
```

## 4. Build

1. **Activate OpenFOAM** for your shell:

   ```bash
   # OpenFOAM.com install (example):
   source /opt/openfoam2412/etc/bashrc
   # ETSI install (example):
   # source ~/openfoam/OpenFOAM-v2412/etc/bashrc
   ```

   This defines `WM_PROJECT_DIR`, `FOAM_APP_BIN`, `FOAM_APP_LIB`,
   `FOAM_APP_CASES`, and sets up the compiler wrapper `wmake`.

2. **Compile everything:**

   ```bash
   make            # => compile + doc
   ```

   or just the code, skipping Doxygen:

   ```bash
   make compile    # => ./Allwmake
   make clean      # => ./Allwclean
   ```

   `make compile` builds, in order: `customFunctions`,
    `customTurbulenceModels`, `aerosolThermo`,
    `aerosolModels`, the two solvers, the two utilities
    (`setSaturatedMixture`, `setLogNormal`), and the test app. Each library and
   solver lands in your OpenFOAM user-app area (`$FOAM_APP_LIB`,
    `$FOAM_APP_BIN`).

3. **Verify the build:**

   ```bash
   which aerosolEulerFoam        # should print the built solver
   aerosolEulerFoam -help 2>&1 | head  # (prints its header / git info)
   ```

> **On macOS / BSD:** the build depends on GNU coreutils for parts of the
> OpenFOAM toolchain itself. AeroSolved's own build-time helper
> `libraries/aerosolModels/generateGitInfo.sh` was made portable (no GNU-only
> `find -mmin` / `touch -d`), but the wider OpenFOAM toolchain is Linux-first.
> Use the OpenCFD container for reproducible CI or macOS builds.

## 5. Run a case

Each folder under `cases/` is self-contained and ships an `Allrun` script that
prepares the case, generates the mesh, runs the solver, and (for some cases)
posts the result. Example — the bent-pipe deposition case, using the full Stokes
drift model and the fixed-sectional aerosol model:

```bash
cd cases/bentPipe
./Allrun fullStokes sectional
```

| Argument 1 (inertial drift) | Argument 2 (aerosol model) |
|-----------------------------|----------------------------|
| `fullStokes`                | `sectional`                |
| `Manninen`                  | `moment`                   |

Run `./Allrun` without arguments to see the usage, and `./Allclean` to remove
generated artefacts (mesh, `0`, `postProcessing`).

## 6. Common problems

### 6.1 `Could not find python3 binary`

`checkPython3` in `scripts/AeroSolvedRunFunctions` hard-exits when `python3` is
not on `PATH`. Ensure a `python3` (with `numpy`) is available, as in §2.3.

### 6.2 `ModuleNotFoundError: No module named 'numpy._core._multiarray_umath'`

This does **not** mean numpy is broken. It almost always means a numpy built for
a *different Python version or ABI* is being imported by the wrong interpreter.
Typical causes:

- **Mixing Python versions.** E.g. a `PYTHONPATH` (or `VIRTUAL_ENV`) points at a
  virtualenv built for **Python 3.11**, while the `python3` on your `PATH` is
  **Python 3.14** (or 3.12/3.13). The 3.11 C-extension cannot load into the
  3.14 interpreter.
- **A stale venv** whose numpy was installed against a Python that no longer
  exists on the box.

Diagnosis:

```bash
python3 -c "import sys; print(sys.version)"          # which interpreter
python3 -c "import numpy; print(numpy.__version__)"  # does it import at all?
echo "$PYTHONPATH" "$VIRTUAL_ENV"                    # are they pointing elsewhere?
```

Fix: use one coherent interpreter. Either

```bash
# unset the foreign environment and use the system python
unset PYTHONPATH VIRTUAL_ENV
python3 -m venv .venv-aerosolved && source .venv-aerosolved/bin/activate
pip install -r requirements.txt
```

or run the scripts through the venv python explicitly:

```bash
.venv-aerosolved/bin/python3 scripts/CAG.py
```

### 6.3 `m4: command not found`

Install `m4` (§2.2). The `.m4` template files in `cases/` and
`constant/`/`system/` directories are expanded by `setMacros` into the
runtime configuration files.

### 6.4 Build fails with `wmake: command not found`

OpenFOAM was not activated (step 1 of §4). `wmake` / `wmakeLnInclude` /
`wclean` are only available after sourcing the OpenFOAM `bashrc`.

## 7. Documentation build (optional)

If [Doxygen](https://www.stack.nl/~dimitri/doxygen) and
[Graphviz](http://graphviz.org/) are installed:

```bash
make doc
# HTML API reference under doc/output/html/
```
