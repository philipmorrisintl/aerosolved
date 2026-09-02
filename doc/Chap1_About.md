# About

_Navigation_

1. [Installation](Chap0_Installation.md)
2. [About](Chap1_About.md)
3. [Model](Chap2_Model.md)
4. [Tutorial](Chap3_Tutorial.md)
5. [Cases](Chap4_Cases.md)
6. [Nomenclature](Chap5_Nomenclature.md)
7. [Classes](Chap6_Classes.md)
8. [References](Chap7_References.md)

## Welcome to AeroSolved

AeroSolved is an OpenFOAM library offering tools to model aerosols in an Eulerian way. The library relies fully on OpenFOAM but implements most key aspects of the models and methods by itself.

Originally, AeroSolved was developed jointly by Philip Morris International R&D (PMI R&D) and the Department of Applied Mathematics, University of Twente (UT), The Netherlands.

The main goal of AeroSolved is to offer a platform for the simulation of aerosol dynamics, including:

* **Aerosol formation** through supersaturation and subsequent nucleation
* **Aerosol evolution** through condensation/evaporation and coalescence
* **Aerosol deposition** due to dispersed diffusion or dispersed inertial drift

A key part of AeroSolved is that it models the particle size distribution, such that detailed information on particle size can be predicted locally. This is done using two independent methods:

* A **moment model**, based on two moments of the log-normal size distribution
* A **fixed sectional model**, fully discretizing the size space

The most important parts of AeroSolved are:

* The **aerosolEulerFoam solver**: a solver based on reactingFoam and incorporating Eulerian aerosol models
* The **aerosolModels library**: contains the implementation of different aerosol models such as the fixedSectional and twoMomentLogNormal models. The main purpose of aerosolModels library is the modeling of the particle size distribution. It relies on various submodels such as nucleation, condensation, and coalescence
* The **aerosolThermo library**: a thermo package that is based on psiThermo, and contains two separate thermo libraries (which are each based on rhoThermo) for the continuous and dispersed phases. The purpose of the aerosolThermo library is to combine the continuous and dispersed thermo libraries in order to create a mixture thermo library. This startegy is following the twoPhaseMixtureThermo library of OpenFOAM's standard compressibleInterFoam solver.

## Installation

AeroSolved is an OpenFOAM library package. Full installation, build and
post-processing instructions are in the
[Installation](Chap0_Installation.md) chapter. In short: install a supported
OpenFOAM (see the [Dependencies](#dependencies) section), then

```bash
git clone https://github.com/philipmorrisintl/aerosolved.git
cd aerosolved
make            # compile the libraries, solvers and utilities
```

Each case under `cases/` is self-contained and is run with `./Allrun` (for
example `./Allrun fullStokes sectional` in `cases/bentPipe`).

## Dependencies

AeroSolved has no special dependencies other than OpenFOAM. It is developed and tested against **OpenFOAM-v2406** and **OpenFOAM-v2412** (see the [README](../README.md) and [Installation](#installation)). Some cases require python3 and numpy for generating a post-processing plot (see the [attribution note](../AttributionNote)). The availability of python3 is tested in the Allrun scripts.

## Documentation

From the header of this page you can navigate the model, cases, nomenclature, and classes documentation. As with standard OpenFOAM applications and libraries, further documentation is available in the headers of important `.H` files of the source code. This source code documentation is also parsed by Doxygen. We recommend generating the Doxygen documentation as explained in the Installation section.
