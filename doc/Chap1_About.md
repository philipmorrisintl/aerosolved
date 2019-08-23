# About

_Navigation_

1. [About](Chap1_About.md)
2. [Model](Chap2_Model.md)
3. [Tutorial](Chap3_Tutorial.md)
4. [Cases](Chap4_Cases.md)
5. [Nomenclature](Chap5_Nomenclature.md)
6. [Classes](Chap6_Classes.md)
7. [References](Chap7_References.md)

## Welcome to AeroSolved

AeroSolved is an OpenFOAM library offering tools to model aerosols in an Eulerian way. The library relies fully on OpenFOAM but implements most key aspects of the models and methods by itself.

Originally, AeroSolved was developed jointly by Philip Morris International R&D (PMI R&D) and the Department of Applied Mathematics, University of Twente (UT), The Netherlands.

The main goal of AeroSolved is to offer a platform for the simulation of aerosol dynamics, including:

* **Aerosol formation** through supersaturation and subsequent nucleation
* **Aerosol evolution** through condensation/evaporation and coalescence
* **Aerosol deposition** due to Brownian motion or inertial drift

A key part of AeroSolved is that it models the particle size distribution, such that detailed information on particle size can be predicted locally. This is done using two independent methods:

* A **moment model**, based on two moments of the log-normal size distribution
* A **fixed sectional model**, fully discretizing the size space

The most important parts of AeroSolved are:

* The **aerosolEulerFoam solver**: a solver based on reactingFoam and incorporating Eulerian aerosol models
* The **aerosolModels library**: contains the implementation of different aerosol models such as the fixedSectional and twoMomentLogNormal models. The main purpose of aerosolModels library is the modeling of the particle size distribution. It relies on various submodels such as nucleation, condensation, and coalescence
* The **aerosolThermo library**: a thermo package that is based on psiThermo, and contains two separate thermo libraries (which are each based on rhoThermo) for the continuous and dispersed phases. The purpose of the aerosolThermo library is to combine the continuous and dispersed thermo libraries in order to create a mixture thermo library. This startegy is following the twoPhaseMixtureThermo library of OpenFOAM's standard compressibleInterFoam solver.

## Installation

Once a copy of AeroSolved has been obtained, it can be built with:

    make

This will compile all libraries and executables. It also generates documentation in doc/output/html if [DoxyGen](https://www.stack.nl/~dimitri/doxygen) and [Graphviz](http://graphviz.org/) are installed.

The `cases` directory contains a number of cases, which can each be run by using

    ./Allrun

This prepares the case, generates a mesh, and runs aerosolEulerFoam, the main solver of AeroSolved. For some cases, the Allrun script also performs a post-processing step and generates a plot.

## Dependencies

AeroSolved has no special dependencies other than OpenFOAM. It is developed and tested against OpenFOAM-v1812 and OpenFOAM-v1906. Some cases require python3 and numpy for generating a post-processing plot. The availability of python3 is tested in the Allrun scripts.

## Documentation

From the header of this page you can navigate the model, cases, nomenclature, and classes documentation. As with standard OpenFOAM applications and libraries, further documentation is available in the headers of important `.H` files of the source code. This source code documentation is also parsed by Doxygen. We recommend generating the Doxygen documentation as explained in the Installation section.
