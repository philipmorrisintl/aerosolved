# Classes

_Navigation_

1. [About](Chap1_About.md)
2. [Model](Chap2_Model.md)
3. [Tutorial](Chap3_Tutorial.md)
4. [Cases](Chap4_Cases.md)
5. [Nomenclature](Chap5_Nomenclature.md)
6. [Classes](Chap6_Classes.md)
7. [References](Chap7_References.md)

## Introduction

This page gives a brief overview of all the classes which are present in AeroSolved. The overview is hierarchal, following the directory structure of AeroSolved

## The aerosolModel class

This class contains the implementation of different aerosol models, which are responsible for the modeling of the Population Balance Equation (PBE). Ultimately, the goal of solving the PBE (either in its true form, or in terms of its moments) is to produce size information on the distribution of particles, locally. In turn, this size information provides sub-models such as the condensationModel and the driftFluxModel with relevant average diameters. The aerosolModel class is the most important class in AeroSolved, and also contains the aerosolThermo class, which describes the thermodynamical state of the two-phase (i.e., continuous phase and dispersed phase) mixture.

The following key member functions are contained in the aerosolModel class:

* `correct()`: corrects coefficients, rates and velocities which are part of the aerosol model
* `correctModel()`: called by `correct()`, implemented by the selected aerosol model and responsible for the correction of coefficients, rates and velocities which are specific to the selected aerosol model
* `meanDiameter(p,q)`: implemented by the selected aerosol model and provides the mean diameter which is related to moment $p$ and moment $q$, and which carries the definition
$$
    d_{p,q} = \left(\frac{\int_0^\infty d^p n(d)\mathrm{d}d}{\int_0^\infty d^q n(d)\mathrm{d}d}\right)^{1/(p-q)}
$$
* `medianDiameter(p)`: implemented by the selected aerosol model and provides the median diameter which is related to moment $p$. For example, with $p$ = 0, the cound median diameter (CMD) is provided, while with $p$ = 3, the mass (or volume) median diameter (MMD) is provided
* `Qdot()`: implemented by the selected aerosol model and provides the heat release rate associated with aerosol-related processes
* `R(Y)`: implemented by the selected aerosol model and provides the the right-hand side aerosol-related source term for species mass fraction `Y`
* `solvePost()`: implemented by the selected aerosol model and responsible for executing the part of the solution algorithm to the PBE which should be run after the solution of the mass transport equations
* `solvePre()`: implemented by the selected aerosol model and responsible for executing the part of the solution algorithm to the PBE which should be run before the solution of the mass transport equations
* `thermo()`: returns a (const) reference to the aerosolThermo object

### Models

The actual implementation of the aerosolModel is specified inside the aerosolProperties file. Models can be chosen from the following list:

* **twoMomentLogNormal**. Solves the PBE by assuming a log-normal distribution the width of which is fixed. The distribution is closed by solving the number concentration transport equation. Explicit right-hand side source terms are provided for the $Y_j$ and $Z_j$-equations
* **twoMomentLogNormalAnalytical**. A copy of the twoMomentLogNormal model, but does not provide explicit source terms for the $Y_j$ and $Z_j$ equations because these source terms are solved analytically in the `solvePost()` step. Generally, the twoMomentLogNormalAnalytical moment is more stable than the twoMomentLogNormal model and, therefore, is recommended for use
* **fixedSectional**. Solves the PBE by using a sectional discretization, in which the sections (specified in terms of particle mass) are fixed in time and space. The fixedSectional object relies on the fixedSectionalSystem, which, in turn, provides the sectional distribution and interpolation functionalities
* **noAerosol** (can be selected with 'none'). Provides an empty implementation of the aerosolModel class

### Sub-models

The following sub-models are implemented in the aerosolModel class:

* **coalescenceModel**. Provides models for the computation of the coalescence kernel. The coalescence kernel is provided in polynomial form, via the `rate(...)` member function which returns a `coaData` object 
* **condensationModel**. Provides models for the computation of the condensation rate. The condensation rate source and sink coefficients are returned via the `rate(...)` member function using the `conData` object
* **nucleationModel**. Provides models for the computation of the nucleation rate. The nucleation rate, size and composition are returned via the `rate(...)` member function using the `nucData` object
* **driftFluxModel**. Is a simple container classes which computes the corrective drift flux and drift stress tensor, based on the following 'sub-sub-models':
    - _BrownianModel_. Provides the Brownian diffusivity given a droplet size
    - _diffusionModel_. Provides the vapor diffusivity for a given species index $j$
    - _inertialModel_. Provides the inertial drift velocity given a size $d$. A size name is also provided, such that the drift velocity field can be stored and re-used for more advanced non-algebraic models

### functionObjects

The aerosolModel class provides the following functionObjects, which can be configured inside controlDict:

* **countMeanDiameter**: computes the count mean diameter field
* **Knudsen number**: computes the Knudsen number given the local particle size
* **massFlux**: evaluates all the mass fluxes, per phase, per species and per provided patch or faceZoneSet
* **meanDiameter**: calls the `meanDiameter(p,q)` member function of the aerosolModel
* **medianDiameter**: calls the `medianDiameter(p)` member function of the aerosolModel
* **Qdot**: calls the `Qdot()` member function of the aerosolModel
* **sectionalFlux**: evaluates all the number fluxes per provided patch or faceZoneSet
* **twoMomentFlux**: evaluates the total number flux per provided patch or faceZoneSet

### derivedFvPatchFields

The aerosolModel class provides the following derivedFvPatchFields, which can be used as boundary conditions:

* **zeroGradientAbsorbingWall**: Provides a zero-gradient condition when evaluating the value and a zero-value condition when evaluating the surface-normal gradient. This is useful for capturing both inertial drift and diffusion deposition across a patch
* **zeroGradientDepositionVelocity**: Provides a zero-gradient condition when the direction is wall-ward, and a zero-value condition when the direction if the velocity would be wall-outward
* **saturatedMixture**: Provides a fixed-value boundary condition specifying a vapor mixture which is at thermodynamic equilibrium
* **sectionalConstant**: Provides a boundary condition which specifies a section-independent number density while respecting the total dispersed mass fraction at the patch, for each $M_i$ field
* **sectionalLogNormal**: provides a boundary condition for the sectional number concentration fields $M_i$, such that a log-normal size distribution is specified while respecting the total dispersed mass fraction at the patch
* **sectionalSubGridDepositionVelocity**: uses the analytical solution to the one-dimensional particle equation of motion to predict deposition velocity at the patch
* **twoMomentLogNormal**: provides a boundary condition for the total number concentration field, such that a desired log-normal distribution is specified

## The aerosolThermo class

This class provides the thermodynamic state of the two-phase continuous-dispersed mixture. It does so, by combining two thermo objects: one for the continuous phase and another for the dispersed phase. These two thermo objects are mixed together using phase mixing laws which can be specified in a case for each mixture property. The aerosolThermo class is derived from the psiThermo class, which is a standard OpenFOAM library. This means that the aerosolThermo class works together with other standard OpenFOAM tools, such as turbulence models.

The following key member functions are contained in the aerosolThermo class:

* `correct()`: updates all mixture properties
* `correctThermo()`: updates the contained thermo objects (for each phase)
* `Cp()`: computes and returns the mixture heat capacity at constant pressure field [J/K]
* `Cv()`: computes and returns the mixture heat capacity at constant volume field [J/K]
* `nu()`: computes and returns the mixture kinematic viscosity field [m<sup>2</sup>/s]
* `rho()`: computes and returns the mixture density field [kg/m<sup>3</sup>]
* `sumY()`: computes and returns the sum of continuous mass fraction fields
* `sumZ()`: computes and returns the sum of dispersed mass fraction fields
* `thermoCont()`: returns a (const) reference to the continuous thermo object
* `thermoDisp()`: returns a (const) reference to the dispersed thermo object
* `W()`": computes and returns the mixture molecular weight [kmol/kg]
* `WMix()`": alias of `W()`
* `Y()`": returns a (const) reference to the list of continuous mass fraction fields
* `Z()`": returns a (const) reference to the list of dispersed mass fraction fields

### rhoAerosolPhaseThermo

Each phase thermo object (for the continuous and the dispersed phase) is of type rhoAerosolPhaseThermo. This class is derived from the rhoReactionThermo class, which is standard and OpenFOAM and which, in turn, derives from the rhoThermo class. The rhoAerosolPhaseThermo class contains the `WMix()` member function, which is pure virtual. It is implemented by the heAerosolRhoThermo class which is also part of aerosolThermo. This means that the rhoAerosolPhaseThermo can only be used as a template to heAerosolRhoThermo, and that, therefore, these two classes go hand in hand. In any case, the fact that rhoAerosolPhaseThermo derives from the standard rhoThermo class allows the phase thermo objects to be used by other standard OpenFOAM tools.

### The icoFunction equationOfState class

In addition to the standard equationOfState classes from OpenFOAM, aerosolThermo implements the icoFunction equationOfState. This allows the equationOfState to be specified as a Function1 object. This makes the rhoAerosolPhaseThermo object rather versatile, however, the call to the Function1 object for each cell can be quite expensive

### Sub-models

The following sub-models are implemented in the aerosolThermo class:

* **mixtureDiffusivityModel**. Provides the effective diffusivity of a species $j$ in the current mixture
* **phaseMixingModel**. Provides the mixing laws for the mixture viscosity, conductivity and heat capacities. The mixing laws are specified in the thermophysicalProperties case file. The implemented laws are: 'continuousPhase', 'dispersedPhase', 'harmonic', 'mass', 'mole' and 'volume'.

### functionObjects

The aerosolThermo class provides the following functionObjects, which can be configured inside controlDict:

* **thermoField**: computes a thermodynamic field. The following fields can be selected: 'Cp', 'Cv', 'gamma', 'Cpv', 'CpByCpv', 'hc', 'nu', 'kappa', 'sumY', 'sumZ' and 'lambda' (mean free path)

## The customFunctions class

This class implements custom Function1 functions. The currently implemented functions are: exponential, NSRDS and VDI. They can be used to provide parameter input in the thermophysicalProperties.continues/dispersed files.

## The customTurbulenceModels class

This class provides custom LES and RANS models, by extending the standard turbulenceModel object. Currently, the Vreman LES model is implemented.
