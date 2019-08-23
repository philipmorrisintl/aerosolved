# Nomenclature

_Navigation_

1. [About](Chap1_About.md)
2. [Model](Chap2_Model.md)
3. [Tutorial](Chap3_Tutorial.md)
4. [Cases](Chap4_Cases.md)
5. [Nomenclature](Chap5_Nomenclature.md)
6. [Classes](Chap6_Classes.md)
7. [References](Chap7_References.md)

## Introduction

This document provides the nomenclature of important variables that are used inside AeroSolved. The nomenclature of the top-level solver is given first, followed by those of the aerosolModel and the aerosolThermo classes.

## General AeroSolved definitions

* `continuous`: reference to the continuous phase (i.e., the gas/vapor phase in an aerosol)
* `dispersed`: reference to the dispersed phase (i.e., the particle phase in an aerosol)
* `U`: mixture velocity $\mathbf{u}_m$ [m/s]
* `V`: relative velocity of the dispersed phase to the continuous phase [m/s] 
* `Y`: list of continuous mass fraction fields for each species [-]
* `Z`: list of dispersed mass fraction fields for each species [-]  
* `M`: mass-based number concentration of dispersed phase $M = N /\rho_m$ [#/kg]

## Top-level aerosolEulerFoam solver

* `aerosol`: reference to the aerosolModel object3
* `dpdt`: material time derivative of the pressure field [kg/m/s<sup>3</sup>]
* `fvOptions`: reference to the fv:options object
* `inertIndex`: index of the inert species in species lists
* `K`: kinetic energy field [m<sup>2</sup>/s<sup>2</sup>]
* `MRF`: IOMRFZoneList object
* `p`: pressure [kg/m/s<sup>2</sup>]
* `phi`: mixture phase flux field [kg/s]
* `phiEffY`: list of effective continuous phase flux fields, for each species [kg/s]
* `phiEffZ`: list of effective dispersed phase flux fields, for each species [kg/s]
* `psi`: compressibility field [s<sup>2</sup>/m<sup>2</sup>]
* `rho`: mixture density $\rho_m$ [kg/m<sup>3</sup>]
* `T`: mixture temperature field [K]
* `thermo`: reference to the aerosolThermo object
* `thermoCont`: reference to the rhoAerosolPhaseThermo object belonging to the continuous phase
* `thermoDisp`: reference to the rhoAerosolPhaseThermo object belonging to the dispersed phase
* `turbulence`: pointer to the compressible::turbulenceModel object


These variables have a local scope in YEqn.H, which is part of the top-level aerosolEulerFoam solver:

* `mut`: turbulent viscosity [kg/m/s]
* `mvPhi`: reference to the multi-variate convection scheme relating to flux `phi`
* `mvPhiBrownian`: reference to the multi-variate convection scheme relating to flux `phiBrownian`
* `mvPhiDrift`: reference to the multi-variate convection scheme relating to flux `phiDrift`
* `mvPhiInertial`: reference to the multi-variate convection scheme relating to flux `phiInertial`
* `phiBrownian`: reference to the Brownian (additional) drift flux field [kg/s]
* `phiDrift`: reference the corrective drift flux field, corresponding to $\mathbf{V}_d$ [kg/s]
* `phiInertial`: reference to the inertial drift flux field [kg/s]
* `Yt`: sum of all the mass fractions

## aerosolModel library

The following key variables are protected inside the aerosolModel class:

* `coalescence_`: pointer to the coalescenceModel object
* `condensation_`: pointer to the condensationModel object
* `DCont_`: continuous mass fraction diffusivity fields, for each species [m<sup>2</sup>/s]
* `DDisp_`: dispersed mass fraction diffusivity field [m<sup>2</sup>/s]
* `dMax_`: maximum allowable diameter
* `dMin_`: minimum allowable diameter
* `drift_`: pointer to the driftFluxModel object
* `mesh_`: const reference to the fvMesh object
* `mvPhi_`:  multi-variate convection scheme relating to the mixture flux
* `mvPhiBrownian_`: multi-variate convection scheme relating to the Brownian (additional) drift flux
* `mvPhiDrift_`: multi-variate convection scheme relating to the corrective flux
* `mvPhiInertial_`: multi-variate convection scheme relating to the inertial drift flux
* `nucleation_`: pointer to the nucleationModel object
* `phiBrownian_`: Brownian (additional) drift flux field [kg/s]
* `phiDrift_`: corrective drift flux field, corresponding to $\mathbf{V}_d$ [kg/s]
* `phiEff_`: list of effective particle number flux fields [#/s]
* `phiInertial_`: inertial drift flux field [kg/s]
* `tauDrift_`: mixture drift stress tensor [kg/m/s<sup>2</sup>]
* `thermo_`: aerosolThermo object
* `turbulencePtr_`: pointer to the compressibleTurbulenceModel object (set externally by the top-level solver)

## aerosolThermo library

* `activeSpecies_`: table with all active species names
* `activeSpeciesMap_`: index map of all active species inside `species_`
* `contSpeciesMap_`: index map of all continuous species inside `species_`
* `dispSpeciesMap_`: index map of all dispersed species inside `species_`
* `diffusivity_`: pointer to the mixtureDiffusivityModel object
* `inactiveSpecies_`: table with all inactive species names
* `inactiveSpeciesMap_`: index map of all inactive species inside `species_`
* `inertSpecie_`: index of the inert species inside `species_`
* `mesh_`: const reference to the fvMesh object
* `phaseMix_`: pointer to the phaseMixing object
* `species_`: table with all species names
* `thermoCont_`: the rhoAerosolPhaseThermo object belonging to the continuous phase
* `thermoDisp_`: the rhoAerosolPhaseThermo object belonging to the dispersed phase
