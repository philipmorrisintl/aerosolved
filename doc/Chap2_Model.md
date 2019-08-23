# Model

_Navigation_

1. [About](Chap1_About.md)
2. [Model](Chap2_Model.md)
3. [Tutorial](Chap3_Tutorial.md)
4. [Cases](Chap4_Cases.md)
5. [Nomenclature](Chap5_Nomenclature.md)
6. [Classes](Chap6_Classes.md)
7. [References](Chap7_References.md)

## Introduction

This page describes the aerosol model that is solved by AeroSolved. The description is given in a compact manner. For more details, the reader is referred to the PhD thesis of Edo Frederix, see [Frederix (2016)](http://dx.doi.org/10.3990/1.9789036542289). First, the governing set of equations is discussed. Second, the implementation of the solution algorithm for the governing set of equations is discussed, step by step.

## Governing set of equations

This section describes the governing set of equations. We discuss mass, momentum and energy transport. Finally, we introduce the Population Balance Equation (PBE).

### Mass transport

AeroSolved adopts the mixture model for describing the aerosol mixture. The aerosol mixture consists of a continuous phase and a dispersed phase (i.e., particles). Each phase, in turn, consists of a number of chemical species.  We let $\rho_m$ be the mixture mass density, i.e., the total mixture mass per unit volume. We defined $Y_j$ as the mass fraction of the  continuous phase of species $j$ and $Z_j$ as the mass fraction of dispersed phase of species $j$.

The continuity equation for the mixture is given by

$$
  \partial_t\rho_m + \nabla\cdot(\rho_m\mathbf{u}_m) = 0
$$

with $\mathbf{u}_m$, the mixture velocity, being defined as the velocity of the average mixture mass. The continuous mass fractions adhere to the species transport equation

$$
  \partial_t (Y_j \rho_m) + \nabla\cdot(Y_j\rho_m\mathbf{u}_m) + \nabla\cdot(Y_j\rho_m\mathbf{v}_j) = \nabla\cdot( D^{eff}_j \rho_m \nabla Y_j)  -\Gamma_j + \nabla\cdot(Y_j\rho_m\mathbf{V}_d)
$$

with $\mathbf{v}_j$ being the relative velocity between the continuous phase of species $j$ and the total continuous phase; $D_j$, the diffusion coefficient of species $j$; $D^{eff}_j$, the effective diffusion $D^{eff}_j = D_j + \nu^{turb}$;  $\Gamma_j$, an inter-phase transfer term; and $\mathbf{V}_d$, a corrective velocity which will be defined later. The dispersed mass fractions adhere to the species transport equation

$$
  \partial_t (Z_j \rho_m) + \nabla\cdot(Z_j\rho_m\mathbf{u}_m) + \nabla\cdot(Z_j\rho_m\mathbf{w}_j) = \nabla\cdot( D^{eff}_z \rho_m \nabla Z_j) +\nabla\cdot(Z_j \nabla(\rho_m D_z) ) +\Gamma_j + \nabla\cdot(Z_j\rho_m\mathbf{V}_d)
$$

with $\mathbf{w}_j$ being relative velocity between the dispersed phase of species $j$ and the total continuous phase; $D_z$, the equivalent Brownian diffusion coefficient for the dispersed phase mass fraction; and $D^{eff}_z$, the effective diffusion coefficient $D^{eff}_z = D_z + \nu^{turb}$.

The relative velocities $\mathbf{v}_j$ captures turbulent and molecular diffusion for the continuous phase and $\mathbf{w}_j$ captures  turbulent diffusion, Brownian diffusion, and inertial drift for the dispersed phase. The underlying models are commonly defined in terms of a relative motion with respect to the continuous 'carrier' phase. Therefore, $\mathbf{v}_j$ and $\mathbf{w}_j$ are also defined as relative velocities with respect to the continuous phase. However, the continuity equation is defined with respect to the mixture velocity requiring an additional corrective term to appear in the species transport equation. This term is based on the corrective velocity defined as

$$
  \mathbf{V}_d = \sum_j (Y_j \mathbf{v}_j + Z_j \mathbf{w}_j).
$$

It can be shown that the mixture continuity equation is obtained when the $Y_j$ and $Z_j$-equations are added together and summed over $j$. This guarantees that the consistency relation

$$
  \sum_j (Y_j + Z_j) = 1
$$

holds.

### Momentum transport

The mixture momentum transport equation is given by

$$
  \partial_t(\rho_m\mathbf{u}_m) + \nabla\cdot(\rho_m \mathbf{u}_m \mathbf{u}_m) = -\nabla p + \nabla\cdot(\tau_m + \tau_m^\mathrm{turb} + \tau_m^\mathrm{drift}) + \rho_m\mathbf{g},
$$

with $p$ representing pressure; $\tau_m$, mixture viscous stress tensor; $\tau_m^\mathbf{turb}$,  mixture turbulent stress tensor; $\tau_m^\mathrm{drift}$, drift stress tensor; and  $\mathbf{g}$, gravitational acceleration vector. The drift stress tensor accounts for the distribution of momentum due to the drift of vapors and droplets.

### Energy transport

The mixture energy transport equation is formulated in terms of the temperature $T$. This is more convenient when specifying boundary conditions, which are usually adiabatic or isothermal. The energy transport equation is given by:

$$
  c_v[\partial_t (\rho_m T) + \nabla\cdot(\rho_m\mathbf{u}_m T)] = \nabla\cdot(\kappa_m \nabla T) + \nabla\cdot(\kappa_m^\mathrm{turb} \nabla T) - \nabla\cdot(\mathbf{u}p) - \mathrm{D}_t K + \dot{Q},
$$

with $c_v$ representing heat capacity at constant volume; $\kappa_m$, mixture thermal conductivity; $\kappa_m^\mathrm{turb}$, mixture turbulent thermal conductivity; $K=\frac{1}{2}|\mathbf{u}|^2$, kinetic energy; $\dot{Q}$, and heat source .

### Population balance equation  

In addition to the mass fractions $Z_j$, the dispersed phase is also described by particle size distribution $n(d)$, with diameter $d$ in m. The size distribution is defined in such a way that $n(d)\mathrm{d}d$ gives the total number of particles per unit of volume with size $[d,d+\mathrm{d}d]$. The $\gamma$th length moment of the size distribution is defined as

$$
  \mathcal{N}_\gamma = \int_0^\infty d^\gamma n(d)\mathrm{d}d
$$

The zeroth moment $\mathcal{N}_0$ corresponds to the particle number concentration, in units of #/m<sup>3</sup>, and is denoted as $N$. However, in AeroSolved we usually work with $M=N/\rho_m$ (i.e., a mass-based total number concentration). The reason for this is that the moment transport equation can be written in terms of

$$
  \partial_t(\rho_m M) + \nabla\cdot(\rho_m\mathbf{u}_m M) = \ldots,
$$

in which the product $\rho_m\mathbf{u}$ is readily obtained from the solution of the momentum equation.

The size distribution $n(d)$ adheres to the Population Balance Equation (PBE):

$$
  \partial_t n(d) + \nabla\cdot[\mathbf{u}_d(d)n(d)] + \partial_d [I(d)n(d)] = J(d),
$$

with velocity $\mathbf{u}_d(d)$ representing the size-dependent droplet; $I(d)$, growth rate  (e.g., due to condensation); and $J(d)$, source term  (e.g., accounting for nucleation and coalescence). The solution to the PBE provides AeroSolved with information on the size of particles. In turn, sub-models such as the condensation and coalescence models can use this information. AeroSolved contains two models that are designed to provide an (approximate) solution to the PBE. These are:

* The **log-normal moment model**. This model solves, instead of the full PBE, the zeroth moment of the PBE. From the zeroth moment (i.e., $M$) and the mass fractions of dispersed phase $Z_j$, the particle size distribution $n(d)$ can be reconstructed by assuming that $n(d)$ has a log-normal shape with an assumed geometric standard deviation $\sigma_g$
`./libraries/aerosolModels/twoMomentLogNormal/`

* The **fixed sectional model**. This model discretizes the size space in a number of sections. For each section $i$, a scalar transport equation is solved for the total number concentration of that section $N_i$, which is defined as $N_i = \int_{d_i}^{d_{i+1}} n(d) \mathrm{d}d$, with $d_i$ being the lower bound of section $i$ and $d_{i+1}$ being the upper bound of section $i$. For the solution of the PBE subject to sectional discretization and conversational growth, the Characteristics-Based Sectional Method (CBSM) is used, see [Frederix et al. (2016)](https://doi.org/10.1016/j.jcp.2016.09.005)
`./libraries/aerosolModels/fixedSectional/`

## Aerosol sub-models

`./libraries/aerosolModels/submodels/` 

The transport equations for mass and number contain coefficients and source terms for condensation, nucleation, coalescence, and drift. These coefficients and source terms are determined by a number of sub-models that are present in AeroSolved. This section discusses these sub-models.  

### Condensation

`./libraries/aerosolModels/submodels/condensationModels/coupledCondensation/`

The PBE subject only to condensational growth is given by

$$
  \partial_t n(d) = - \partial_d [I(d)n(d)]
$$

in which $I(d)$ is the condensational growth rate in m/s, i.e., the rate of change of the diameter $d$ of a particle. This equation is solved using the Characteristics-Based Sectional Method (CBSM) inside the fixed sectional model. Moreover, for the log-normal moment model it can be shown that the condensation rate does not affect the total number concentration. The source term $\Gamma_j$ simply follows from taking the third moment of the right-hand side of the PBE.

The form of the implemented condensation model (i.e., 'coupled') is given by

$$
  I^m_j(d) = \frac{2\pi d \beta_j D_j \rho_j}{p} f(\xi_j) [\exp(\xi_j) p _j^v - p_j^\mathrm{surf}]
$$

with transition regime correction $\beta_j$, vapor diffusivity $D_j$, vapor mass density $\rho_j$, pressure $p$, vapor pressure $p_j^v$, surface vapor pressure $p_j^\mathrm{surf}$ and subscript $j$ denoting the species. Note that $I^m_j(d)$ is the condensational growth rate in units of kg/s. Moreover, we have the 'Stefan flow coupling' terms

$$
  f(\xi_j) = \frac{\xi_j}{1-\exp(\xi_j)} \qquad\mathrm{and}\qquad
  \xi_j = \frac{D_I^\star}{D_j} \log\left[\frac{p-\sum_k p_k^\mathrm{surf}}{p-\sum_k p_k^v}\right]
$$

in which $\sum_k$ goes over all active species and $D_j^I$ is the inert gas diffusivity. A characteristic feature of this model is that it couples the condensation or evaporation flux of all species together, leading to a situation in which a supersaturated species can evaporate and an undersaturated vapor can condense due to Stefan flow.


### Nucleation

`./libraries/aerosolModels/submodels/nucleationModels/coupledNucleation/`

The PBE subject only to nucleation is given by

$$
  \partial_t n(d) = J_\mathrm{nuc}
$$

The nucleation rate $J_\mathrm{nuc}$ is given by [Winkelmann et al. (2018)](https://doi.org/10.1007/s10665-017-9918-6). The mass transfer rate $\Gamma_j$ is obtained by multiplying the nucleation rate with the critical mass $s_j^\star$, with which species $j$ contributes to the critical cluster of nucleation.


### Coalescence

`./libraries/aerosolModels/submodels/coalescenceModels/`

The PBE subject to coalescence is given by

$$
  \partial_t n(d) = \frac{1}{2} \int_0^d \beta(d-\xi,\xi)n(d-\xi)n(\xi)\mathrm{d}\xi - n(d)\int_0^\infty \beta(d,\xi)n(\xi)\mathrm{d}\xi,
$$

with coalescence rate $\beta(a,b)$, giving the rate at which particles of diameters $a$ and $b$ coalesce. For the fixed sectional model, this equation is solved as described by [Kumar & Ramkrishna (1996)](https://doi.org/10.1016/0009-2509(96)88489-2). For the log-normal moment model, we use [Lee & Chen (1984)](https://doi.org/10.1080/02786828408959020). The following three coalescence kernels are implemented:

* Continuum regime coalescence (Kn < 0.1):

$$
  \beta(a,b) = K(a+b)(1/a + 1/b)
$$

* Gas-slip regime (0.1 < Kn < 1):

$$
  \beta(a,b) = K(a+b)(1/a+1/b + 2A\lambda/a^2 + 2A\lambda/b^2)
$$

* Free-molecule regime (Kn > 10):

$$
  \beta(a,b) = \tilde{K} (a+b)^2 (1/a^{3/2} + 1/b^{3/2})
$$

In these equations $Kn$ is the Knudsen number; $K$ and $\tilde{K}$ are coalescence coefficients $K=2k_B T/3\mu$ and $\tilde{K} = b3\sqrt{3}/2\sqrt{\mu^2/(\rho_\ell k_B T)}K$  with the Boltzmann's constant $k_B$,  temperature $T$,  gas viscosity $\mu$, liquid mass density $\rho_\ell$,  mean free path $\lambda$, coefficients $A=1.591$, and $b$ (see [Lee & Chen (1984)](https://doi.org/10.1080/02786828408959020)).

In order to capture the transition from one regime to the other, also a harmonic mean of the gas-slip and free-molecule coalescence kernel is implemented.

### Inertial drift

Inertial drift is captured inside the continuum--particle relative velocity $\mathbf{w}_j(d)$. Since particles are assumed to be internally mixed (i.e., particle composition is independent of size), all species have, given a particle size $d$, the same inertial drift velocity. This implies that $\mathbf{w}_j(d)=\mathbf{w}(d)$. There are two inertial drift models implemented in AeroSolved, which will be discussed now.

#### Algebraic drift model

`./libraries/aerosolModels/submodels/driftFluxModel/inertialModels/Manninen/`

One can adopt the so-called 'local equilibrium assumption', see Manninen et al. (1996). In that case, the inertial drift velocity can be modeled algebraically, by

$$
  \mathbf{w} = \tau [(1-\gamma)\mathbf{g}-\mathrm{d}\mathbf{u}_m/\mathrm{d}t],
$$

with $\gamma=\rho_c/\rho_d$, gravitational acceleration vector $\mathbf{g}$ and Stokes particle relaxation time
$$
  \tau = \frac{\rho_d d^2}{18\mu}.
$$
Note that the algebraic drift model is only valid for small particle Stokes numbers. It was shown that the algebraic model fails for large particle Stokes numbers  [Frederix (2016)](http://dx.doi.org/10.3990/1.9789036542289).


#### Full Stokes drift model

`./libraries/aerosolModels/submodels/driftFluxModel/inertialModels/fullStokes/`

The local equilibrium assumption allows to reduce the full PDE for the particulate velocity to an algebraic formula. If this assumption is not valid, the PDE should be solved by

$$
  \partial_t\mathbf{w} + \partial_t\mathbf{u}_m + (\mathbf{w}+\mathbf{u}_m)\cdot\nabla(\mathbf{w}+\mathbf{u}_m) = -\frac{1}{\tau} \mathbf{w} + (1-\gamma)\mathbf{g}
$$

This equation can be solved for $\mathbf{w}$ subject to appropriate boundary conditions.

### Brownian drift

`./libraries/aerosolModels/submodels/driftFluxModel/BrownianModels/StokesEinstein/`

The particulate phase may diffuse because of Brownian motion. This is modeled inside AeroSolved by using an additional diffusion term. The Brownian diffusivity is computed by using the Stokes-Einstein relationship, which is given by

$$
  D(d) = \frac{k_B T C_c}{3\pi\mu d}
$$

with Cunningham correction factor $C_c$.

## AeroSolved Boundary Conditions 

### sectionalSubGridDepositionVelocity

`libraries/aerosolModels/fixedSectional/derivedFvPatchFields/sectionalSubGridDepositionVelocity/`

This boundary condition uses the analytical solution to one-dimensional particle motion subject to Stokes drag to compute the particle velocity at the wall ([Frederix (2018)](https://doi.org/10.1016/j.compfluid.2017.09.018)). In scaled form, the particle equation of motion is given by
$$
\ddot{x}(t) + \dot{x}(t) + ux(t) = g
$$

The particle velocity at the wall is then given by $v = \dot{x}(t^*)$, with $t^*$ representing the impaction time at the wall of the particle. 

## The aerosolEulerFoam solver

`applications/solvers/aerosolEulerFoam/`

The aerosolEulerFoam solver is based on the standard reactingFoam OpenFOAM solver. It can be seen as a reactingFoam solver in which the reaction model is replaced by an aerosol model. The solver is a transient one and solves the evolution of a two-phase aerosol mixture by using OpenFOAM's PIMPLE algorithm. The top-level structure of the solver is as follows.

For each time step, do:

* Correct the aerosol model by using `aerosol->correct()`. This aerosolModel member function has the task of updating all aerosol-related coefficients and fields. These are: the molecular diffusivities; inertial drift velocity; Brownian diffusivity; corrective drift velocity $\mathbf{V}_d$; nucleation source term; condensation source term; and coalescence source term. The twoMomentLogNormalAnalytical and fixedSectional models solve the 'right-hand side' contributions of the $Y_j$ and $Z_j$ equations as captured by $\Gamma_j$ themselves. Therefore, these two aerosol models do not produce a nucleation, condensation, or coalescence source term. The 'twoMomentLogNormal' model, on the other hand, does:

* Solve the continuity equation. Standard reactingFoam

* For each outer iteration, do:

    * Solve the momentum equation. Standard reactingFoam.
  
    * Solve the pre-species part of the aerosol model by using `aerosol->solvePre()`. This aerosolModel member function is directly implemented by the specific aerosolModel. Its implementation therefore fully depends on the model. Although the current models have no implementation,  future models may.
  
    * Solve the species equations. For each species $j$ and for each phase (continuos or dispersed) the corresponding mass transport equation is solved. This happens in a similar way as that with which the species equations in the reactingFoam solver are solved.
  
    * Solve the post-species part of the aerosol model using `aerosol->solvePost()`. This aerosolModel member function is also directly implemented by the specific aerosolModel. For the twoMomentLogNormalAnalytical model, the $M$-equation is solved and the 'right-hand side' contributions of the $Y_j$ and $Z_j$ equations are solved analytically in an additional fractional step. For the fixedSectional model, all $M_i$ transport equations without the 'right-hand side' contributions are solved first. Next, the 'right-hand side' contributions of the $M_i$, $Y_j$ and $Z_j$ equations (i.e., nucleation, condensation and coalescence) are solved analytically in an additional fractional step.
  
    * For each inner iteration, do:

        * Solve the pressure equation. Standard reactingFoam

This completes the short description of the model and the corresponding implementation of the model in the top-level solver.
