#!/usr/bin/env python3

# Script to compute saturated mixture flow

# Imports --------------------------------------------------------

import numpy as np
import os
import thermophysicalFunctions as tpf

# Parameters -----------------------------------------------------

T = 293.15                  # Temperature [K]
p = 1e5                     # Pressure [Pa]

quality = 0.01              # Quality of the mixture (dispersed
                            # mass flow over total mass flow)

# Order of species: PG, VG, water, air

alphat = [0.4, 0.4, 0.2, 0] # Precursor liquid volume fractions
betat = [0, 0, 0, 1]        # Dilution gas volume fractions
gamma = [1, 1, 1, 1]        # Activity coefficients

# Specie groups

A = [0,1,2]
I = [3]

# Iterative parameters

tolerance = 1e-12
maxIter = 100

# Script ---------------------------------------------------------

# Specific densities

rhog = np.array([tpf.rhoPGg(T), tpf.rhoVGg(T), tpf.rhoWaterg(T), tpf.rhoAirg(T)])
rhol = np.array([tpf.rhoPGl(T), tpf.rhoVGl(T), tpf.rhoWaterl(T), tpf.rhoAirl(T)])

# Molar weights

M = np.array([tpf.MPG, tpf.MVG, tpf.MWater, tpf.MAir])

# Mass fractions w.r.t. initial gas mixture

ym = np.array(betat)*rhog
yt = ym/np.sum(ym)
xt = yt/M/np.sum(yt/M)

# Mass and mole factions w.r.t. initial liquid mixture

zm = np.array(alphat)*rhol
zt = zm/np.sum(zm)
wt = zt/M/np.sum(zt/M)

# Mean molecular weights

Mct = 1/np.sum(yt/M)
Mdt = 1/np.sum(zt/M)

# Quality in mole flow

moleQuality = quality/Mdt / ((1.0-quality)/Mct + quality/Mdt)

# Mass and mole fractions w.r.t. the mixture

Xt = xt*(1.0-moleQuality)
Wt = wt*moleQuality

# Saturation pressures

pSat = np.array([tpf.psPG(T), tpf.psVG(T), tpf.psWater(T), tpf.psAir(T)])

# Saturation 'matrix'

S = gamma*pSat/p

# Equilibrium search

X = Xt
W = Wt

Mt = Xt+Wt

count = 0

while True:

    count = count+1

    X1 = np.sum(X)
    W1 = np.sum(W)

    # Update only active species

    dX = Mt[A]*S[A]*X1/(W1+X1*S[A]) - X[A]

    # Check

    if count > maxIter:

        print("Did not converge")

        break

    elif np.sum(np.abs(dX)) < tolerance:

        # Mixture mass fractions

        Y = X*M/np.sum(X*M+W*M)
        Z = W*M/np.sum(X*M+W*M)

        # Mixture volume fractions

        Rhog = np.maximum(rhog,1e-8)
        Rhol = np.maximum(rhol,1e-8)

        alpha = Y/Rhog/np.sum(Y/Rhog + Z/Rhol)
        beta = Z/Rhol/np.sum(Y/Rhog + Z/Rhol)

        print("Converged in", count, "iterations, with solution")
        print("    Continuous phase mole fractions X = ", X)
        print("    Dispersed phase mole fractions W = ", W)
        print("    Continuous phase mass fractions Y = ", Y)
        print("    Dispersed phase mass fractions Z = ", Z)
        print("    Continuous phase volume fractions alpha = ", alpha)
        print("    Dispersed phase volume fractions beta = ", beta)

        print("    sum(X+W) =", np.sum(X+W))
        print("    sum(Y+Z) =", np.sum(Y+Z))
        print("    sum(alpha+beta) =", np.sum(alpha+beta))

        break

    else:

        # Update

        X[A] = X[A] + dX
        W[A] = W[A] - dX
