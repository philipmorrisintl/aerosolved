#!/usr/bin/env python3

# Script to compute mass fractions of a vapor mixture at saturation S with inert
# species molar composition g, held over a hypothetical pool with liquid molar
# composition f

# Imports ----------------------------------------------------------------------

import numpy as np
import os
import thermophysicalFunctions as tpf

# Parameters -------------------------------------------------------------------

T = 310.15                  # Temperature [K]
p = 1e5                     # Pressure [Pa]

S = 0.99                    # The saturation of the mixture

# Order of species: PG, VG, water, air

f = [0.2, 0.2, 0.6, 0]      # Pool phase mixture mole fractions
g = [0, 0, 0, 1]            # Inert mole fractions
gamma = [10, 10, 1, 1]      # Activity coefficients

# Specie groups

A = [0,1,2]                 # Active
I = [3]                     # Inert

# Script -----------------------------------------------------------------------

# Specific densities

rhog = np.array([tpf.rhoPGg(T), tpf.rhoVGg(T), tpf.rhoWaterg(T), tpf.rhoAirg(T)])

# Molar weights

M = np.array([tpf.MPG, tpf.MVG, tpf.MWater, tpf.MAir])

# Saturation pressures

pSat = np.array([tpf.psPG(T), tpf.psVG(T), tpf.psWater(T), tpf.psAir(T)])

# Convert to numpy arrays

f = np.array(f)
g = np.array(g)
gamma = np.array(gamma)

# Mole fractions using Raoult's law

X = S*f*gamma*pSat/p

# Inert species mole fractions

X[I] = (1.0-np.sum(X[A]))*g[I]

# Mixture mean molar mass

Mm = np.sum(X*M)

# Mass fractions

Y = X*M/Mm

# Volume fractions

alpha = Y/rhog/np.sum(Y/rhog)

# Print

print("Mole fractions X = ", X)
print("Mass fractions Y = ", Y)
print("Volume fractions alpha = ", alpha)
print("sum(X) =", np.sum(X))
print("sum(Y) =", np.sum(Y))
print("sum(alpha) =", np.sum(alpha))
