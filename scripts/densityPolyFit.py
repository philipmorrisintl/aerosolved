#!/usr/bin/python

# Script which plots the PG, VG and water temperature-dependent
# density functions, and the polynomial which is used to
# approximate them

# Imports --------------------------------------------------------

import sys
import numpy as np
import matplotlib as mpl
import thermophysicalFunctions as tps

import matplotlib.pyplot as plt

import figStyle as fs

# Parameters -----------------------------------------------------

T1 = 273.15     # Lowest temperature
T2 = 523.15     # Highest temperature

Np = 3          # Polynormial degree

N = 1024        # Number of points to plot

# Script ---------------------------------------------------------

fs.prep(plt)

fig = plt.figure();

T = np.linspace(T1, T2, N, endpoint=True)

rhoPG = np.zeros(N)
rhoVG = np.zeros(N)
rhoWater = np.zeros(N)

for i in range(0,N):

    rhoPG[i] = tps.rhoPGl(T[i])
    rhoVG[i]= tps.rhoVGl(T[i])
    rhoWater[i] = tps.rhoWaterl(T[i])

pPG = np.polyfit(T, rhoPG, Np)
pVG = np.polyfit(T, rhoVG, Np)
pWater = np.polyfit(T, rhoWater, Np)

rhoPGp = np.polyval(pPG, T)
rhoVGp = np.polyval(pVG, T)
rhoWaterp = np.polyval(pWater, T)

# Print

print("polyfit(rhoPG) =", pPG[::-1])
print("polyfit(rhoVG) =", pVG[::-1])
print("polyfit(rhoWater) =", pWater[::-1])

# Plot

plt.plot(T, rhoPG, label='PG')
plt.plot(T, rhoVG, label='VG')
plt.plot(T, rhoWater, label='Water')

plt.plot(T, rhoPGp, 'o', label='PG polyfit', markevery=int(N/32), mew=0.5, mec='white')
plt.plot(T, rhoVGp, 'o', label='VG polyfit', markevery=int(N/32), mew=0.5, mec='white')
plt.plot(T, rhoWaterp, 'o', label='Water polyfit', markevery=int(N/32), mew=0.5, mec='white')

# Style and save figure

plt.xlabel(r'$T$ [K]')
plt.ylabel(r'$\rho_\ell$ [kg/m$^3$]')

fs.post(fig, plt.legend())

plt.savefig('densityPolyFit.pdf')
