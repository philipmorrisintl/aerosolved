#!/usr/bin/env python3

# Script which computes a polynomial which approximates butanol
# density in a given temperature range

# Imports --------------------------------------------------------

import sys
import numpy as np
import matplotlib as mpl
import thermophysicalFunctions as tps

import matplotlib.pyplot as plt

import figStyle as fs

# Parameters -----------------------------------------------------

T1 = 270        # Lowest temperature
T2 = 300        # Highest temperature

Np = 2          # Polynormial degree

N = 1024        # Number of points to plot

# Script ---------------------------------------------------------

fs.prep(plt)

fig = plt.figure();

T = np.linspace(T1, T2, N, endpoint=True)

muButanol = np.zeros(N)

for i in range(0,N):

    muButanol[i] = tps.muButanol(T[i])

pButanol = np.polyfit(T, muButanol, Np)

muButanolp = np.polyval(pButanol, T)

# Print

print("polyfit(muButanol) =", pButanol[::-1])

# Plot

plt.plot(T, muButanol, label='Butanol')

plt.plot(T, muButanolp, 'o', label='Butanol polyfit', markevery=int(N/32), mew=0.5, mec='white')

# Style and save figure

plt.xlabel(r'$T$ [K]')
plt.ylabel(r'$\mu$ [kg/m/s]')

fs.post(fig, plt.legend())

plt.savefig('butanolViscosityFit.pdf')
