#!/usr/bin/python

# Script which plots a log-normal size distribution given a sectional grid
# This script is useful to find appropriate yMin and yMax

# Imports --------------------------------------------------------

import sys
import numpy as np
import matplotlib as mpl

import matplotlib.pyplot as plt

from functions import sectionalGrid
from functions import logNormal
from functions import logNormalInterval

import figStyle as fs

# Parameters -----------------------------------------------------

N = 16                      # Number of sections

yMin = 1E-24                # Lower sectional grid size [kg]
yMax = 1E-7                 # Upper sectional grid size [kg]

CMD = 1E-6                  # Count median diameter of the log-normal size distribution [m]
sigmag = 4                  # Geometric standard deviation of the log-normal size distribution

rhol = 1000                 # Dispersed mass density [kg/m^3]

gridType = 'logarithmic'    # Type of the sectional grid (linear, logarithmic)

# Script ---------------------------------------------------------

fs.prep(plt)

fig = plt.figure();

# Sectional representation

x, y = sectionalGrid(yMin, yMax, N, gridType)

CMM = 1.0/6.0*np.pi*rhol*CMD**3

dfdlogd = np.zeros(N)

for i in range(0,N):

    dfdlogd[i] = logNormalInterval(y[i], y[i+1], sigmag, CMM)

dx = np.power(x/rhol*6.0/np.pi, 1.0/3.0)

plt.plot(dx, dfdlogd, 'ok', label='Sections')

# Exact log-normal curve

dMin = CMD/sigmag**5
dMax = CMD*sigmag**5

d = np.exp(np.linspace(np.log(dMin), np.log(dMax), 1024))

dfdlogd = logNormal(d, sigmag, CMD)

plt.plot(d, dfdlogd, 'k', label='Exact')

# Useful diameters

countMeanDiameter = CMD*np.exp(0.5*np.square(np.log(sigmag)))
massMedianDiameter = CMD*np.exp(3.0*np.square(np.log(sigmag)))
massMeanDiameter = CMD*np.exp(3.5*np.square(np.log(sigmag)))

d32 = CMD*np.exp(0.5*(3.0+2.0)*np.square(np.log(sigmag)))
d53 = CMD*np.exp(0.5*(5.0+3.0)*np.square(np.log(sigmag)))

plt.plot(CMD, logNormal(CMD, sigmag, CMD), 'o', label='CMD')
plt.plot(countMeanDiameter, logNormal(countMeanDiameter, sigmag, CMD), 'o', label='$\overline{d}$')
plt.plot(d32, logNormal(d32, sigmag, CMD), 'o', label='$d_{3,2}$')
plt.plot(massMedianDiameter, logNormal(massMedianDiameter, sigmag, CMD), 'o', label='MMD')
plt.plot(massMeanDiameter, logNormal(massMeanDiameter, sigmag, CMD), 'o', label='$d_{\mathrm{mm}}$')
plt.plot(d53, logNormal(d53, sigmag, CMD), 'o', label='$d_{5,3}$')

print('Count median diameter =', CMD)
print('Count mean diameter =', countMeanDiameter)
print('Sauter mean diameter =', d32)
print('Mass median diameter =', massMedianDiameter)
print('Mass mean diameter =', massMeanDiameter)
print('Inertial diameter =', d53)

# Style and save figure

plt.xscale('log')

plt.xlabel(r'$d$ [m]')
plt.ylabel(r'$\mathrm{d}f/\mathrm{d}\log d$')

fs.post(fig, plt.legend())

plt.savefig('logNormalSectionalGrid.pdf')
