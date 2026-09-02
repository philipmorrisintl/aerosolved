#!/usr/bin/env python3

import sys
import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, '../../scripts')

import figStyle as fs

# ----------------------------------------------------------------------

# Sectional data

yMin = 1E-22
yMax = 1E-12
N = 16

# Mass density of droplets and air [kg/m3]

rhol = 1000.0
rhog = 1.0

# Dynamic viscosity of air [kg/m/s]

mug = 1.84E-5

# Mole weight of air [g/mol]

Mg = 28.96

# Radius of the deposition plate [m]

R = 4E-3

# Exposure chamber height [m]

gap = 2.9E-3

# Mixture mass flow rate [kg/s]

massFlow = 1.26067454E-09

# Wedge angle [degrees]

angle = 5

# Temperature [K]

T = 293.15

# Pressure [Pa]

p = 1.0E5

# ----------------------------------------------------------------------


# Prepare figure

fs.prep(plt)

fig = plt.figure()


# Plot data

logy = np.linspace(np.log10(yMin), np.log10(yMax), N+1)
x = np.power(10, (logy[1:]+logy[:-1])/2.0)

d = np.power(x/rhol*6.0/np.pi, 1.0/3.0)

inletSectionalFluxData = np.loadtxt('postProcessing/dropletFlux/0/patch.inlet.dat')
lowerSectionalFluxData = np.loadtxt('postProcessing/dropletFlux/0/patch.lower.dat')

phiInlet = -inletSectionalFluxData[:,1:]
phiLower = lowerSectionalFluxData[:,1:]

eta = phiLower/phiInlet

plt.plot(d, eta[-1,:], '-ok', label='AeroSolved')


# Plot theoretical solution

dMin = 7e-9
dMax = 1e-5

d = np.exp(np.linspace(np.log(dMin), np.log(dMax), 1024))

area = np.pi*R**2*angle/360
flow = massFlow/rhog

tau = rhol*d**2/(18.0*mug)
G = 9.81
V = G*tau

kB = 1.380649E-23
NA = 6.02214E23

mg = 0.001*Mg/NA
lamb = np.sqrt(8.0*kB*T/(np.pi*mg)) * 4.0/5.0*mug/p
Kn = lamb/d
C = 1.0 + Kn*(2.34+1.05*np.exp(-0.39/Kn))

diff = kB*T*C/(3.0*np.pi*mug*d)

nug = mug/rhog
x = 2.0/3.0*R
u = flow/(2.0*np.pi*angle/360.0*x*gap)
Pr = nug/diff
delta = 1.0/Pr**(1.0/3.0) * np.sqrt(nug*x/u)

etaTheory = area/flow * (V + diff/delta)
etaTheory[etaTheory>1] = 1

plt.plot(d, etaTheory, '-', color='C0', label='Theory (Lucci et al. (2022))')


# Style/save

plt.xlabel(r'$d$')
plt.ylabel(r'$\eta$')

plt.xscale('log')
plt.yscale('log')

fs.post(fig, plt.legend(loc='best'))

plt.savefig('plot.pdf')
