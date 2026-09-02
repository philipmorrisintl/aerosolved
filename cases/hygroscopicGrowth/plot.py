#!/usr/bin/env python3

import sys
import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, '../../scripts')

import figStyle as fs

# Prepare figure

fs.prep(plt)


# Load data

masses = np.loadtxt('postProcessing/masses/0/volFieldValue.dat')
totalMass = np.loadtxt('postProcessing/totalMass/0/volFieldValue.dat')
diameters = np.loadtxt('postProcessing/meanDiameters/0/volFieldValue.dat')
fluxes = np.loadtxt('postProcessing/massFlux/0/patch.walls.dat')

t = masses[:,0]


# Plot

fig = plt.figure('mass')

plt.plot(t, masses[:,1]/totalMass[:,1], color='C0', label='PG')
plt.plot(t, masses[:,2]/totalMass[:,1], color='C1', label='VG')
plt.plot(t, masses[:,3]/totalMass[:,1], color='C2', label='Water')

plt.plot(t, masses[:,4]/totalMass[:,1], '--', color='C0')
plt.plot(t, masses[:,5]/totalMass[:,1], '--', color='C1')
plt.plot(t, masses[:,6]/totalMass[:,1], '--', color='C2')

plt.plot([], [], '-k', label='continuous')
plt.plot([], [], '--k', label='dispersed')

##

fig = plt.figure('d')

plt.plot(t, diameters[:,1]/1e-6, color='C0', label='dcm')
plt.plot(t, diameters[:,2]/1e-6, color='C1', label='CMD')

##

fig = plt.figure('fluxCont')

plt.plot(t, fluxes[:,1]*1e9, label='PG')
plt.plot(t, fluxes[:,2]*1e9, label='VG')
plt.plot(t, fluxes[:,3]*1e9, label='Water')
plt.plot(t, fluxes[:,4]*1e9, label='Air')

fig = plt.figure('fluxDisp')

plt.plot(t, fluxes[:,5]*1e9, label='PG')
plt.plot(t, fluxes[:,6]*1e9, label='VG')
plt.plot(t, fluxes[:,7]*1e9, label='Water')

# Style/save

fig = plt.figure('mass')

plt.yscale('log')

plt.xlabel(r'$t$ [s]')
plt.ylabel(r'mass fraction')

fs.post(fig, plt.legend())

plt.savefig('mass.pdf')

##

fig = plt.figure('d')

plt.xlabel(r'$t$ [s]')
plt.ylabel(r'scaled volume-averaged diameter')

fs.post(fig, plt.legend())

plt.savefig('d.pdf')

##

fig = plt.figure('fluxCont')

plt.xlabel(r'$t$ [s]')
plt.ylabel(r'continuous mass flux [$\mu$g/s]')

fs.post(fig, plt.legend())

plt.savefig('fluxCont.pdf')

##

fig = plt.figure('fluxDisp')

plt.xlabel(r'$t$ [s]')
plt.ylabel(r'dispersed mass flux [$\mu$g/s]')

fs.post(fig, plt.legend())

plt.savefig('fluxDisp.pdf')
