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

t = masses[:,0]


# Plot

fig = plt.figure()

plt.plot(t, masses[:,1]/totalMass[:,1], color='C0', label='Water1')
plt.plot(t, masses[:,2]/totalMass[:,1], color='C1', label='Water2')
plt.plot(t, masses[:,3]/totalMass[:,1], color='C2', label='Water3')

plt.plot(t, masses[:,4]/totalMass[:,1], '--', color='C0')
plt.plot(t, masses[:,5]/totalMass[:,1], '--', color='C1')
plt.plot(t, masses[:,6]/totalMass[:,1], '--', color='C2')

plt.plot([], [], '-k', label='continuous')
plt.plot([], [], '--k', label='dispersed')

# Style/save

plt.yscale('log')

plt.xlabel(r'$t$ [s]')
plt.ylabel(r'mass fraction')

fs.post(fig, plt.legend())

plt.savefig('plot.pdf')
