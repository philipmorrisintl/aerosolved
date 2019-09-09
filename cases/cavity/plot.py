#!/usr/bin/python

import sys
import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, '../../scripts')

import figStyle as fs

time = '10'


# Prepare figure

fs.prep(plt)


# Load Ghia's data

GhiaU = np.loadtxt('referenceData/GhiaU.txt')
GhiaV = np.loadtxt('referenceData/GhiaV.txt')


# Load sampled data

SimU = np.loadtxt('postProcessing/sampleDict/' + time + '/Y_U.xy')
SimV = np.loadtxt('postProcessing/sampleDict/' + time + '/X_U.xy')


# Plot

fig = plt.figure('plotU')

plt.plot(SimU[:,0], SimU[:,1], label='AeroSolved')
plt.plot(GhiaU[:,0], GhiaU[:,1], 'ok', label='Ghia et al.')

fig = plt.figure('plotV')

plt.plot(SimV[:,0], SimV[:,2], label='AeroSolved')
plt.plot(GhiaV[:,0], GhiaV[:,1], 'ok', label='Ghia et al.')


# Style/save

fig = plt.figure('plotU')

plt.xlabel('$y$')
plt.ylabel('$u$')

fs.post(fig, plt.legend())

plt.savefig('plotU.pdf')

fig = plt.figure('plotV')

plt.xlabel('$x$')
plt.ylabel('$v$')

fs.post(fig, plt.legend())

plt.savefig('plotV.pdf')
