#!/usr/bin/python

import sys
import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, '../../../scripts')

import figStyle as fs


# Prepare figure

fs.prep(plt)


# Read data

data = np.loadtxt('postProcessing/probes/0/dcm')
refData = np.loadtxt('referenceData/Ranz.txt')


# Plot

fig = plt.figure()

plt.plot(data[:,0], np.square(data[:,1])*1E6, label='AeroSolved')
plt.plot(refData[:,0], refData[:,1], 'ok', label='Ranz et al.')


# Style/save

plt.xlabel(r'$t$ [s]')
plt.ylabel(r'$d^2$ [mm${}^2$]')

fs.post(fig, plt.legend())

plt.savefig('plot.pdf')
