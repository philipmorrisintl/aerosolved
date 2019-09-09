#!/usr/bin/python

import sys
import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, '../../scripts')

import figStyle as fs

yMin = 1E-13
yMax = 1E-10
N = 16

rhol = 1000.0
mu = 1E-5
U = 1
R = 5E-3


# Prepare figure

fs.prep(plt)


# Load data

y = np.array([1.5806e-13, 6.4836e-13, 1.3242e-12, 2.1439e-12, 3.0855e-12, 4.1349e-12, 5.2820e-12, 6.5191e-12, 7.8400e-12, 9.2398e-12, 1.0714e-11, 1.2260e-11, 1.3873e-11, 1.5551e-11, 1.7293e-11, 1.9094e-11, 2.0955e-11]);
x = np.array([4.0321e-13, 9.8628e-13, 1.7340e-12, 2.6147e-12, 3.6102e-12, 4.7085e-12, 5.9006e-12, 7.1796e-12, 8.5399e-12, 9.9770e-12, 1.1487e-11, 1.3066e-11, 1.4712e-11, 1.6422e-11, 1.8193e-11, 2.0024e-11]);

d = np.power(x/rhol*6.0/np.pi, 1.0/3.0)

inletSectionalFluxData = np.loadtxt('postProcessing/dropletFlux/0/patch.inlet.dat');
wallsSectionalFluxData = np.loadtxt('postProcessing/dropletFlux/0/patch.walls.dat');

phiMInlet = -inletSectionalFluxData[-1,1:]
phiMWalls = wallsSectionalFluxData[-1,1:]

eta = phiMWalls/phiMInlet

St = rhol*d**2*U/(18*mu*R)


# Load reference data

dataPui = np.loadtxt('referenceData/Pui.txt')
dataChengAndWang = np.loadtxt('referenceData/ChengAndWang.txt')
dataPilou = np.loadtxt('referenceData/Pilou.txt')

StPui = dataPui[:,0]
StChengAndWang = dataChengAndWang[:,0]
StPilou = dataPilou[:,0]

etaPui = dataPui[:,1]/100;
etaChengAndWang = dataChengAndWang[:,1];
etaPilou = dataPilou[:,1];

# Plot

fig = plt.figure();

plt.plot(St, eta, '-ok', label='AeroSolved')

plt.plot(StPui, etaPui, 'o', label='Pui et al.')
plt.plot(StChengAndWang, etaChengAndWang, 'o', label='Cheng & Wang')
plt.plot(StPilou, etaPilou, 'o', label='Pilou et al.')


# Style/save

plt.xlabel(r'St')
plt.ylabel(r'$\eta$')

fs.post(fig, plt.legend())

plt.savefig('plot.pdf')
