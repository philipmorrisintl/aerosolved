#!/usr/bin/python

import sys
import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, '../../scripts')

import figStyle as fs

yMin = 2E-15
yMax = 2E-11
rhol = 1000.0
U = 1.0
mu = 1.81E-5
D = 2E-3


# Prepare figure

fs.prep(plt)


# Belyaev & Levin model

def GBL(R):
    return (2.0 + 0.62/R)

def aBL(R, St):
    return (1 - 1/(1+GBL(R)*St))

def ABLf(R, St):
    return (1 + aBL(R,St)*(R-1))


# Load data

R = np.loadtxt('R.txt')

inletFluxData = np.loadtxt('postProcessing/massFlux/0/patch.inlet.dat');
probeFluxData = np.loadtxt('postProcessing/massFlux/0/faceZoneSet.probe.dat');

inletSectionalFluxData = np.loadtxt('postProcessing/sectionalFlux/0/patch.inlet.dat');
probeSectionalFluxData = np.loadtxt('postProcessing/sectionalFlux/0/faceZoneSet.probe.dat');

phiInlet = np.sum(-inletFluxData[-1,1:])
phiProbe = np.sum(probeFluxData[-1,1:])

phiMInlet = -inletSectionalFluxData[-1,1:]
phiMProbe = probeSectionalFluxData[-1,1:]

M0 = phiMInlet/phiInlet;

A = phiMProbe/(phiProbe*M0)

N = len(A)

y = np.power(10, np.linspace(np.log10(yMin), np.log10(yMax), N+1))
x = np.exp((np.log(y[:-1])+np.log(y[1:]))/2.0)

d = np.power(x/rhol*6.0/np.pi, 1.0/3.0)
tau = np.power(d, 2.0)*rhol/(18.0*mu)

St = tau*U/D

StBL = np.power(10.0, np.linspace(np.log10(np.min(St)), np.log10(np.max(St)), 100, endpoint=True))
ABL = ABLf(R, StBL)


# Plot

fig = plt.figure()

plt.plot(StBL, ABL, label='Belyaev & Levin')
plt.plot(St, A, '-ok', label='AeroSolved')


# Style/save

plt.xscale('log')

plt.xlabel(r'$\mathrm{St}$')
plt.ylabel(r'$A$')

fs.post(fig, plt.legend())

plt.savefig('plot.pdf')
