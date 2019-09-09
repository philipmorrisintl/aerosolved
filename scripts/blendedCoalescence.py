#!/usr/bin/python

# Script which plots the blended coalescence kernel
# This script can be used to optimize the a, b and Kn0 parameters

# Imports --------------------------------------------------------

import sys
import numpy as np
import matplotlib as mpl

import matplotlib.pyplot as plt

import figStyle as fs

# Parameters -----------------------------------------------------

a = 3
offset = 1.5
Kn0 = 5
b = 0.9

W = 28.9596
p = 1e5
T = 293.15
mu = 1.789e-5
rhol = 1e3

kB = 1.38064852e-23
NA = 6.02214086e23

A = 1.591

# Script ---------------------------------------------------------

fs.prep(plt)

fig = plt.figure();

Kns = np.power(10, np.linspace(-4, 4, 1024))

m = W*0.001/NA
lam = mu/p*np.sqrt(np.pi*kB*T/(2*m))

K = 2.0*kB*T/(3.0*mu)

Kt = 3.0*np.sqrt(3.0)/2.0*np.sqrt(mu**2/(rhol*kB*T))*K

wGS = np.array([K, K, 2*A*lam*K, 2*A*lam*K])
pGS = np.array([0, 1, 0, 1])
qGS = np.array([0, -1, -1, -2])

wFM = np.array([b*Kt, b*Kt, b*2.0*Kt])
pFM = np.array([0.5, 2, 1])
qFM = np.array([0, -1.5, -0.5])

fHM = np.zeros(len(Kns))
fKuczaj = np.zeros(len(Kns))

fWeighted = np.zeros(len(Kns))

for i,Kn in enumerate(Kns):

    d = 2.0*lam/Kn

    fGS = sum(wGS*np.power(d, pGS+qGS))
    fFM = sum(wFM*np.power(d, pFM+qFM))

    fHM[i] = 1.0/(1.0/fGS+1.0/fFM)
    fKuczaj[i] = 1.0/(1.0/fGS**2+1.0/fFM**2)**0.5

    phi = fGS/(fFM+fGS)

    fWeighted[i] = phi**2 * fFM + (1-phi)**2 * fGS


plt.plot(Kns, fHM/(2.0*K), label='Harmonic mean')
plt.plot(Kns, fKuczaj/(2.0*K), label='Kuczaj')

plt.plot(Kns, fWeighted/(2.0*K), 'o', markevery=0.05, label='Weighted')


# Style and save figure

plt.xscale('log')
plt.yscale('log')

plt.xlabel(r'Kn')
plt.ylabel(r'$\gamma/2K$')

plt.gca().invert_xaxis()

fs.post(fig, plt.legend())

plt.savefig('blendedCoalescence.pdf')
