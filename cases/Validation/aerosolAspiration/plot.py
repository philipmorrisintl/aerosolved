#!/usr/bin/python
import sys
import numpy as np
import matplotlib as mpl

mpl.use('Agg')

import matplotlib.pyplot as plt

plt.rc('text', usetex=True)

rhol = 3960.0
D = 2E-3
mu = 1.81E-5
rho = 1.0
U = 1.0
Re = D*U*rho/mu

# Compute Belyaev & Levin curve

def GBL(R): 
    return (2.0 + 0.62/R)

def aBL(R, St): 
    return (1 - 1/(1+GBL(R)*St))

def ABLf(R, St):
    return (1 + aBL(R,St)*(R-1))

if len(sys.argv) > 1:

    data = np.loadtxt('simdata.txt')

    x = data[:,1]
    A = data[:,3]

    d = np.power(x/rhol*6.0/np.pi, 1.0/3.0)

    U = Re*mu/rho/D

    tau = np.power(d, 2.0)*rhol/(18.0*mu)
    St = tau*U/D

    R = data[0,2]

    print 'Aspiration rate R =', R

    StBL = np.power(10.0, np.linspace(np.log10(np.min(St)), np.log10(np.max(St)), 10, endpoint=True))
    ABL = ABLf(R, StBL)

    plt.figure()

    plt.plot(St, A, label='Simulation')

    plt.plot(StBL, ABL, 'ok', label='Belyaev \& Levin')

    plt.xscale('log')

    plt.xlabel(r'$St$', fontsize=20)
    plt.ylabel(r'$A$', fontsize=20)

    plt.legend(loc='best')

    plt.savefig('plot.pdf')

else:

    print('Please specify the time to be used')


