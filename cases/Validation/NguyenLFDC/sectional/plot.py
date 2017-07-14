#!/usr/bin/python
import sys
import numpy as np
import matplotlib as mpl

mpl.use('Agg')

import matplotlib.pyplot as plt

plt.rc('text', usetex=True)

# Compute droplet sizes

yMin = 1E-24
yMax = 1E-12
P = 100
rhol = 1.0453e+03
A = 1.089446783e-06

logy = np.linspace(np.log10(yMin), np.log10(yMax), P+1, endpoint=True)

y = np.power(10.0, logy)
x = np.power(10.0, (logy[1:P+1]+logy[0:P])/2.0)

dx = np.power(x/rhol*6.0/np.pi, 1.0/3.0)
dy = np.power(y/rhol*6.0/np.pi, 1.0/3.0)

# Load Nguyen's reference data

data = np.loadtxt('referenceData/Nguyen.txt')

NguyenD = data[:,0]
NguyenN = data[:,1]

if len(sys.argv) > 1:

    rho = np.loadtxt('postProcessing/probes/' + str(sys.argv[1]) + '/rho')[1]

    N = np.zeros(P)

    for i in range(0,P):

        M = np.loadtxt('postProcessing/probes/' + str(sys.argv[1]) + '/M.' + str(i).zfill(2))[1]
        N[i] = M*rho

    print N

    plt.figure()

    dlogd = np.log10(dy[1:P+1])-np.log10(dy[0:P])

    plt.plot(dx, N/dlogd/1E6, label='Simulation (sectional)')

    plt.plot(NguyenD, NguyenN, 'ko', label='Nguyen et al.')

    plt.xlabel(r'$d$ [m]', fontsize=20)
    plt.ylabel(r'$\mathrm{d}N/\mathrm{d}\log d$ [1/cm${}^3$]', fontsize=20)

    plt.legend(loc='best')

    plt.xscale('log')
    plt.yscale('log')

    plt.gca().set_ylim([1E3, 3E7])
    plt.gca().set_xlim([1E-7, 1E-5])

    plt.savefig('plot-sectional.pdf')

else:

    print('Please specify the time to be used')


