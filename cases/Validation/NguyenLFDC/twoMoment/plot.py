#!/usr/bin/python
import sys
import numpy as np
import scipy.special as sps
import matplotlib as mpl

mpl.use('Agg')

import matplotlib.pyplot as plt

plt.rc('text', usetex=True)

rhol = 1.0453e+03
sigmag = 1.33333
sigma = np.log(sigmag)

# Load Nguyen's reference data

data = np.loadtxt('referenceData/Nguyen.txt')

NguyenD = data[:,0]
NguyenN = data[:,1]

if len(sys.argv) > 1:

    rho = np.loadtxt('postProcessing/probes/' + str(sys.argv[1]) + '/rho')[1]
    M = np.loadtxt('postProcessing/probes/' + str(sys.argv[1]) + '/M')[1]
    DBPZ = np.loadtxt('postProcessing/probes/' + str(sys.argv[1]) + '/DBPZ')[1]

    N = M*rho

    W = rho * DBPZ / rhol

    CMD = np.power(6.0*W/(np.pi*N), 1.0/3.0) * np.exp(-3.0/2.0 * np.square(np.log(sigmag)))

    P = 100

    dy = np.power(10.0, np.linspace(-7, -5, P, endpoint=True))

    dx = np.power(10.0, (np.log10(dy[1:P])+np.log10(dy[0:P-1]))*0.5)

    dlogd = np.log10(dy[1:P])-np.log10(dy[0:P-1])

    dN = 0.5 * N * ( sps.erf((np.log(dy[1:P])-np.log(CMD))/(np.sqrt(2)*sigma)) - sps.erf((np.log(dy[0:P-1])-np.log(CMD))/(np.sqrt(2)*sigma)) )

    plt.figure()

    print dN/dlogd

    plt.plot(dx, dN/dlogd/1E6, label='Simulation (moment model)')
    plt.plot(NguyenD, NguyenN, 'ko', label='Nguyen et al.')

    plt.xlabel(r'$d$ [m]', fontsize=20)
    plt.ylabel(r'$\mathrm{d}N/\mathrm{d}\log d$ [1/cm${}^3$]', fontsize=20)

    plt.legend(loc='best')

    plt.xscale('log')
    plt.yscale('log')

    plt.gca().set_ylim([1E3, 3E7])

    plt.savefig('plot-twoMoment.pdf')

else:

    print('Please specify the time to be used')


