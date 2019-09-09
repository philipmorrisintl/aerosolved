#!/usr/bin/python

import sys
import numpy as np
import matplotlib.pyplot as plt
import scipy.special as sps

sys.path.insert(0, '../../scripts')

import figStyle as fs

rhol = 1.0453e+03
sigma = 1.2
P = 128


# Prepare figure

fs.prep(plt)


# Load simulation data and plot

fig = plt.figure()

if sys.argv[1] == "sectional":

    data = np.loadtxt('params.txt')

    yMin = data[0]
    yMax = data[1]
    P = int(data[2])

    logy = np.linspace(np.log10(yMin), np.log10(yMax), P+1, endpoint=True)

    y = np.power(10.0, logy)
    x = np.power(10.0, (logy[1:P+1]+logy[0:P])/2.0)

    dx = np.power(x/rhol*6.0/np.pi, 1.0/3.0)
    dy = np.power(y/rhol*6.0/np.pi, 1.0/3.0)

    rho = np.loadtxt('postProcessing/probes/0/rho')[-1,1]

    N = np.zeros(P)

    for i in range(1,P):

        fileName = 'M.' + str(i).zfill(int(np.log10(P))+1);
        M = np.loadtxt('postProcessing/probes/0/' + fileName)[-1,1]
        N[i] = M*rho

    # Plot

    dlogd = np.log10(dy[1:P+1])-np.log10(dy[0:P])

    plt.plot(dx, N/dlogd/1E6, label='AeroSolved (sectional)')


else:

    rho = np.loadtxt('postProcessing/probes/0/rho')[-1,1]
    N = np.loadtxt('postProcessing/probes/0/M')[-1,1]*rho
    Z = np.loadtxt('postProcessing/probes/0/DBP.dispersed')[-1,1]
    d = np.loadtxt('postProcessing/probes/0/dcm')[-1,1]

    CMD = d*np.exp(-3.0/2.0*np.square(np.log(sigma)))

    P = 128

    dy = np.power(10.0, np.linspace(-7, -5, P, endpoint=True))

    dx = np.power(10.0, (np.log10(dy[1:P])+np.log10(dy[0:P-1]))*0.5)

    dlogd = np.log10(dy[1:P])-np.log10(dy[0:P-1])

    dN = 0.5 * N * ( sps.erf((np.log(dy[1:P])-np.log(CMD))/(np.sqrt(2)*np.log(sigma))) - sps.erf((np.log(dy[0:P-1])-np.log(CMD))/(np.sqrt(2)*np.log(sigma))) )

    # Plot

    plt.plot(dx, dN/dlogd/1E6, label='AeroSolved (moment model)')


# Load Nguyen's reference data and plot

data = np.loadtxt('referenceData/Nguyen.txt')

NguyenD = data[:,0]
NguyenN = data[:,1]

plt.plot(NguyenD, NguyenN, 'ko', label='Nguyen et al.')


# Style/save

plt.xlabel(r'$d$ [m]')
plt.ylabel(r'$\mathrm{d}N/\mathrm{d}\log d$ [1/cm${}^3$]')

plt.xscale('log')
plt.yscale('log')

plt.gca().set_ylim([1E1, 3E7])
plt.gca().set_xlim([1E-7, 1E-5])

fs.post(fig, plt.legend())

plt.savefig('plot.pdf')
