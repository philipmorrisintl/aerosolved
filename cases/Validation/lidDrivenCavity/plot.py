#!/usr/bin/python
import sys
import numpy as np
import matplotlib as mpl

mpl.use('Agg')

import matplotlib.pyplot as plt

plt.rc('text', usetex=True)

# Load Ghia's data

GhiaU = np.loadtxt('referenceData/GhiaU.txt')
GhiaV = np.loadtxt('referenceData/GhiaV.txt')

# Load sampled data

if len(sys.argv) > 1:

    SimU = np.loadtxt('postProcessing/sets/' + str(sys.argv[1]) + '/Y_U.xy')
    SimV = np.loadtxt('postProcessing/sets/' + str(sys.argv[1]) + '/X_U.xy')

    plt.figure()

    plt.plot(SimU[:,0], SimU[:,1], label='Simulation')
    plt.plot(GhiaU[:,0], GhiaU[:,-1], 'ok', label='Ghia et al.')

    plt.xlabel(r'$y$', fontsize=20)
    plt.ylabel(r'$u$', fontsize=20)

    plt.legend(loc='best')

    plt.savefig('plotU.pdf')

    plt.figure()

    plt.plot(SimV[:,0], SimV[:,2], label='Simulation')
    plt.plot(GhiaV[:,0], GhiaV[:,-1], 'ok', label='Ghia et al.')

    plt.xlabel(r'$x$', fontsize=20)
    plt.ylabel(r'$v$', fontsize=20)

    plt.legend(loc='best')

    plt.savefig('plotV.pdf')

else:

    print('Please specify the time to be used')


