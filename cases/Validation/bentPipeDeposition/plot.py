#!/usr/bin/python
import sys
import numpy as np
import matplotlib as mpl

mpl.use('Agg')

import matplotlib.pyplot as plt

plt.rc('text', usetex=True)

rho = 1.189799936415828
rhol = 1E3
mu = 1.81e-5
D = 4E-3
Re = 1000
U = Re*mu/D

x = np.array([4.213279e-15, 2.189284e-14, 4.710589e-14, 7.803102e-14, 1.137585e-13, 1.537125e-13, 1.974855e-13, 2.447694e-13, 2.953205e-13, 3.489399e-13, 4.054611e-13, 4.64742e-13, 5.266599e-13, 5.911067e-13, 6.579868e-13, 7.272149e-13, 7.987137e-13, 8.724133e-13, 9.482499e-13, 1.026165e-12, 1.106104e-12, 1.188018e-12, 1.271859e-12, 1.357585e-12, 1.445155e-12, 1.53453e-12, 1.625676e-12, 1.718558e-12, 1.813145e-12, 1.909406e-12, 2.007312e-12, 2.106837e-12])

d = np.power(x/rhol*6.0/np.pi, 1.0/3.0)
tau = np.square(d)*rhol/(18.0*mu)
St = tau*U/(D/2.0)

print St

# Load Pui's reference data

data = np.loadtxt('referenceData/Pui.txt')

PuiSt = data[:,0]
PuiDE = data[:,1]

if len(sys.argv) > 1:

    data = np.loadtxt('simdata.txt')
    DE = -data[:,2]/data[:,0]*100

    plt.figure()

    plt.plot(St, DE, label='Simulation')

    plt.plot(PuiSt, PuiDE, 'ko', label='Pui et al.')

    plt.xlabel(r'$St$', fontsize=20)
    plt.ylabel(r'$DE$ [\%]', fontsize=20)

    plt.legend(loc='best')

    plt.savefig('plot.pdf')

else:

    print('Please specify the time to be used')


