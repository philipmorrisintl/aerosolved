#!/usr/bin/env python3

import sys
import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, '../../scripts')

import figStyle as fs

yMin = 1E-24
yMax = 1E-10
N = 10
rhol = 1000

# Prepare figure

fs.prep(plt)


# Load data

topNumberFlux = np.loadtxt('postProcessing/numberFlux/0/patch.top.dat')
bottomNumberFlux = np.loadtxt('postProcessing/numberFlux/0/patch.bottom.dat')

topMassFlux = np.loadtxt('postProcessing/massFlux/0/patch.top.dat')
bottomMassFlux = np.loadtxt('postProcessing/massFlux/0/patch.bottom.dat')

times = topNumberFlux[:,0]

# Plot

logy = np.linspace(np.log(yMin), np.log(yMax), N+1, endpoint=True)

x = np.exp((logy[:-1]+logy[1:])/2.0)
d = (x/rhol*6.0/np.pi)**(1/3)

fig = plt.figure('number')

i1 = int(len(times)/3)-1
i2 = int(len(times)/3*2)-1
i3 = int(len(times))-1

plt.plot(d, np.abs(topNumberFlux[i1,1:]), '-o', color='C0', label='influx, $t =' + str(times[i1]) + '$ s', mew=0.25, mec='white', ms=4)
plt.plot(d, np.abs(bottomNumberFlux[i1,1:]), '-o', color='C1', label='outflux, $t =' + str(times[i1]) + '$ s', mew=0.25, mec='white', ms=4)

plt.plot(d, np.abs(topNumberFlux[i2,1:]), '-d', color='C0', label='influx, $t =' + str(times[i2]) + '$ s', mew=0.25, mec='white', ms=4)
plt.plot(d, np.abs(bottomNumberFlux[i2,1:]), '-d', color='C1', label='outflux, $t =' + str(times[i2]) + '$ s', mew=0.25, mec='white', ms=4)

plt.plot(d, np.abs(topNumberFlux[i3,1:]), '-s', color='C0', label='influx, $t =' + str(times[i3]) + '$ s', mew=0.25, mec='white', ms=4)
plt.plot(d, np.abs(bottomNumberFlux[i3,1:]), '-s', color='C1', label='outflux, $t =' + str(times[i3]) + '$ s', mew=0.25, mec='white', ms=4)

##

fig = plt.figure('mass')

topMassFluxAr = np.array([topMassFlux[-1,1], topMassFlux[-1,2], topMassFlux[-1,3], np.sum(topMassFlux[-1,1:])])
bottomMassFluxAr = np.array([bottomMassFlux[-1,1], bottomMassFlux[-1,2], bottomMassFlux[-1,3], np.sum(bottomMassFlux[-1,1:])])

ind = np.arange(4)

plt.bar(ind, topMassFluxAr*1E9, 0.4, yerr=0.001, label='top')
plt.bar(ind+4, bottomMassFluxAr*1E9, 0.4, yerr=0.001, label='bottom')

# Style/save

fig = plt.figure('number')

plt.xscale('log')
plt.yscale('log')

plt.xlabel(r'$d$ [m]')
plt.ylabel(r'particle flux [#/s]')

fs.post(fig, plt.legend())

plt.savefig('number.pdf')

##

fig = plt.figure('mass')

plt.ylabel(r'mass flux [$\mathrm{\mu}$g/s]')
plt.xticks(np.arange(8), ['vapor', 'air', 'droplet', 'sum', 'vapor', 'air', 'droplet', 'sum'], fontsize=6)

plt.axhline(0, color='black', lw=0.5)

fs.post(fig, plt.legend())

plt.savefig('mass.pdf')
