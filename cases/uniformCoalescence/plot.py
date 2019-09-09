#!/usr/bin/python

import sys
import numpy as np
import matplotlib.pyplot as plt
import scipy.special as sps

sys.path.insert(0, '../../scripts')

import figStyle as fs

sigmaInf = 1.320
sigmaInfFM = 1.355
A = 1.591


# Prepare figure

fs.prep(plt)

colors = ['#1f77b4', '#ff7f0e']


# Load parameters

data = np.loadtxt('params.txt')

Kn0 = data[0]
K0 = data[1]
Kfm = data[2]
M0 = data[3]
tau = data[4]
T = data[5]
l = data[6]
sigma0 = data[7]
yMin = data[8]
yMax = data[9]
Nsec = int(data[10])
rhol = data[11]

# Load simulation data and plot

fig = plt.figure('N')

data = np.loadtxt('postProcessing/probes/0/M')

t = data[:,0]/tau
f = data[:,1]/M0

plt.plot(t, f, label='AeroSolved')


fig = plt.figure('dist')

plotTimes = [0, -1]

if sys.argv[1] == "sectional":

    logy = np.linspace(np.log10(yMin), np.log10(yMax), Nsec+1, endpoint=True)

    y = np.power(10.0, logy)
    x = np.power(10.0, (logy[1:Nsec+1]+logy[0:Nsec])/2.0)

    dx = np.power(x/rhol*6.0/np.pi, 1.0/3.0)
    dy = np.power(y/rhol*6.0/np.pi, 1.0/3.0)

    rho = np.loadtxt('postProcessing/probes/0/rho')[:,1]
    time = np.loadtxt('postProcessing/probes/0/rho')[:,0]

    ii = 0

    for i in plotTimes:

        N = np.zeros(Nsec)

        for j in range(1,Nsec):

            fileName = 'M.' + str(j).zfill(int(np.log10(Nsec))+1);
            M = np.loadtxt('postProcessing/probes/0/' + fileName)[i,1]
            N[j] = M*rho[i]

        # Plot

        dlogd = np.log10(dy[1:Nsec+1])-np.log10(dy[0:Nsec])

        plt.plot(dx, N/dlogd/sum(N), color=colors[ii], label=r'$t/\tau=' + str(round(time[i]/tau)) + '$')

        ii = ii+1


else:

    rho = np.loadtxt('postProcessing/probes/0/rho')[:,1]
    N = np.loadtxt('postProcessing/probes/0/M')[:,1]*rho
    Z = np.loadtxt('postProcessing/probes/0/Water.dispersed')[:,1]
    CMD = np.loadtxt('postProcessing/probes/0/CMD')[:,1]
    time = np.loadtxt('postProcessing/probes/0/rho')[:,0]

    ii = 0

    for i in plotTimes:

        P = 1024

        dy = np.power(10.0, np.linspace(-10, -4, P, endpoint=True))

        dx = np.power(10.0, (np.log10(dy[1:P])+np.log10(dy[0:P-1]))*0.5)

        dlogd = np.log10(dy[1:P])-np.log10(dy[0:P-1])

        df = 0.5 * ( sps.erf((np.log(dy[1:P])-np.log(CMD[i]))/(np.sqrt(2)*np.log(sigma0))) - sps.erf((np.log(dy[0:P-1])-np.log(CMD[i]))/(np.sqrt(2)*np.log(sigma0))) )

        # Plot

        plt.plot(dx, df/dlogd, color=colors[ii], label=r'$t/\tau=' + str(round(time[i]/tau)) + '$')

        ii = ii+1


# Compute analytical solution of Park et al. (1999)

def Park(s, p, q, f):

    return \
        12.0*s/5.0 * (f**(-5.0/6.0)-1.0) \
      + 3.0/p \
      * ( \
            1.0/3.0*(f**(-1)-1.0) \
          - (q/p)*0.5*(f**(-2.0/3.0)-1.0) \
          + (q/p)**2*(f**(-1.0/3.0)-1.0) \
          - (q/p)**3*np.log((p+q*f**(1.0/3.0))/((p+q)*f**(1.0/3.0))) \
        )

Z0 = np.log(sigma0)**2
Zinf = np.log(sigmaInf)**2
ZinfFM = np.log(sigmaInfFM)**2

p = 1.0 + np.exp(Z0)
q = A*Kn0*np.exp(1.5*(Zinf-Z0))*(np.exp(0.5*Z0)+np.exp(2.5*Z0))
H = np.exp(25.0*Z0/8.0)+2.0*np.exp(5.0*Z0/8.0)+np.exp(Z0/8.0)
b0 = 1.0+1.2*np.exp(-2.0*sigma0)-0.646*np.exp(-0.35*sigma0**2)

dg0 = 2.0*l/Kn0
vg0 = np.pi/6.0*dg0**3

s = K0/(2.0*b0*Kfm*H*vg0**(1.0/6.0))

f = np.power(10, np.linspace(0, -3, 1024))

t = Park(s, p, q, f)

fig = plt.figure('N')

plt.plot(t, f, '-ok', markevery=0.05, lw=0.5, label='Park et al. (1999)')


fig = plt.figure('dist')

ii = 0

for i in plotTimes:

    ti = time[i]/tau

    idx = np.abs(t-ti).argmin()
    fi = f[idx]

    d = \
        ( \
            2.0*s*fi**(1.0/6.0) \
          + 1.0/(p+q*fi**(1.0/3.0)) \
        ) \
      / ( \
            2.0*s*fi**(1.0/6.0)*np.exp(-1.5*ZinfFM) \
          + 1.0/(p+q*fi**(1.0/3.0)*np.exp(-3.0*Zinf)) \
        )

    Z = 1.0/9.0*np.log(2.0*d+fi*(np.exp(9.0*Z0)-2.0*d))

    vg = np.exp(4.5*Z0)/fi/np.sqrt(2*d+fi*(np.exp(9.0*Z0)-2.0*d))*vg0

    dg = (vg*6.0/np.pi)**(1.0/3.0)

    P = 1024

    dy = np.power(10.0, np.linspace(-10, -4, P, endpoint=True))
    dx = np.power(10.0, (np.log10(dy[1:P])+np.log10(dy[0:P-1]))*0.5)
    dlogd = np.log10(dy[1:P])-np.log10(dy[0:P-1])

    df = 0.5 * ( sps.erf((np.log(dy[1:P])-np.log(dg))/(np.sqrt(2*Z))) - sps.erf((np.log(dy[0:P-1])-np.log(dg))/(np.sqrt(2*Z))) )

    # Plot

    plt.plot(dx, df/dlogd, '-o', lw=0.5, color=colors[ii], markevery=0.05)

    ii = ii+1

plt.plot([], [], '-k', label='AeroSolved')
plt.plot([], [], '-ok', lw=0.5, label='Park et al. (1999)')


# Style/save

fig = plt.figure('N')

plt.xlabel(r'$t/\tau$')
plt.ylabel(r'$N/N_0$')

plt.xscale('log')
plt.yscale('log')

fs.post(fig, plt.legend())

plt.savefig('N.pdf')


fig = plt.figure('dist')

plt.xlabel(r'$d$ [m]')
plt.ylabel(r'$\mathrm{d}f/\mathrm{d}\log d$')

plt.xscale('log')

fs.post(fig, plt.legend())

plt.savefig('dist.pdf')
