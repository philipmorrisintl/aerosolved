#!/usr/bin/env python3

import sys
import csv
import numpy as np
import matplotlib.pyplot as plt
import subprocess
import glob

def prep(p):

    p.rc('font', family='serif', size=9, serif='STIXGeneral')
    p.rc('mathtext', fontset='stix')
    p.rc('axes', labelsize=9, titlesize=9)
    p.rc('lines', dash_capstyle='round', linewidth=1, color='black', markersize=3)

    p.rc('xtick.major', size=4, width=0.25, pad=2)
    p.rc('xtick.minor', size=2, width=0.25, pad=2)
    p.rc('ytick.major', size=4, width=0.25, pad=2)
    p.rc('ytick.minor', size=2, width=0.25, pad=2)

    p.rc('legend', fontsize=6, handlelength=2, numpoints=1, labelspacing=0.2, borderpad=0.2, fancybox=False, )

    p.rc('savefig', format='pdf')

    p.rc('figure', max_open_warning=0)

    p.autoscale(tight=True)

    p.rc('figure', figsize=(2.75, 2.55))


def post(f, l=False):

    ax = f.gca()

    from matplotlib.ticker import ScalarFormatter

    for pos in ['bottom', 'top', 'right', 'left']:
        ax.spines[pos].set_linewidth(0.5)

    if l:
        l.get_frame().set_linewidth(0.25)
        l.get_frame().set_edgecolor('black')

    figSize = f.get_size_inches()

    ax.set_position([0.25, 0.16, 0.72, 0.72*figSize[0]/figSize[1]])

    if ax.get_yscale() != 'log':
        fmt = ScalarFormatter()
        fmt.set_powerlimits((-2,2))
        ax.yaxis.set_major_formatter(fmt)

    return


def readFile(filename):
    data = []
    header = []
    with open(filename) as csvDataFile:
        csvReader = csv.reader(csvDataFile)
        header = next(csvReader,None)
        for row in csvReader:
#            row=reformat_row(row)
            for i in range(len(row)): row[i]=float(row[i])
            data.append(row)
    return header,data


# Prepare figure

prep(plt)

rhol=997

x = np.array([ 5.2203E-19,4.1762E-18,1.4095E-17,3.3410E-17,6.5253E-17,1.1276E-16,1.7906E-16,2.6728E-16,3.8056E-16,5.2203E-16,4.1762E-15,1.4095E-14,3.3410E-14,6.5253E-14,1.1276E-13,1.7906E-13,2.6728E-13,3.8056E-13,5.2203E-13,4.1762E-12,1.4095E-11,3.3410E-11,6.5253E-11,1.1276E-10,1.7906E-10,2.6728E-10,3.8056E-10,5.2400E-10,4.1900E-09,1.4140E-08,3.3500E-08,6.5500E-08,1.1310E-07,1.7960E-07,2.6800E-07,3.8200E-07,5.2400E-07,9.0500E-07,1.4370E-06,2.1400E-06,3.0500E-06,4.1900E-06,5.5800E-06,7.2400E-06,9.2000E-06,1.1490E-05,1.4140E-05   ]);

d = np.power(x/rhol*6.0/np.pi, 1.0/3.0)*1E6
print(d)


# Read data
header,refData = readFile('referenceData/TerminalVelocity-Gunn1949-table2.csv')
refData=np.array(refData)

vs=[]

i=0

Vlist=glob.glob("postProcessing/probes/0/V.*")
subprocess.call(["sed -i -e 's/(/ /g'  postProcessing/probes/0/V.* "], shell=True)
subprocess.call(["sed -i -e 's/)/ /g'  postProcessing/probes/0/V.* "], shell=True)

Vlist.sort()

for v in Vlist:
    filedata = np.loadtxt(v)
    vs.append(float(-filedata[-1,2]))
    i+=1

vs=np.array(vs)

# Plot

fig = plt.figure()


plt.plot(refData[:,0], refData[:,1]/100, 'ok', label='Gunn & Kinzer (1949)')
plt.plot(d[:], vs[:],'.r', label='AeroSolved')

# Style/save

plt.xlim(0.1,3000 ) 
plt.ylim(0.0,10.0 )

plt.xlabel(r'$d$ [$\mu$m]')
plt.ylabel(r'$v$ [m/s]')

post(fig, plt.legend())

plt.savefig('plot.pdf')



fig2 = plt.figure()

v_settling=np.power(d,2)/18*996*9.81/(1.8e-5)*1E-12 

plt.plot(d[:], v_settling[:],label=r'$g \rho d^2 /  \mu $')
plt.plot(refData[:,0], refData[:,1]/100, 'ok', label='Gunn & Kinzer (1949)')
plt.plot(d[:], vs[:], '.r', label='AeroSolved')
# Style/save

plt.xlim(0.1,3000) 

plt.xlabel(r'$d$ [$\mu$m]')
plt.ylabel(r'$v$ [m/s]')

post(fig2, plt.legend())
plt.xscale('log')
plt.yscale('log')
plt.savefig('plotLog.pdf')
