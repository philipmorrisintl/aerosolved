#!/usr/bin/python

import sys
import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, '../../scripts')
import figStyle as fs

# Prepare figure
fs.prep(plt)

# Read data
data = np.loadtxt('postProcessing/probes/0/dcm')

# Time
t_aero = data[:,0]

# d² in µm²
d2_um2 = np.square(data[:,1]) * 1e12

# -----------------------------
# SAVE PLOT DATA
# -----------------------------
np.savetxt(
    "AeroSolved_plotData_um2.txt",
    np.column_stack((t_aero, d2_um2)),
    header="t [s]    d^2 [micron^2]",
    comments=''
)

# -----------------------------
# Plot
# -----------------------------
fig = plt.figure()

plt.plot(t_aero, d2_um2, label='AeroSolved')

# Labels
plt.xlabel(r'$t$ [s]')
plt.ylabel(r'$d^2$ [$\mu$m$^2$]')

# Style/save
fs.post(fig, plt.legend())

# PNG + PDF
plt.savefig('plot.png', dpi=600, bbox_inches='tight')
plt.savefig('plot.pdf', bbox_inches='tight')
