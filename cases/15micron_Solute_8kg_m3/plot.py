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

# Prepare plot data
t_aero = data[:,0]

# d^2 in micron^2
d2_um2 = np.square(data[:,1]) * 1e12     # µm^2

# d^2 in micron^2 × 10^4
d2_um2_1e4 = d2_um2 / 1e4                 # (µm^2 × 10^4)

# -----------------------------
# SAVE PLOT DATA
# -----------------------------
np.savetxt(
    "AeroSolved_plotData_um2.txt",
    np.column_stack((t_aero, d2_um2)),
    header="t [s]    d^2 [micron^2]",
    comments=''
)

np.savetxt(
    "AeroSolved_plotData_um2_x1e4.txt",
    np.column_stack((t_aero, d2_um2_1e4)),
    header="t [s]    d^2 [micron^2 x 1e4]",
    comments=''
)

# -----------------------------
# Plot (choose which one)
# -----------------------------
fig = plt.figure()

# Option A: direct micron^2
# plt.plot(t_aero, d2_um2, label='AeroSolved')

# Option B: micron^2 × 10^4 
plt.plot(t_aero, d2_um2_1e4, label='AeroSolved')

# Style/save
plt.xlabel(r'$t$ [s]')
plt.ylabel(r'$d^2$ [$\mu$m$^2 \times 10^{4}$]')

fs.post(fig, plt.legend())

# PNG + PDF
plt.savefig('plot.png', dpi=600, bbox_inches='tight')
plt.savefig('plot.pdf', bbox_inches='tight')

