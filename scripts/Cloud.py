#!/usr/bin/env python3

# Script to compute the dispersed inlet velocity in the Cloud case

# Imports --------------------------------------------------------

import numpy as np
import os
import thermophysicalFunctions as tpf

# Parameters -----------------------------------------------------

Z = 0.1     # Dispersed inlet mass fraction
phi = 0.4   # Dispersed volumetric flow rate [ml/min]
R = 0.01    # Inlet radius [m]
T = 293.15  # Temperature [K]

# Script ---------------------------------------------------------

A = np.pi*R**2

phiSI = phi/60*1E-6

rhod = tpf.rhoWaterl(T)
rhoc = tpf.rhoAirg(T)

Y = 1.0-Z

alphad = Z/rhod/(Z/rhod+Y/rhoc)

Ud = phiSI/A/alphad

# Print ----------------------------------------------------------

print('dispersed inlet mass fraction Z =', Z)
print('dispersed inlet void fraction alpha =', alphad)
print('dispersed inlet velocity Ud =', Ud, 'm/s')
print('mixture inlet velocity U =', Ud*Z, 'm/s')
