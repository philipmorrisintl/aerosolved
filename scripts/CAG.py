#!/usr/bin/python

# Script to compute variables for the CAG case

# Imports --------------------------------------------------------

import numpy as np
import os
import thermophysicalFunctions as tpf

# Parameters -----------------------------------------------------

theta = 5               # Wedge angle [degrees]

alphaPG = 0.7           # Liquid mixture PG volume fraction
alphaVG = 0.2           # Liquid mixture VG volume fraction
alphaWater = 0.1        # Liquid mixture water volume fraction

flowLiquid = 0.5E-3     # Liquid mixture flow rate [l/min]

TCold = 273.15 + 20     # Cold temperature [K]
THot = 273.15 + 250     # Hot temperature [K]

flowAirHot = 2.0        # Hot air flow rate [l/min]
flowAirCold = 10.0      # Cold air flow rate [l/min

r1 = 0.33E-3            # Vapor inlet pipe radius [m]

# Script ---------------------------------------------------------

# Wedge fraction

f = theta/360

# Compute mass flow rates [kg/s]

massFlowPG = flowLiquid*1E-3/60*alphaPG*tpf.rhoPGl(TCold)*f
massFlowVG = flowLiquid*1E-3/60*alphaVG*tpf.rhoVGl(TCold)*f
massFlowWater = flowLiquid*1E-3/60*alphaWater*tpf.rhoWaterl(TCold)*f

massFlowAirHot = flowAirHot*1E-3/60*tpf.rhoAirg(TCold)*f
massFlowAirCold = flowAirCold*1E-3/60*tpf.rhoAirg(TCold)*f

massFlowVapor = massFlowPG + massFlowVG + massFlowWater

# Compute vapor mass fractions

YPG = massFlowPG/massFlowVapor
YVG = massFlowVG/massFlowVapor
YWater = massFlowWater/massFlowVapor

# Compute vapor volume flow rates [m^3/s]

volumeFlowPG = massFlowPG/tpf.rhoPGg(THot)
volumeFlowVG = massFlowVG/tpf.rhoVGg(THot)
volumeFlowWater = massFlowWater/tpf.rhoWaterg(THot)

volumeFlowVapor = volumeFlowPG + volumeFlowVG + volumeFlowWater

# Compute vapor jet Reynolds number

nu = YPG*tpf.muPGg(THot)/tpf.rhoPGg(THot) + YVG*tpf.muVGg(THot)/tpf.rhoVGg(THot) + YWater*tpf.muWaterg(THot)/tpf.rhoWaterg(THot)

D = 2.0*r1
A = np.pi/4.0*D**2
U = volumeFlowVapor/f/A

ReVaporJet = D*U/nu

# Compute estimate of outlet temperature and densities

massFlow = massFlowAirHot + massFlowAirHot + massFlowVapor

TMix = (massFlowAirCold*TCold + massFlowAirHot*THot + massFlowVapor*THot)/massFlow

rhoAirMix = tpf.rhoAirg(TMix)

rhoPGgMix = tpf.rhoPGg(TMix)
rhoVGgMix = tpf.rhoVGg(TMix)
rhoWatergMix = tpf.rhoWaterg(TMix)

rhoPGlMix = tpf.rhoPGl(TMix)
rhoVGlMix = tpf.rhoVGl(TMix)
rhoWaterlMix = tpf.rhoWaterl(TMix)

# Print

print('cold temperature =', TCold)
print('hot temperature =', THot);

print('mass flow hot air =', massFlowAirHot)
print('mass flow cold air =', massFlowAirCold)
print('mass flow vapor =', massFlowVapor)

print('mass fraction PG vapor =', YPG)
print('mass fraction VG vapor =', YVG)
print('mass fraction Water vapor =', YWater)

print('vapor jet Reynolds number =', ReVaporJet)

print('mixture temperature =', TMix)

print('density Air at Tmix =', rhoAirMix)

print('density PG gas at TMix =', rhoPGgMix)
print('density VG gas at TMix =', rhoVGgMix)
print('density Water gas at TMix =', rhoWatergMix)

print('density PG liquid at TMix =', rhoPGlMix)
print('density VG liquid at TMix =', rhoVGlMix)
print('density Water liquid at TMix =', rhoWaterlMix)

print('flowLiquid/massFlowVapor =', flowLiquid/massFlowVapor)
print('flowAirCold/massFlowAirCold =', flowAirCold/massFlowAirCold)
