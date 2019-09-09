#!/usr/bin/python

# Script to compute variables for the uniformNucleation case

# Imports --------------------------------------------------------

import numpy as np
import os
import thermophysicalFunctions as tpf

# Parameters -----------------------------------------------------

alphaPG = 0.7           # Liquid mixture PG volume fraction
alphaVG = 0.2           # Liquid mixture VG volume fraction
alphaWater = 0.1        # Liquid mixture water volume fraction

flowLiquid = 0.5E-3     # Liquid mixture flow rate [l/min]

TCold = 273.15 + 20     # Cold temperature [K]
THot = 273.15 + 250     # Hot temperature [K]

flowAirHot = 2.0        # Hot air flow rate [l/min]
flowAirCold = 10.0      # Cold air flow rate [l/min]

# Script ---------------------------------------------------------

# Mixture temperature based on air only

TMix = (flowAirCold*TCold + flowAirHot*THot)/(flowAirHot+flowAirCold)

# Compute mass flow rates [g/min]

massFlowPG = flowLiquid*alphaPG*tpf.rhoPGl(TCold)
massFlowVG = flowLiquid*alphaVG*tpf.rhoVGl(TCold)
massFlowWater = flowLiquid*alphaWater*tpf.rhoWaterl(TCold)

massFlowAirHot = flowAirHot*tpf.rhoAirg(TCold)
massFlowAirCold = flowAirCold*tpf.rhoAirg(TCold)

totalMassFlow = massFlowPG+massFlowVG+massFlowWater+massFlowAirHot+massFlowAirCold

# Compute mass fractions of supersaturated vapor mixture

YPG = massFlowPG/totalMassFlow
YVG = massFlowVG/totalMassFlow
YWater = massFlowWater/totalMassFlow
YAir = (massFlowAirHot+massFlowAirCold)/totalMassFlow

# Print

print('THot =', THot)
print('TMix =', TMix)
print('PG mass fraction =', YPG)
print('VG mass fraction =', YVG)
print('Water mass fraction =', YWater)
print('Air mass fraction =', YAir)
