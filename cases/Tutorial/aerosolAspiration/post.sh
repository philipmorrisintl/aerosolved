#!/bin/bash

reconstructPar -latestTime

sample -latestTime

LATESTTIME=$(foamInfoExec -latestTime)

#python plot.py $LATESTTIME

