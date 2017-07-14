#!/bin/bash

. $WM_PROJECT_DIR/bin/tools/RunFunctions

NPROC=$(getNumberOfProcessors)

runParallel aerosolEulerFoam $NPROC
