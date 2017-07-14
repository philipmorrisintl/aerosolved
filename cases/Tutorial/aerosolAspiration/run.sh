#!/bin/bash

. $WM_PROJECT_DIR/bin/tools/RunFunctions

application=`getApplication`

NPROC=$(getNumberOfProcessors)

runParallel $application $NPROC

#runApplication $application
