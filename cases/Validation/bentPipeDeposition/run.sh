#!/bin/bash

source $FOAM_SRC/../bin/tools/RunFunctions
source $FOAM_SRC/../bin/tools/CleanFunctions

APP=$(getApplication)
NPROC=$(getNumberOfProcessors)

runParallel $APP $NPROC
