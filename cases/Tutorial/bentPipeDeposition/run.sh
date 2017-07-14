#!/bin/bash

source $FOAM_SRC/../bin/tools/RunFunctions
source $FOAM_SRC/../bin/tools/CleanFunctions

NPROC=$(getNumberOfProcessors)
APP=$(getApplication)

runParallel $APP $NPROC
