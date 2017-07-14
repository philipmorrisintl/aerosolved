#!/bin/bash

source $FOAM_SRC/../bin/tools/RunFunctions
source $FOAM_SRC/../bin/tools/CleanFunctions

runApplication reconstructPar -latestTime

sampleDropletFlux -latestTime | grep -P '^\t' | cut -c 2- > sample.txt
