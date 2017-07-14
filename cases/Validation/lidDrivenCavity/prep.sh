#!/bin/bash

source $FOAM_SRC/../bin/tools/RunFunctions
source $FOAM_SRC/../bin/tools/CleanFunctions

cp -r 0.org 0

runApplication blockMesh

runApplication decomposePar
