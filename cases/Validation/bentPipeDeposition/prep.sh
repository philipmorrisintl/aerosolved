#!/bin/bash

source $FOAM_SRC/../bin/tools/RunFunctions
source $FOAM_SRC/../bin/tools/CleanFunctions

r=$(echo "print 0.5*4E-3*5.7" | python)

runApplication blockMesh

cp -r 0.org 0

runApplication bentPipeDeform $r

runApplication decomposePar
