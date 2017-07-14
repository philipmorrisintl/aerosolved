#!/bin/bash





source $FOAM_SRC/../bin/tools/RunFunctions
source $FOAM_SRC/../bin/tools/CleanFunctions

r=$(echo "print 0.5*4E-3*5" | python)

NCELLS=$(cat constant/polyMesh/blockMeshDict | grep "Number of cells:" | cut -c 22-)

echo D = 4E-3, RS = 5, r = $r

runApplication blockMesh

cp -r 0.org 0

runApplication bentPipeDeform $r

runApplication decomposePar
