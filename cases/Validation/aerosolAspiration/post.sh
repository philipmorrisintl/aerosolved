#!/bin/bash

if [ -d "processor0" ]; then

    echo Reconstructing latest time...

    reconstructPar -latestTime > log.reconstructPar 2>&1
fi

TIME=$(foamInfoExec -latestTime)

if [ ! -z "$TIME" ]; then

    echo Sampling at time = $TIME...

    sampleAspirationWedge -time $TIME | grep Section | cut -c 10- > simdata.txt

    echo Writing plots...

    python plot.py $TIME > log.python 2>&1

fi

echo Done
