#!/bin/bash

. $WM_PROJECT_DIR/bin/tools/RunFunctions

APP=`getApplication`

NPROC=$(getNumberOfProcessors)

runApplication $APP
