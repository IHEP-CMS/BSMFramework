#!/bin/bash

# Lxplus Batch Job Script
# > Don't forget . . . whitespace matters ;)

CMSSW_SRC=$(echo ${CMSSW_BASE})/src
PROJECTDIR="BSMFramework/BSM3G_TNT_Maker"

export CMSSW_BASE=$CMSSW_SRC
cd $CMSSW_SRC
eval `scramv1 runtime -sh`

cd $PROJECTDIR
echo $PROJECTDIR

#====== cmsRun Config File =======
cmsRun $CMSSW_SRC/$PROJECTDIR/python/miniAOD_RunEF_DL.py
cmsRun $CMSSW_SRC/$PROJECTDIR/python/miniAOD_RunEF_SL.py

#===== Output file =======
# N.B: must be the same as output filename set in config file code
# ========================
