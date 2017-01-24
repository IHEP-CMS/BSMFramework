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
cmsRun $CMSSW_SRC/$PROJECTDIR/python/miniAOD_RunH_DL.py
cmsRun $CMSSW_SRC/$PROJECTDIR/python/miniAOD_RunH_SL.py

#===== Output file =======
# N.B: must be the same as output filename set in config file code
# ========================
#mv Data_RunsB-G_SL.root $CMSSW_SRC/$PROJECTDIR/Data_RunBCD_SL.root
#mv Data_RunsB-G_SL.root $CMSSW_SRC/$PROJECTDIR/Data_RunBCD_DL.root
#mv Data_RunH_SL.root $CMSSW_SRC/$PROJECTDIR/Data_RunEF_SL.root
#mv Data_RunH_DL.root $CMSSW_SRC/$PROJECTDIR/Data_RunEF_DL.root
#mv Data_RunH_SL.root $CMSSW_SRC/$PROJECTDIR/Data_RunG_SL.root
#mv Data_RunH_DL.root $CMSSW_SRC/$PROJECTDIR/Data_RunG_DL.root
#mv Data_RunH_SL.root $CMSSW_SRC/$PROJECTDIR/Data_RunH_SL.root
#mv Data_RunH_DL.root $CMSSW_SRC/$PROJECTDIR/Data_RunH_DL.root
