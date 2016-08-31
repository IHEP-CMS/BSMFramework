#!/bin/bash
#Specify needed variables
varType=int 
varList=(BJetness_isSingleLepton BJetness_isDoubleLepton)
#HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ Muon HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ HLT_IsoMu20 HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL TTHbb Electron HLT_Ele105_CaloIdVT_GsfTrkIdT HLT_Ele27_eta2p1_WP75_Gsf HLT_Ele27_WP85_Gsf HLT_Ele27_eta2p1_WPLoose_Gsf Muon HLT_Mu45_eta2p1 HLT_Mu50 HLT_IsoMu17_eta2p1 HLT_IsoMu24_eta2p1 HLT_IsoMu18 TTHLep Electron HLT_Ele23_WPLoose_Gsf Data HLT_Ele23_CaloIdL_TrackIdL_IsoVL MC Muon HLT_IsoTkMu20 CrossEle-Mu HLT_DiMu9_Ele9_CaloIdL_TrackIdL HLT_Mu8_DiEle12_CaloIdL_TrackIdL HLT_TripleMu_12_10_5 HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL Analysis TTHbb_SL TTHbb_DL TTHLep_2Mu TTHLep_2Ele TTHLep_MuEle TTHLep_3L4L)
#Print info
echo " "
#Declare variables
#echo -e " $varType \c"
pos=0
for count in ${varList[@]}; 
do
  #echo " ${varList[$pos]}, \c"
  echo "  $varType ${varList[$pos]};"
  let pos=pos+1
done
echo " "
#Initialise
pos=0
for count in ${varList[@]}; 
do
  #echo "    ${varList[$pos]} = "
  echo "    uint "${varList[$pos]}"_v(trigNames.triggerIndex(("\"${varList[$pos]}"_v\"+string(buffer)).c_str()));"
  echo "    if("${varList[$pos]}"_v<triggerBits->size()) "${varList[$pos]}" = triggerBits->accept("${varList[$pos]}"_v);"  
  let pos=pos+1
done
echo " "
#Set branches
pos=0
for count in ${varList[@]}; 
do
  echo "  AddBranch(&${varList[$pos]}               ,\"${varList[$pos]}\");"
  let pos=pos+1
done
echo " "
#Set clear 
pos=0
for count in ${varList[@]}; 
do
  echo "  ${varList[$pos]} = -9999;"
  let pos=pos+1
done
echo " "
