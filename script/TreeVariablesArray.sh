#!/bin/bash
#Specify needed variables
#int
#varType=int
#varTYPE=INT
#capLetter=I
#varList=(jet_sumndaus jet_sumchtrks jet_sumchtrkspv jet_sumchtrksnpv jet_sumchtrkspvtt jet_sumchtrksnpvtt)
#varLast=jet_sumchtrksnpvtt
#varNum=jet_num
#varCount=jct
#double
varType=double
varTYPE=DOUBLE
capLetter=D
varList=(pftauchhads_IP3Ddv_val pftauchhads_IP3Ddv_err pftauchhads_IP3Ddv_sig pftauchhads_IP2Ddv_val pftauchhads_IP2Ddv_err pftauchhads_IP2Ddv_sig pftauchhads_sIP3Ddv_val pftauchhads_sIP3Ddv_err pftauchhads_sIP3Ddv_sig pftauchhads_sIP2Ddv_val pftauchhads_sIP2Ddv_err pftauchhads_sIP2Ddv_sig)
#dR_recojettau_recotau dR_tautrk_recotau)
#(pftauchhads_IP1D_val pftauchhads_IP1D_err pftauchhads_IP1D_sig pftauchhads_sIP1D_val pftauchhads_sIP1D_err pftauchhads_sIP1D_sig)
#pftauchhads_pt pftauchhads_eta pftauchhads_phi pftauchhads_en pftauchhads_IP3D_val pftauchhads_IP3D_err pftauchhads_IP3D_sig pftauchhads_IP2D_val pftauchhads_IP2D_err pftauchhads_IP2D_sig pftauchhads_sIP3D_val pftauchhads_sIP3D_err pftauchhads_sIP3D_sig pftauchhads_sIP2D_val pftauchhads_sIP2D_err pftauchhads_sIP2D_sig)
varLast=pftauchhads_sIP2Ddv_sig
varNum='pftauchhads_numv'
varCount=p
#Print info
echo " "
#Declare variables
echo -e " $varType \c"
pos=0
for count in ${varList[@]}; 
do
  if [ "${varList[$pos]}" != "$varLast" ] 
  then
   echo -e "${varList[$pos]}[DEF_SIZE1D], \c"
  else
   echo "${varList[$pos]}[DEF_SIZE1D];"
  fi
  let pos=pos+1
done
echo " "
#Initialise
pos=0
for count in ${varList[@]}; 
do
  echo "  INIT_1DARRAY(${varList[$pos]},DEF_SIZE1D,DEF_VAL_$varTYPE);"
  let pos=pos+1
done
echo " "
#Set branches
pos=0
for count in ${varList[@]}; 
do
  echo "  tree->Branch(\"${varList[$pos]}\", &${varList[$pos]}, \"${varList[$pos]}[$varNum]/$capLetter\");"
  let pos=pos+1
done
echo " "
#Set branch address
pos=0
for count in ${varList[@]}; 
do
  echo "  tree->SetBranchAddress(\"${varList[$pos]}\", &${varList[$pos]});"
  let pos=pos+1
done
echo " "
#Analyzer
pos=0
for count in ${varList[@]}; 
do
  echo "tree->${varList[$pos]}[$varCount] = ;"
  let pos=pos+1
done
echo " "
