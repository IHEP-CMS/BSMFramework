#!/bin/bash
#Specify needed variables
varType=double 
varList=(
    BJetnessValues_jetpt0
    BJetnessValues_numleps
    BJetnessValues_npvTrkOVcollTrk
    BJetnessValues_avip3dval
    BJetnessValues_avip3dsig
    BJetnessValues_avsip3dsig
    BJetnessValues_avip1dsig
)
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
  echo "    ${varList[$pos]} = "
  #echo "    uint "${varList[$pos]}"_v(trigNames.triggerIndex(("\"${varList[$pos]}"_v\"+string(buffer)).c_str()));"
  #echo "    if("${varList[$pos]}"_v<triggerBits->size()) "${varList[$pos]}" = triggerBits->accept("${varList[$pos]}"_v);"  
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
