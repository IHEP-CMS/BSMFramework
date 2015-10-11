#!/bin/bash
#Specify needed variables
varType=double
varList=(BJetness_num_pdgid_eles BJetness_num_pdgid_mus BJetness_num_soft_leps BJetness_num_soft_eles BJetness_num_soft_mus BJetness_num_vetonoipnoiso_leps BJetness_num_vetonoipnoiso_eles BJetness_num_loosenoipnoiso_leps BJetness_num_loosenoipnoiso_eles BJetness_num_loose_mus)
varLast=BJetness_npvPtOVcollPt
varCount=p
#Print info
echo " "
#Declare variables
echo -e " vector<$varType> \c"
pos=0
for count in ${varList[@]}; 
do
  if [ "${varList[$pos]}" != "$varLast" ] 
  then
   echo -e "${varList[$pos]}, \c"
  else
   echo "${varList[$pos]};"
  fi
  let pos=pos+1
done
echo " "
#Initialise
pos=0
for count in ${varList[@]}; 
do
  echo " ${varList[$pos]}.push_back();"
  let pos=pos+1
done
echo " "
#Set branches
pos=0
for count in ${varList[@]}; 
do
  echo " AddBranch(&${varList[$pos]}               ,\"${varList[$pos]}\");"
  let pos=pos+1
done
echo " "
#Set clear 
pos=0
for count in ${varList[@]}; 
do
  echo " ${varList[$pos]}.clear();"
  let pos=pos+1
done
echo " "
