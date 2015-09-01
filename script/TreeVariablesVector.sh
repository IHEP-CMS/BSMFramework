#!/bin/bash
#Specify needed variables
varType=double
varList=(Muon_isTracker Muon_isMedium)
varLast=Muon_eta
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
