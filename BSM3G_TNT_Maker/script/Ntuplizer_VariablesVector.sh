#!/bin/bash
#Specify needed variables
obj=
varType=double
varList=(
Tau_packedLeadTauCand_dxy Tau_packedLeadTauCand_dz
)
varLast=
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
   echo -e "$obj${varList[$pos]}, \c"
  else
   echo "$obj${varList[$pos]};"
  fi
  let pos=pos+1
done
echo " "
#Initialise
pos=0
for count in ${varList[@]}; 
do
  echo "    $obj${varList[$pos]}.push_back();" #-999);"
  let pos=pos+1
done
echo " "
#Set branches
pos=0
for count in ${varList[@]}; 
do
  echo "  AddBranch(&$obj${varList[$pos]}               ,\"$obj${varList[$pos]}\");"
  let pos=pos+1
done
echo " "
#Set clear 
pos=0
for count in ${varList[@]}; 
do
  echo "  $obj${varList[$pos]}.clear();"
  let pos=pos+1
done
echo " "
