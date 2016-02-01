#!/bin/bash
#Specify needed variables
varType=varC
vartype=double #double, int, float
obj=
varList=(charge cosdphi pzeta pzetavis masseff)
#Declare variables for reading
if [ "$varType" == varA ] || [ "$varType" == varB ]
then
  echo " "
  pos=0
  for count in $obj${varList[@]};
  do
    echo "  $vartype r$obj${varList[$pos]}; r$obj${varList[$pos]} = 0; TBranch* b_r$obj${varList[$pos]} = 0; readingtree->SetBranchAddress(\"$obj${varList[$pos]}\",&r$obj${varList[$pos]},&b_r$obj${varList[$pos]});"
  let pos=pos+1
  done
fi

if [ "$varType" == varA ] || [ "$varType" == varC ] 
then
  echo " "
  pos=0
  for count in $obj${varList[@]};
  do
   echo "  $vartype $obj${varList[$pos]}; newtree->Branch(\"$obj${varList[$pos]}\",&$obj${varList[$pos]});";
  let pos=pos+1
  done
fi 

#Read branches
if [ "$varType" == varA ] || [ "$varType" == varB ]
then 
  echo " "
  pos=0
  for count in $obj${varList[@]};
  do
    echo "   b_r$obj${varList[$pos]}->GetEntry(en);"
    let pos=pos+1
  done
fi

#Initialise new branches
if [ "$varType" == varA ] || [ "$varType" == varC ]
then
  echo " "
  pos=0
  for count in $obj${varList[@]};
  do
    echo "   $obj${varList[$pos]} = -999;"  
  let pos=pos+1
  done
fi

#Saving branches for varC variables
if [ "$varType" == varA ] || [ "$varType" == varC ]
then
 echo " "
 pos=0
 for count in $obj${varList[@]};
 do
  echo "   $obj${varList[$pos]} = ;"
  let pos=pos+1
 done
fi
