#!/bin/bash
#Specify needed variables
varType=varC #varA: RW, varB: R, varC: W
vartype=double #double, int, float
obj=
varList=(BJetnessFV_num_leps BJetnessFV_npvTrkOVcollTrk BJetnessFV_avip3d_val BJetnessFV_avip3d_sig BJetnessFV_avsip3d_sig BJetnessFV_avip1d_sig)
#Declare variables for reading
if [ "$varType" == varA ] || [ "$varType" == varC ]
then 
 echo " "
 pos=0
 for count in $obj${varList[@]};
 do
   echo "$vartype r$obj${varList[$pos]}; TBranch* b_r$obj${varList[$pos]} = 0;" 
 let pos=pos+1
 done
fi

#Initialize variables to be written
if [ "$varType" == varA ] || [ "$varType" == varC ] 
then
  echo " "
  pos=0
  for count in $obj${varList[@]};
  do
    echo "$vartype $obj${varList[$pos]};" 
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
   #echo " $obj${varList[$pos]} = r$obj${varList[$pos]};"
   echo " $obj${varList[$pos]} = r$obj${varList[$pos]}->at(0);"
   let pos=pos+1
  done
fi

#Variable handling
if [ "$varType" == varA ] || [ "$varType" == varB ]
then 
  echo " "
  pos=0
  for count in $obj${varList[@]};
  do
    echo " readingtree->SetBranchAddress(\"$obj${varList[$pos]}\",&r$obj${varList[$pos]},&b_r$obj${varList[$pos]});"
    let pos=pos+1
  done
fi

#Variable handling
if [ "$varType" == varA ] || [ "$varType" == varC ]
then 
  echo " "
  pos=0
  for count in $obj${varList[@]};
  do
    echo " newtree->Branch(\"$obj${varList[$pos]}\",&$obj${varList[$pos]});"
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
    echo " $obj${varList[$pos]} = -999;"  
  let pos=pos+1
  done
fi

#Get entry
if [ "$varType" == varA ] || [ "$varType" == varB ]
then 
  echo " "
  pos=0
  for count in $obj${varList[@]};
  do
    echo " b_r$obj${varList[$pos]}->GetEntry(tentry);"
    let pos=pos+1
  done
fi
