#!/bin/bash
#Specify needed variables
vartype=double #double, int, float
obj=  #rMet_type1PF_pt #Tau_, Jet_, Muon_
ob= #met
varList=(BJetnessFV_num_leps BJetnessFV_npvTrkOVcollTrk BJetnessFV_avip3d_val BJetnessFV_avip3d_sig BJetnessFV_avsip3d_sig BJetnessFV_avip1d_sig)
#_Declare variables for reading
echo " "
pos=0
for count in ${varList[@]};
do
 echo "vector<$vartype>* r$obj${varList[$pos]}; TBranch* b_r$obj${varList[$pos]} = 0;"
 let pos=pos+1
done

#Declare variables for writing
echo " "
pos=0
for count in ${varList[@]};
do
 echo "vector<$vartype>* $obj${varList[$pos]} = new std::vector<$vartype>;"
 let pos=pos+1
done

#Save new branches content
echo " "
pos=0
for count in ${varList[@]};
do
 echo " $obj${varList[$pos]}->push_back(r$obj${varList[$pos]}->at(0));"
 let pos=pos+1
done

#Reading tree
echo " "
pos=0
for count in ${varList[@]};
do
 echo " readingtree->SetBranchAddress(\"$obj${varList[$pos]}\",&r$obj${varList[$pos]},&b_r$obj${varList[$pos]});"
 let pos=pos+1
done

#New tree
echo " "
pos=0
for count in ${varList[@]};
do
 echo " newtree->Branch(\"$obj${varList[$pos]}\",&$obj${varList[$pos]});"
 let pos=pos+1
done

#Clear new branches
echo " "
pos=0
for count in ${varList[@]};
do
 echo " $obj${varList[$pos]}->clear();"  
 let pos=pos+1
done

#Read branches
echo " "
pos=0
for count in ${varList[@]};
do
 echo " b_r$obj${varList[$pos]}->GetEntry(tentry);"
 let pos=pos+1
done
