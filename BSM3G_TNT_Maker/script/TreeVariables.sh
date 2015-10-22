#!/bin/bash
#Specify needed variables
#int
#varType=int
#varTYPE=INT
#capLetter=I
#varList=(rtau_pt rtau_eta rtau_phi rtau_en)
#varLast=jet_chtrksnpvtt
#varNum=jn
#varCount=jn
#double
obj=Muon_
varType=double
varTYPE=DOUBLE
capLetter=D
varList=(p)
varLast=
varNum=ele_num
varCount=gl
#Print info
echo " "
#Declare variables
echo -e " $varType \c"
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
  echo "  $obj${varList[$pos]} = DEF_VAL_$varTYPE;"
  let pos=pos+1
done
echo " "
#Set branches
pos=0
for count in ${varList[@]}; 
do
  echo "  tree->Branch(\"$obj${varList[$pos]}\", &$obj${varList[$pos]}, \"$obj${varList[$pos]}/$capLetter\");"
  let pos=pos+1
done
echo " "
#Set branch address
pos=0
for count in ${varList[@]}; 
do
  echo "  tree->SetBranchAddress(\"$obj${varList[$pos]}\", &$obj${varList[$pos]});"
  let pos=pos+1
done
echo " "
#Analyzer
pos=0
for count in ${varList[@]}; 
do
  echo "tree->$obj${varList[$pos]} = ;"
  let pos=pos+1
done
echo " "
#pftauchhads_AEIP1D_x_val_trk0 pftauchhads_AEIP1D_y_val_trk0 pftauchhads_IP3DvalAE_trk0 pftauchhads_IP3DvalTE_trk0 pftauchhads_IP2DvalAE_trk0 pftauchhads_IP2DvalTE_trk0 pftauchhads_IP1DvalAE_trk0 pftauchhads_IP1DvalTE_trk0)
#dR_tautrk_recojettau_trk0 dR_tautrk_recojettau_trk1 dR_tautrk_recojettau_trk2)
#tau_pt_DIV_recojettau_pt tau_pt_DIV_recojettau_en tau_trk0_pt_DIV_recojettau_pt tau_trk1_pt_DIV_recojettau_pt tau_trk2_pt_DIV_recojettau_pt tau_trk0_pt_DIV_recojettau_en tau_trk1_pt_DIV_recojettau_en tau_trk2_pt_DIV_recojettau_en)
#pftauchhads_AEIP1D_val_trk0 pftauchhads_AEIP1D_val_trk1 pftauchhads_AEIP1D_val_trk2 pftauchhads_AEIP1D_sig_trk0 pftauchhads_AEIP1D_sig_trk1 pftauchhads_AEIP1D_sig_trk2 pftauchhads_sAEIP1D_val_trk0 pftauchhads_sAEIP1D_val_trk1 pftauchhads_sAEIP1D_val_trk2 pftauchhads_sAEIP1D_sig_trk0 pftauchhads_sAEIP1D_sig_trk1 pftauchhads_sAEIP1D_sig_trk2)
#pv_avfbs_nchi2 unbpv_avfbs_nchi2 diff_pv_avfbs_nchi2_unbpv_avfbs_nchi2)
#pvsv_dist3d_val pvsv_dist2d_val pvsv_dist3d_sig pvsv_dist2d_sig)
#recojettau_pt recojettau_eta recojettau_phi recojettau_en)
#unbpv_KVF_nv unbpv_KVFbs_nv unbpv_AVF_nv unbpv_AVFbs_nv)
#vtxKVF_x vtxKVF_y vtxKVF_z vtxKVFbs_x vtxKVFbs_y vtxKVFbs_z vtxAVF_x vtxAVF_y vtxAVF_z vtxAVFbs_x vtxAVFbs_y vtxAVFbs_z)
#genditauvtx_x genditauvtx_y genditauvtx_z)
#pvtauvtx_genditauvtx_x pvtauvtx_genditauvtx_y pvtauvtx_genditauvtx_z unbpv_genditauvtx_x unbpv_genditauvtx_y unbpv_genditauvtx_z)

