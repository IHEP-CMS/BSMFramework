#include "BSMFramework/BSM3G_TNT_Maker/interface/METSelector.h"
METSelector::METSelector(std::string name, TTree* tree, bool debug, const pset& iConfig):baseTree(name,tree,debug){
  if(debug) std::cout<<"in METSelector constructor"<<std::endl;
  metToken_ = iConfig.getParameter<edm::InputTag>("mets");
  puppi_metToken_ = iConfig.getParameter<edm::InputTag>("metsPUPPI");
  if(debug) std::cout<<"in pileup constructor: calling SetBrances()"<<std::endl;
  _is_data   = iConfig.getParameter<bool>("is_data");
  _super_TNT = iConfig.getParameter<bool>("super_TNT");
  SetBranches();
}
METSelector::~METSelector(){
  delete tree_;
}
void METSelector::Fill(const edm::Event& iEvent){
  if(debug_) std::cout<<"getting met info"<<std::endl;
  Clear();
  /////
  //   Recall collections
  /////  
  edm::Handle<pat::METCollection> mets;
  iEvent.getByLabel(metToken_, mets);
  edm::Handle<pat::METCollection> puppimets;
  iEvent.getByLabel(puppi_metToken_, puppimets);
  if(debug_) std::cout<<"Filling met branches"<<std::endl;
  /////
  //   Get muon information
  ///// 
  ////slimmedMETs
  const pat::MET &met = mets->front();
  //Kinematic
  Met_type1PF_pt = met.pt();
  Met_type1PF_px = met.px();
  Met_type1PF_py = met.py();
  Met_type1PF_pz = met.pz();
  Met_type1PF_phi   = met.phi();
  Met_type1PF_sumEt = met.sumEt();  
  //Corrections/Systematics
  Met_type1PF_shiftedPtUp   = met.shiftedPt(pat::MET::JetEnUp);
  Met_type1PF_shiftedPtDown = met.shiftedPt(pat::MET::JetEnDown);
  //MC
  if(!_is_data) Gen_type1PF_Met = met.genMET()->pt();
  ////slimmedMETsPUPPI
  const pat::MET &puppimet = puppimets->front();
  //Kinematic
  Met_puppi_pt = puppimet.pt();
  Met_puppi_px = puppimet.px();
  Met_puppi_py = puppimet.py();
  Met_puppi_pz = puppimet.pz();
  Met_puppi_phi   = puppimet.phi();
  Met_puppi_sumEt = puppimet.sumEt();  
  //Corrections/Systematics
  Met_puppi_shiftedPtUp   = puppimet.shiftedPt(pat::MET::JetEnUp);
  Met_puppi_shiftedPtDown = puppimet.shiftedPt(pat::MET::JetEnDown);
  //MC
  if(!_is_data) Gen_puppi_Met = puppimet.genMET()->pt();
  if(debug_) std::cout<<"got MET info"<<std::endl;
}
void METSelector::SetBranches(){
  if(debug_) std::cout<<"setting branches: calling AddBranch of baseTree"<<std::endl;
  ////slimmedMETs
  //Kinematic  
  AddBranch(&Met_type1PF_pt,            "Met_type1PF_pt");
  AddBranch(&Met_type1PF_px,            "Met_type1PF_px");
  AddBranch(&Met_type1PF_py,            "Met_type1PF_py");
  AddBranch(&Met_type1PF_pz,            "Met_type1PF_pz");
  AddBranch(&Met_type1PF_phi,           "Met_type1PF_phi");
  AddBranch(&Met_type1PF_sumEt,         "Met_type1PF_sumEt");
  //Corrections/Systematics
  AddBranch(&Met_type1PF_shiftedPtUp,   "Met_type1PF_shiftedPtUp");
  AddBranch(&Met_type1PF_shiftedPtDown, "Met_type1PF_shiftedPtDown");
  //MC
  AddBranch(&Gen_type1PF_Met,           "Gen_type1PF_Met");
  ////slimmedMETsPUPPI
  //Kinematic  
  AddBranch(&Met_puppi_pt,            "Met_puppi_pt");
  AddBranch(&Met_puppi_px,            "Met_puppi_px");
  AddBranch(&Met_puppi_py,            "Met_puppi_py");
  AddBranch(&Met_puppi_pz,            "Met_puppi_pz");
  AddBranch(&Met_puppi_phi,           "Met_puppi_phi");
  AddBranch(&Met_puppi_sumEt,         "Met_puppi_sumEt");
  //Corrections/Systematics
  AddBranch(&Met_puppi_shiftedPtUp,   "Met_puppi_shiftedPtUp");
  AddBranch(&Met_puppi_shiftedPtDown, "Met_puppi_shiftedPtDown");
  //MC
  AddBranch(&Gen_puppi_Met,           "Gen_puppi_Met");
  if(debug_) std::cout<<"set branches"<<std::endl;
}
void METSelector::Clear(){
  ////slimmedMETs
  //Kinematic  
  Met_type1PF_pt            = -9999;
  Met_type1PF_px            = -9999;
  Met_type1PF_py            = -9999;
  Met_type1PF_pz            = -9999;
  Met_type1PF_phi           = -9999;
  Met_type1PF_sumEt         = -9999; 
  //Corrections/Systematics
  Met_type1PF_shiftedPtUp   = -9999;
  Met_type1PF_shiftedPtDown = -9999;
  //MC
  Gen_type1PF_Met           = -9999;
  ////slimmedMETsPUPPI
  //Kinematic  
  Met_puppi_pt            = -9999;
  Met_puppi_px            = -9999;
  Met_puppi_py            = -9999;
  Met_puppi_pz            = -9999;
  Met_puppi_phi           = -9999;
  Met_puppi_sumEt         = -9999; 
  //Corrections/Systematics
  Met_puppi_shiftedPtUp   = -9999;
  Met_puppi_shiftedPtDown = -9999;
  //MC
  Gen_puppi_Met           = -9999;
}
