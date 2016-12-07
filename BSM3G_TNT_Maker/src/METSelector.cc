#include "BSMFramework/BSM3G_TNT_Maker/interface/METSelector.h"
METSelector::METSelector(std::string name, TTree* tree, bool debug, const pset& iConfig, edm::ConsumesCollector && ic):
  baseTree(name,tree,debug)
{
  if(debug) std::cout<<"in METSelector constructor"<<std::endl;
  mets_      = ic.consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("mets"));
  puppimets_ = ic.consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("metsPUPPI"));
  if(debug) std::cout<<"in pileup constructor: calling SetBrances()"<<std::endl;
  _is_data   = iConfig.getParameter<bool>("is_data");
  _super_TNT = iConfig.getParameter<bool>("super_TNT");
  _MiniAODv2 = iConfig.getParameter<bool>("MiniAODv2");
  _PuppiVar = iConfig.getParameter<bool>("PuppiVar");
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
  iEvent.getByToken(mets_, mets);
  edm::Handle<pat::METCollection> puppimets;
  iEvent.getByToken(puppimets_, puppimets);
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
  if(!_is_data){
    Gen_type1PF_Met    = met.genMET()->pt();
    Gen_type1PF_Metpx  = met.genMET()->px();
    Gen_type1PF_Metpy  = met.genMET()->py();
    Gen_type1PF_Metpz  = met.genMET()->pz();
    Gen_type1PF_Meteta = met.genMET()->eta();
    Gen_type1PF_Metphi = met.genMET()->phi();
    Gen_type1PF_Meten  = met.genMET()->energy();
  }
  /////
  //   TTH variables
  /////
  //cout<<met.pt()<<setw(20)<<met.phi()<<endl;
  ////slimmedMETsPUPPI
  if(_PuppiVar){
    const pat::MET &puppimet = puppimets->front();
    //Kinematic
    Met_puppi_pt = puppimet.pt();
    Met_puppi_px = puppimet.px();
    Met_puppi_py = puppimet.py();
    Met_puppi_pz = puppimet.pz();
    Met_puppi_phi   = puppimet.phi();
    Met_puppi_sumEt = puppimet.sumEt();  
    //Corrections/Systematics
    if(_MiniAODv2){
      Met_puppi_shiftedPtUp   = puppimet.shiftedPt(pat::MET::JetEnUp);
      Met_puppi_shiftedPtDown = puppimet.shiftedPt(pat::MET::JetEnDown);
    }else{
      Met_puppi_shiftedPtUp   = -999;
      Met_puppi_shiftedPtDown = -999;
    }
    //MC
    if(!_is_data){
      Gen_puppi_Met    = puppimet.genMET()->pt();
      Gen_puppi_Metpx  = puppimet.genMET()->px();
      Gen_puppi_Metpy  = puppimet.genMET()->py();
      Gen_puppi_Metpz  = puppimet.genMET()->pz();
      Gen_puppi_Meteta = puppimet.genMET()->eta();
      Gen_puppi_Metphi = puppimet.genMET()->phi();
      Gen_puppi_Meten  = puppimet.genMET()->energy();
    }
  }
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
  AddBranch(&Gen_type1PF_Metpx         ,"Gen_type1PF_Metpx");
  AddBranch(&Gen_type1PF_Metpy         ,"Gen_type1PF_Metpy");
  AddBranch(&Gen_type1PF_Metpz         ,"Gen_type1PF_Metpz");
  AddBranch(&Gen_type1PF_Meteta        ,"Gen_type1PF_Meteta");
  AddBranch(&Gen_type1PF_Metphi        ,"Gen_type1PF_Metphi");
  AddBranch(&Gen_type1PF_Meten         ,"Gen_type1PF_Meten");
  ////slimmedMETsPUPPI
  if(_PuppiVar){
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
    AddBranch(&Gen_puppi_Metpx         ,"Gen_puppi_Metpx");
    AddBranch(&Gen_puppi_Metpy         ,"Gen_puppi_Metpy");
    AddBranch(&Gen_puppi_Metpz         ,"Gen_puppi_Metpz");
    AddBranch(&Gen_puppi_Meteta        ,"Gen_puppi_Meteta");
    AddBranch(&Gen_puppi_Metphi        ,"Gen_puppi_Metphi");
    AddBranch(&Gen_puppi_Meten         ,"Gen_puppi_Meten");
  }
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
  Gen_type1PF_Metpx         = -9999;
  Gen_type1PF_Metpy         = -9999;
  Gen_type1PF_Metpz         = -9999;
  Gen_type1PF_Meteta        = -9999;
  Gen_type1PF_Metphi        = -9999;
  Gen_type1PF_Meten         = -9999;
  ////slimmedMETsPUPPI
  if(_PuppiVar){
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
    Gen_puppi_Metpx         = -9999;
    Gen_puppi_Metpy         = -9999;
    Gen_puppi_Metpz         = -9999;
    Gen_puppi_Meteta        = -9999;
    Gen_puppi_Metphi        = -9999;
    Gen_puppi_Meten         = -9999;
  }
}
