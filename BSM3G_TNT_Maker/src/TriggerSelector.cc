#include "BSMFramework/BSM3G_TNT_Maker/interface/TriggerSelector.h"
#include <string>       // std::string
#include <iostream>     // std::cout
#include <sstream>      // std::stringstream, std::stringbuf
TriggerSelector::TriggerSelector(std::string name, TTree* tree, bool debug, const pset& iConfig):baseTree(name,tree,debug){
  triggerBits_       = iConfig.getParameter<edm::InputTag>("bits");
  _maxtriggerversion = iConfig.getParameter<double>("maxtriggerversion");
  _reduced = iConfig.getParameter<bool>("reduced");
  SetBranches();
}
TriggerSelector::~TriggerSelector(){
  delete tree_;
}
void TriggerSelector::Fill(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  if(debug_)    std::cout<<"getting met info"<<std::endl;
  Clear();
  //Trigget paths  
  edm::Handle<edm::TriggerResults> triggerBits;
  iEvent.getByLabel(triggerBits_, triggerBits);
  const edm::TriggerNames &trigNames = iEvent.triggerNames(*triggerBits);
  for(double tv=0.; tv<=_maxtriggerversion; tv++){ 
    char buffer[14]; sprintf(buffer,"%g",tv);
    uint HLT_DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg_v(trigNames.triggerIndex(("HLT_DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg_v"+string(buffer)).c_str()));
    if(HLT_DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg_v<triggerBits->size()) HLT_DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg = triggerBits->accept(HLT_DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg_v);
    uint HLT_DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_v(trigNames.triggerIndex(("HLT_DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_v"+string(buffer)).c_str()));
    if(HLT_DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_v<triggerBits->size()) HLT_DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg = triggerBits->accept(HLT_DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_v);
    uint HLT_DoubleTightChargedIsoPFTau40_Trk1_eta2p1_Reg_v(trigNames.triggerIndex(("HLT_DoubleTightChargedIsoPFTau40_Trk1_eta2p1_Reg_v"+string(buffer)).c_str()));
    if(HLT_DoubleTightChargedIsoPFTau40_Trk1_eta2p1_Reg_v<triggerBits->size()) HLT_DoubleTightChargedIsoPFTau40_Trk1_eta2p1_Reg = triggerBits->accept(HLT_DoubleTightChargedIsoPFTau40_Trk1_eta2p1_Reg_v);
    uint HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1_v(trigNames.triggerIndex(("HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1_v"+string(buffer)).c_str()));
    if(HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1_v<triggerBits->size()) HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1 = triggerBits->accept(HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1_v);
    uint HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTau30_eta2p1_CrossL1_v(trigNames.triggerIndex(("HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTau30_eta2p1_CrossL1_v"+string(buffer)).c_str()));
    if(HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTau30_eta2p1_CrossL1_v<triggerBits->size()) HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTau30_eta2p1_CrossL1 = triggerBits->accept(HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTau30_eta2p1_CrossL1_v);
    uint HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTau30_eta2p1_CrossL1_v(trigNames.triggerIndex(("HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTau30_eta2p1_CrossL1_v"+string(buffer)).c_str()));
    if(HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTau30_eta2p1_CrossL1_v<triggerBits->size()) HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTau30_eta2p1_CrossL1 = triggerBits->accept(HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTau30_eta2p1_CrossL1_v);
    uint HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1_v(trigNames.triggerIndex(("HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1_v"+string(buffer)).c_str()));
    if(HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1_v<triggerBits->size()) HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1 = triggerBits->accept(HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1_v);
    uint HLT_IsoMu20_eta2p1_MediumChargedIsoPFTau27_eta2p1_CrossL1_v(trigNames.triggerIndex(("HLT_IsoMu20_eta2p1_MediumChargedIsoPFTau27_eta2p1_CrossL1_v"+string(buffer)).c_str()));
    if(HLT_IsoMu20_eta2p1_MediumChargedIsoPFTau27_eta2p1_CrossL1_v<triggerBits->size()) HLT_IsoMu20_eta2p1_MediumChargedIsoPFTau27_eta2p1_CrossL1 = triggerBits->accept(HLT_IsoMu20_eta2p1_MediumChargedIsoPFTau27_eta2p1_CrossL1_v);
    uint HLT_IsoMu20_eta2p1_TightChargedIsoPFTau27_eta2p1_CrossL1_v(trigNames.triggerIndex(("HLT_IsoMu20_eta2p1_TightChargedIsoPFTau27_eta2p1_CrossL1_v"+string(buffer)).c_str()));
    if(HLT_IsoMu20_eta2p1_TightChargedIsoPFTau27_eta2p1_CrossL1_v<triggerBits->size()) HLT_IsoMu20_eta2p1_TightChargedIsoPFTau27_eta2p1_CrossL1 = triggerBits->accept(HLT_IsoMu20_eta2p1_TightChargedIsoPFTau27_eta2p1_CrossL1_v);
    uint HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_1pr_v(trigNames.triggerIndex(("HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_1pr_v"+string(buffer)).c_str()));
    if(HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_1pr_v<triggerBits->size()) HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_1pr = triggerBits->accept(HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_1pr_v);
    uint HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_v(trigNames.triggerIndex(("HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_v"+string(buffer)).c_str()));
    if(HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_v<triggerBits->size()) HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1 = triggerBits->accept(HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_v);
    uint HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET100_v(trigNames.triggerIndex(("HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET100_v"+string(buffer)).c_str()));
    if(HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET100_v<triggerBits->size()) HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET100 = triggerBits->accept(HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET100_v);
    uint HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET90_v(trigNames.triggerIndex(("HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET90_v"+string(buffer)).c_str()));
    if(HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET90_v<triggerBits->size()) HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET90 = triggerBits->accept(HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET90_v);
    uint HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_v(trigNames.triggerIndex(("HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_v"+string(buffer)).c_str()));
    if(HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_v<triggerBits->size()) HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr = triggerBits->accept(HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_v);
    uint HLT_PFHT500_PFMET100_PFMHT100_IDTight_v(trigNames.triggerIndex(("HLT_PFHT500_PFMET100_PFMHT100_IDTight_v"+string(buffer)).c_str()));
    if(HLT_PFHT500_PFMET100_PFMHT100_IDTight_v<triggerBits->size()) HLT_PFHT500_PFMET100_PFMHT100_IDTight = triggerBits->accept(HLT_PFHT500_PFMET100_PFMHT100_IDTight_v);
    uint HLT_PFHT700_PFMET85_PFMHT85_IDTight_v(trigNames.triggerIndex(("HLT_PFHT700_PFMET85_PFMHT85_IDTight_v"+string(buffer)).c_str()));
    if(HLT_PFHT700_PFMET85_PFMHT85_IDTight_v<triggerBits->size()) HLT_PFHT700_PFMET85_PFMHT85_IDTight = triggerBits->accept(HLT_PFHT700_PFMET85_PFMHT85_IDTight_v);
    uint HLT_PFHT800_PFMET75_PFMHT75_IDTight_v(trigNames.triggerIndex(("HLT_PFHT800_PFMET75_PFMHT75_IDTight_v"+string(buffer)).c_str()));
    if(HLT_PFHT800_PFMET75_PFMHT75_IDTight_v<triggerBits->size()) HLT_PFHT800_PFMET75_PFMHT75_IDTight = triggerBits->accept(HLT_PFHT800_PFMET75_PFMHT75_IDTight_v);
    uint HLT_CaloJet500_NoJetID_v(trigNames.triggerIndex(("HLT_CaloJet500_NoJetID_v"+string(buffer)).c_str()));
    if(HLT_CaloJet500_NoJetID_v<triggerBits->size()) HLT_CaloJet500_NoJetID = triggerBits->accept(HLT_CaloJet500_NoJetID_v);
    uint HLT_PFHT1050_v(trigNames.triggerIndex(("HLT_PFHT1050_v"+string(buffer)).c_str()));
    if(HLT_PFHT1050_v<triggerBits->size()) HLT_PFHT1050 = triggerBits->accept(HLT_PFHT1050_v);
    uint HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight_v(trigNames.triggerIndex(("HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight"+string(buffer)).c_str()));
    if(HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight_v<triggerBits->size()) HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight = triggerBits->accept(HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight_v);
    uint HLT_PFMET120_PFMHT120_IDTight_v(trigNames.triggerIndex(("HLT_PFMET120_PFMHT120_IDTight_v"+string(buffer)).c_str()));
    if(HLT_PFMET120_PFMHT120_IDTight_v<triggerBits->size()) HLT_PFMET120_PFMHT120_IDTight = triggerBits->accept(HLT_PFMET120_PFMHT120_IDTight_v);
    uint HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v(trigNames.triggerIndex(("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v"+string(buffer)).c_str()));
    if(HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v<triggerBits->size()) HLT_PFMETNoMu120_PFMHTNoMu120_IDTight = triggerBits->accept(HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v);
    uint HLT_PFJet500_v(trigNames.triggerIndex(("HLT_PFJet500_v"+string(buffer)).c_str()));
    if(HLT_PFJet500_v<triggerBits->size()) HLT_PFJet500 = triggerBits->accept(HLT_PFJet500_v);
    uint HLT_AK8PFJet500_v(trigNames.triggerIndex(("HLT_AK8PFJet500_v"+string(buffer)).c_str()));
    if(HLT_AK8PFJet500_v<triggerBits->size()) HLT_AK8PFJet500 = triggerBits->accept(HLT_AK8PFJet500_v);
    uint HLT_Ele27_eta2p1_WPTight_Gsf_v(trigNames.triggerIndex(("HLT_Ele27_eta2p1_WPTight_Gsf_v"+string(buffer)).c_str()));
    if(HLT_Ele27_eta2p1_WPTight_Gsf_v<triggerBits->size()) HLT_Ele27_eta2p1_WPTight_Gsf = triggerBits->accept(HLT_Ele27_eta2p1_WPTight_Gsf_v);
    uint HLT_Ele27_WPTight_Gsf_v(trigNames.triggerIndex(("HLT_Ele27_WPTight_Gsf_v"+string(buffer)).c_str()));
    if(HLT_Ele27_WPTight_Gsf_v<triggerBits->size()) HLT_Ele27_WPTight_Gsf = triggerBits->accept(HLT_Ele27_WPTight_Gsf_v);
    uint HLT_Ele25_eta2p1_WPTight_Gsf_v(trigNames.triggerIndex(("HLT_Ele25_eta2p1_WPTight_Gsf_v"+string(buffer)).c_str()));
    if(HLT_Ele25_eta2p1_WPTight_Gsf_v<triggerBits->size()) HLT_Ele25_eta2p1_WPTight_Gsf = triggerBits->accept(HLT_Ele25_eta2p1_WPTight_Gsf_v);
    uint HLT_Ele115_CaloIdVT_GsfTrkIdT_v(trigNames.triggerIndex(("HLT_Ele115_CaloIdVT_GsfTrkIdT_v"+string(buffer)).c_str()));
    if(HLT_Ele115_CaloIdVT_GsfTrkIdT_v<triggerBits->size()) HLT_Ele115_CaloIdVT_GsfTrkIdT = triggerBits->accept(HLT_Ele115_CaloIdVT_GsfTrkIdT_v);
    uint HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v(trigNames.triggerIndex(("HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v"+string(buffer)).c_str()));
    if(HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v<triggerBits->size()) HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf = triggerBits->accept(HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v);
    uint HLT_DoubleEle33_CaloIdL_MW_v(trigNames.triggerIndex(("HLT_DoubleEle33_CaloIdL_MW_v"+string(buffer)).c_str()));
    if(HLT_DoubleEle33_CaloIdL_MW_v<triggerBits->size()) HLT_DoubleEle33_CaloIdL_MW = triggerBits->accept(HLT_DoubleEle33_CaloIdL_MW_v);
    uint HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v(trigNames.triggerIndex(("HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v"+string(buffer)).c_str()));
    if(HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v<triggerBits->size()) HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW = triggerBits->accept(HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v);
    uint HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v(trigNames.triggerIndex(("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"+string(buffer)).c_str()));
    if(HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v<triggerBits->size()) HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ = triggerBits->accept(HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v);
    uint HLT_IsoMu22_v(trigNames.triggerIndex(("HLT_IsoMu22_v"+string(buffer)).c_str()));
    if(HLT_IsoMu22_v<triggerBits->size()) HLT_IsoMu22 = triggerBits->accept(HLT_IsoMu22_v);
    uint HLT_IsoTkMu22_v(trigNames.triggerIndex(("HLT_IsoTkMu22_v"+string(buffer)).c_str()));
    if(HLT_IsoTkMu22_v<triggerBits->size()) HLT_IsoTkMu22 = triggerBits->accept(HLT_IsoTkMu22_v);
    uint HLT_IsoMu24_v(trigNames.triggerIndex(("HLT_IsoMu24_v"+string(buffer)).c_str()));
    if(HLT_IsoMu24_v<triggerBits->size()) HLT_IsoMu24 = triggerBits->accept(HLT_IsoMu24_v);
    uint HLT_IsoTkMu24_v(trigNames.triggerIndex(("HLT_IsoTkMu24_v"+string(buffer)).c_str()));
    if(HLT_IsoTkMu24_v<triggerBits->size()) HLT_IsoTkMu24 = triggerBits->accept(HLT_IsoTkMu24_v);
    uint HLT_IsoMu22_eta2p1_v(trigNames.triggerIndex(("HLT_IsoMu22_eta2p1_v"+string(buffer)).c_str()));
    if(HLT_IsoMu22_eta2p1_v<triggerBits->size()) HLT_IsoMu22_eta2p1 = triggerBits->accept(HLT_IsoMu22_eta2p1_v);
    uint HLT_IsoTkMu22_eta2p1_v(trigNames.triggerIndex(("HLT_IsoTkMu22_eta2p1_v"+string(buffer)).c_str()));
    if(HLT_IsoTkMu22_eta2p1_v<triggerBits->size()) HLT_IsoTkMu22_eta2p1 = triggerBits->accept(HLT_IsoTkMu22_eta2p1_v);
    uint HLT_Mu50_v(trigNames.triggerIndex(("HLT_Mu50_v"+string(buffer)).c_str()));
    if(HLT_Mu50_v<triggerBits->size()) HLT_Mu50 = triggerBits->accept(HLT_Mu50_v);
    uint HLT_TkMu50_v(trigNames.triggerIndex(("HLT_TkMu50_v"+string(buffer)).c_str()));
    if(HLT_TkMu50_v<triggerBits->size()) HLT_TkMu50 = triggerBits->accept(HLT_TkMu50_v);
    uint HLT_DoubleMu33NoFiltersNoVtx_v(trigNames.triggerIndex(("HLT_DoubleMu33NoFiltersNoVtx_v"+string(buffer)).c_str()));
    if(HLT_DoubleMu33NoFiltersNoVtx_v<triggerBits->size()) HLT_DoubleMu33NoFiltersNoVtx = triggerBits->accept(HLT_DoubleMu33NoFiltersNoVtx_v);
    uint HLT_DoubleMu23NoFiltersNoVtxDisplaced_v(trigNames.triggerIndex(("HLT_DoubleMu23NoFiltersNoVtxDisplaced_v"+string(buffer)).c_str()));
    if(HLT_DoubleMu23NoFiltersNoVtxDisplaced_v<triggerBits->size()) HLT_DoubleMu23NoFiltersNoVtxDisplaced = triggerBits->accept(HLT_DoubleMu23NoFiltersNoVtxDisplaced_v);
    uint HLT_Mu30_TkMu11_v(trigNames.triggerIndex(("HLT_Mu30_TkMu11_v"+string(buffer)).c_str()));
    if(HLT_Mu30_TkMu11_v<triggerBits->size()) HLT_Mu30_TkMu11 = triggerBits->accept(HLT_Mu30_TkMu11_v);
    uint HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v(trigNames.triggerIndex(("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v"+string(buffer)).c_str()));
    if(HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v<triggerBits->size()) HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ = triggerBits->accept(HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v);
  }
}

void TriggerSelector::SetBranches(){
  if(debug_)    std::cout<<"setting branches: calling AddBranch of baseTree"<<std::endl;
  AddBranch(&HLT_PFMET120_PFMHT120_IDTight             ,"HLT_PFMET120_PFMHT120_IDTight");
  AddBranch(&HLT_PFMETNoMu120_PFMHTNoMu120_IDTight     ,"HLT_PFMETNoMu120_PFMHTNoMu120_IDTight");
  AddBranch(&HLT_Mu50				      ,"HLT_Mu50");
  AddBranch(&HLT_TkMu50				      ,"HLT_TkMu50");
  AddBranch(&HLT_DoubleEle33_CaloIdL_MW		      ,"HLT_DoubleEle33_CaloIdL_MW");
  if(!_reduced){
    AddBranch(&HLT_DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg           ,"HLT_DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg");
    AddBranch(&HLT_DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg            ,"HLT_DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg");
    AddBranch(&HLT_DoubleTightChargedIsoPFTau40_Trk1_eta2p1_Reg                    ,"HLT_DoubleTightChargedIsoPFTau40_Trk1_eta2p1_Reg");
    AddBranch(&HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1  ,"HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1");
    AddBranch(&HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTau30_eta2p1_CrossL1 ,"HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTau30_eta2p1_CrossL1");
    AddBranch(&HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTau30_eta2p1_CrossL1  ,"HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTau30_eta2p1_CrossL1");
    AddBranch(&HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1            ,"HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1");
    AddBranch(&HLT_IsoMu20_eta2p1_MediumChargedIsoPFTau27_eta2p1_CrossL1           ,"HLT_IsoMu20_eta2p1_MediumChargedIsoPFTau27_eta2p1_CrossL1");
    AddBranch(&HLT_IsoMu20_eta2p1_TightChargedIsoPFTau27_eta2p1_CrossL1            ,"HLT_IsoMu20_eta2p1_TightChargedIsoPFTau27_eta2p1_CrossL1");
    AddBranch(&HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_1pr       ,"HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_1pr");
    AddBranch(&HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1           ,"HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1");
    AddBranch(&HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET100                 ,"HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET100");
    AddBranch(&HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET90                  ,"HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET90");
    AddBranch(&HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr                        ,"HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr");
    AddBranch(&HLT_PFHT500_PFMET100_PFMHT100_IDTight     ,"HLT_PFHT500_PFMET100_PFMHT100_IDTight");
    AddBranch(&HLT_PFHT700_PFMET85_PFMHT85_IDTight       ,"HLT_PFHT700_PFMET85_PFMHT85_IDTight");
    AddBranch(&HLT_PFHT800_PFMET75_PFMHT75_IDTight       ,"HLT_PFHT800_PFMET75_PFMHT75_IDTight");
    AddBranch(&HLT_CaloJet500_NoJetID                    ,"HLT_CaloJet500_NoJetID");
    AddBranch(&HLT_PFHT1050                              ,"HLT_PFHT1050");
    AddBranch(&HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight   ,"HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight");
    AddBranch(&HLT_PFJet500               ,"HLT_PFJet500");
    AddBranch(&HLT_AK8PFJet500      ,"HLT_AK8PFJet500");
    AddBranch(&HLT_Ele27_eta2p1_WPTight_Gsf             ,"HLT_Ele27_eta2p1_WPTight_Gsf");
    AddBranch(&HLT_Ele27_WPTight_Gsf		      ,"HLT_Ele27_WPTight_Gsf");
    AddBranch(&HLT_Ele25_eta2p1_WPTight_Gsf	      ,"HLT_Ele25_eta2p1_WPTight_Gsf");
    AddBranch(&HLT_Ele115_CaloIdVT_GsfTrkIdT	      ,"HLT_Ele115_CaloIdVT_GsfTrkIdT");
    AddBranch(&HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf    ,"HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf");
    AddBranch(&HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW    ,"HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW");
    AddBranch(&HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ ,"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ");
    AddBranch(&HLT_IsoMu22			,"HLT_IsoMu22");
    AddBranch(&HLT_IsoTkMu22		    ,"HLT_IsoTkMu22");
    AddBranch(&HLT_IsoMu24			,"HLT_IsoMu24");
    AddBranch(&HLT_IsoTkMu24		          ,"HLT_IsoTkMu24");
    AddBranch(&HLT_IsoMu22_eta2p1			  ,"HLT_IsoMu22_eta2p1");
    AddBranch(&HLT_IsoTkMu22_eta2p1		  ,"HLT_IsoTkMu22_eta2p1");
    AddBranch(&HLT_DoubleMu33NoFiltersNoVtx	      ,"HLT_DoubleMu33NoFiltersNoVtx");
    AddBranch(&HLT_DoubleMu23NoFiltersNoVtxDisplaced    ,"HLT_DoubleMu23NoFiltersNoVtxDisplaced");
    AddBranch(&HLT_Mu30_TkMu11			      ,"HLT_Mu30_TkMu11");
    AddBranch(&HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ      ,"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ");
  }
  if(debug_)    std::cout<<"set branches"<<std::endl;
}

void TriggerSelector::Clear(){
  HLT_DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg = -9999;
  HLT_DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg = -9999;
  HLT_DoubleTightChargedIsoPFTau40_Trk1_eta2p1_Reg = -9999;
  HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1 = -9999;
  HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTau30_eta2p1_CrossL1 = -9999;
  HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTau30_eta2p1_CrossL1 = -9999;
  HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1 = -9999;
  HLT_IsoMu20_eta2p1_MediumChargedIsoPFTau27_eta2p1_CrossL1 = -9999;
  HLT_IsoMu20_eta2p1_TightChargedIsoPFTau27_eta2p1_CrossL1 = -9999;
  HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_1pr = -9999;
  HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1 = -9999;
  HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET100 = -9999;
  HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET90 = -9999;
  HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr = -9999;
  HLT_PFHT500_PFMET100_PFMHT100_IDTight = -9999;
  HLT_PFHT700_PFMET85_PFMHT85_IDTight = -9999;
  HLT_PFHT800_PFMET75_PFMHT75_IDTight = -9999;
  HLT_CaloJet500_NoJetID = -9999;
  HLT_PFHT1050 = -9999;
  HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight = -9999;
  HLT_PFMET120_PFMHT120_IDTight = -9999;
  HLT_PFMETNoMu120_PFMHTNoMu120_IDTight = -9999;
  HLT_PFJet500 = -9999;
  HLT_AK8PFJet500 = -9999;
  HLT_Ele27_eta2p1_WPTight_Gsf = -9999;
  HLT_Ele27_WPTight_Gsf = -9999;
  HLT_Ele25_eta2p1_WPTight_Gsf = -9999;
  HLT_Ele115_CaloIdVT_GsfTrkIdT = -9999;
  HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf = -9999;
  HLT_DoubleEle33_CaloIdL_MW = -9999;
  HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW = -9999;
  HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ = -9999;
  HLT_IsoMu22 = -9999;
  HLT_IsoTkMu22 = -9999;
  HLT_IsoMu24 = -9999;
  HLT_IsoTkMu24 = -9999;
  HLT_IsoMu22_eta2p1 = -9999;
  HLT_IsoTkMu22_eta2p1 = -9999;
  HLT_Mu50 = -9999;
  HLT_TkMu50 = -9999;
  HLT_DoubleMu33NoFiltersNoVtx = -9999;
  HLT_DoubleMu23NoFiltersNoVtxDisplaced = -9999;
  HLT_Mu30_TkMu11 = -9999;
  HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ = -9999;

}

void TriggerSelector::startTrigger(edm::EventSetup const& iSetup, edm::Run const & iRun){
  bool changed(true);
  hltConfig_.init(iRun,iSetup,"HLT",changed);
}
