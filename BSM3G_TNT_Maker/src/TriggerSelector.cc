#include "BSMFramework/BSM3G_TNT_Maker/interface/TriggerSelector.h"
#include <string>       // std::string
#include <iostream>     // std::cout
#include <sstream>      // std::stringstream, std::stringbuf
TriggerSelector::TriggerSelector(std::string name, TTree* tree, bool debug, const pset& iConfig):baseTree(name,tree,debug){
  triggerBits_       = iConfig.getParameter<edm::InputTag>("bits");
  _maxtriggerversion = iConfig.getParameter<double>("maxtriggerversion");
  _is_data = iConfig.getParameter<bool>("is_data");
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
    char buffer[10]; sprintf(buffer,"%g",tv);
    //Common
    //Electron
    uint HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v(trigNames.triggerIndex(("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"+string(buffer)).c_str()));
    if(HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v<triggerBits->size()) HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ = triggerBits->accept(HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v);
    //Muon
    uint HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v(trigNames.triggerIndex(("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v"+string(buffer)).c_str()));
    if(HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v<triggerBits->size()) HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ = triggerBits->accept(HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v);
    uint HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v(trigNames.triggerIndex(("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v"+string(buffer)).c_str()));
    if(HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v<triggerBits->size()) HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ = triggerBits->accept(HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v);
    uint HLT_IsoMu20_v(trigNames.triggerIndex(("HLT_IsoMu20_v"+string(buffer)).c_str()));
    if(HLT_IsoMu20_v<triggerBits->size()) HLT_IsoMu20 = triggerBits->accept(HLT_IsoMu20_v);
    uint HLT_IsoMu22_v(trigNames.triggerIndex(("HLT_IsoMu22_v"+string(buffer)).c_str()));
    if(HLT_IsoMu22_v<triggerBits->size()) HLT_IsoMu22 = triggerBits->accept(HLT_IsoMu22_v);
    uint HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v(trigNames.triggerIndex(("HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v"+string(buffer)).c_str()));
    if(HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v<triggerBits->size()) HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL = triggerBits->accept(HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v);
    uint HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v(trigNames.triggerIndex(("HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v"+string(buffer)).c_str()));
    if(HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v<triggerBits->size()) HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL = triggerBits->accept(HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v);

    uint HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v(trigNames.triggerIndex(("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v"+string(buffer)).c_str()));
    if(HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v<triggerBits->size()) HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL = triggerBits->accept(HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v);
    uint HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v(trigNames.triggerIndex(("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v"+string(buffer)).c_str()));
    if(HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v<triggerBits->size()) HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL = triggerBits->accept(HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v);


    //TTHbb
    //Electron
    uint HLT_Ele115_CaloIdVT_GsfTrkIdT_v(trigNames.triggerIndex(("HLT_Ele115_CaloIdVT_GsfTrkIdT_v"+string(buffer)).c_str()));
    if(HLT_Ele115_CaloIdVT_GsfTrkIdT_v<triggerBits->size()) HLT_Ele115_CaloIdVT_GsfTrkIdT = triggerBits->accept(HLT_Ele115_CaloIdVT_GsfTrkIdT_v);
    uint HLT_Ele27_eta2p1_WP75_Gsf_v(trigNames.triggerIndex(("HLT_Ele27_eta2p1_WP75_Gsf_v"+string(buffer)).c_str()));
    if(HLT_Ele27_eta2p1_WP75_Gsf_v<triggerBits->size()) HLT_Ele27_eta2p1_WP75_Gsf = triggerBits->accept(HLT_Ele27_eta2p1_WP75_Gsf_v);
    uint HLT_Ele27_WP85_Gsf_v(trigNames.triggerIndex(("HLT_Ele27_WP85_Gsf_v"+string(buffer)).c_str()));
    if(HLT_Ele27_WP85_Gsf_v<triggerBits->size()) HLT_Ele27_WP85_Gsf = triggerBits->accept(HLT_Ele27_WP85_Gsf_v);
    uint HLT_Ele27_eta2p1_WPLoose_Gsf_v(trigNames.triggerIndex(("HLT_Ele27_eta2p1_WPLoose_Gsf_v"+string(buffer)).c_str()));
    if(HLT_Ele27_eta2p1_WPLoose_Gsf_v<triggerBits->size()) HLT_Ele27_eta2p1_WPLoose_Gsf = triggerBits->accept(HLT_Ele27_eta2p1_WPLoose_Gsf_v);
    uint HLT_Ele27_eta2p1_WPTight_Gsf_v(trigNames.triggerIndex(("HLT_Ele27_eta2p1_WPTight_Gsf_v"+string(buffer)).c_str()));
    if(HLT_Ele27_eta2p1_WPTight_Gsf_v<triggerBits->size()) HLT_Ele27_eta2p1_WPTight_Gsf = triggerBits->accept(HLT_Ele27_eta2p1_WPTight_Gsf_v);
    //Muon
    uint HLT_Mu45_eta2p1_v(trigNames.triggerIndex(("HLT_Mu45_eta2p1_v"+string(buffer)).c_str()));
    if(HLT_Mu45_eta2p1_v<triggerBits->size()) HLT_Mu45_eta2p1 = triggerBits->accept(HLT_Mu45_eta2p1_v);
    uint HLT_Mu50_v(trigNames.triggerIndex(("HLT_Mu50_v"+string(buffer)).c_str()));
    if(HLT_Mu50_v<triggerBits->size()) HLT_Mu50 = triggerBits->accept(HLT_Mu50_v);
    uint HLT_TkMu50_v(trigNames.triggerIndex(("HLT_TkMu50_v"+string(buffer)).c_str()));
    if(HLT_TkMu50_v<triggerBits->size()) HLT_TkMu50 = triggerBits->accept(HLT_TkMu50_v);
    uint HLT_IsoMu17_eta2p1_v(trigNames.triggerIndex(("HLT_IsoMu17_eta2p1_v"+string(buffer)).c_str()));
    if(HLT_IsoMu17_eta2p1_v<triggerBits->size()) HLT_IsoMu17_eta2p1 = triggerBits->accept(HLT_IsoMu17_eta2p1_v);
    uint HLT_IsoMu24_eta2p1_v(trigNames.triggerIndex(("HLT_IsoMu24_eta2p1_v"+string(buffer)).c_str()));
    if(HLT_IsoMu24_eta2p1_v<triggerBits->size()) HLT_IsoMu24_eta2p1 = triggerBits->accept(HLT_IsoMu24_eta2p1_v);
    uint HLT_IsoMu18_v(trigNames.triggerIndex(("HLT_IsoMu18_v"+string(buffer)).c_str()));
    if(HLT_IsoMu18_v<triggerBits->size()) HLT_IsoMu18 = triggerBits->accept(HLT_IsoMu18_v);
    //Cross
    uint HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v(trigNames.triggerIndex(("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"+string(buffer)).c_str()));
    if(HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v<triggerBits->size()) HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ = triggerBits->accept(HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v);
    uint HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v(trigNames.triggerIndex(("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v"+string(buffer)).c_str()));
    if(HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v<triggerBits->size()) HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL = triggerBits->accept(HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v);
    uint HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v(trigNames.triggerIndex(("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v"+string(buffer)).c_str()));
    if(HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v<triggerBits->size()) HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL = triggerBits->accept(HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v);
    uint HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v(trigNames.triggerIndex(("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"+string(buffer)).c_str()));
    if(HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v<triggerBits->size()) HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ = triggerBits->accept(HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v);
    uint HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v(trigNames.triggerIndex(("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v"+string(buffer)).c_str()));
    if(HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v<triggerBits->size()) HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ = triggerBits->accept(HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v);
    //TTHLep
    //Electron
    uint HLT_Ele27_WPTight_Gsf_v(trigNames.triggerIndex(("HLT_Ele27_WPTight_Gsf_v"+string(buffer)).c_str()));
    if(HLT_Ele27_WPTight_Gsf_v<triggerBits->size()) HLT_Ele27_WPTight_Gsf = triggerBits->accept(HLT_Ele27_WPTight_Gsf_v);
    uint HLT_Ele25_eta2p1_WPTight_Gsf_v(trigNames.triggerIndex(("HLT_Ele25_eta2p1_WPTight_Gsf_v"+string(buffer)).c_str()));
    if(HLT_Ele25_eta2p1_WPTight_Gsf_v<triggerBits->size()) HLT_Ele25_eta2p1_WPTight_Gsf = triggerBits->accept(HLT_Ele25_eta2p1_WPTight_Gsf_v);
    uint HLT_Ele23_WPLoose_Gsf_v(trigNames.triggerIndex(("HLT_Ele23_WPLoose_Gsf_v"+string(buffer)).c_str())); //Data
    if(HLT_Ele23_WPLoose_Gsf_v<triggerBits->size()) HLT_Ele23_WPLoose_Gsf = triggerBits->accept(HLT_Ele23_WPLoose_Gsf_v);
    uint HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v(trigNames.triggerIndex(("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v"+string(buffer)).c_str())); //MC
    if(HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v<triggerBits->size()) HLT_Ele23_CaloIdL_TrackIdL_IsoVL = triggerBits->accept(HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v);
    //Muon
    uint HLT_IsoMu22_eta2p1_v(trigNames.triggerIndex(("HLT_IsoMu22_eta2p1_v"+string(buffer)).c_str()));
    if(HLT_IsoMu22_eta2p1_v<triggerBits->size()) HLT_IsoMu22_eta2p1 = triggerBits->accept(HLT_IsoMu22_eta2p1_v);
    uint HLT_IsoTkMu22_eta2p1_v(trigNames.triggerIndex(("HLT_IsoTkMu22_eta2p1_v"+string(buffer)).c_str()));
    if(HLT_IsoTkMu22_eta2p1_v<triggerBits->size()) HLT_IsoTkMu22_eta2p1 = triggerBits->accept(HLT_IsoTkMu22_eta2p1_v);
    uint HLT_IsoTkMu20_v(trigNames.triggerIndex(("HLT_IsoTkMu20_v"+string(buffer)).c_str()));
    if(HLT_IsoTkMu20_v<triggerBits->size()) HLT_IsoTkMu20 = triggerBits->accept(HLT_IsoTkMu20_v);
    uint HLT_IsoTkMu22_v(trigNames.triggerIndex(("HLT_IsoTkMu22_v"+string(buffer)).c_str()));
    if(HLT_IsoTkMu22_v<triggerBits->size()) HLT_IsoTkMu22 = triggerBits->accept(HLT_IsoTkMu22_v);
    uint HLT_IsoMu24_v(trigNames.triggerIndex(("HLT_IsoMu24_v"+string(buffer)).c_str()));
    if(HLT_IsoMu24_v<triggerBits->size()) HLT_IsoMu24 = triggerBits->accept(HLT_IsoMu24_v);
    uint HLT_IsoTkMu24_v(trigNames.triggerIndex(("HLT_IsoTkMu24_v"+string(buffer)).c_str()));
    if(HLT_IsoTkMu24_v<triggerBits->size()) HLT_IsoTkMu24 = triggerBits->accept(HLT_IsoTkMu24_v);
    //CrossEle-Mu
    uint HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v(trigNames.triggerIndex(("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v"+string(buffer)).c_str()));
    if(HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v<triggerBits->size()) HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL = triggerBits->accept(HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v);
    uint HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_v(trigNames.triggerIndex(("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_v"+string(buffer)).c_str()));
    if(HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_v<triggerBits->size()) HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ = triggerBits->accept(HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_v);
    uint HLT_DiMu9_Ele9_CaloIdL_TrackIdL_v(trigNames.triggerIndex(("HLT_DiMu9_Ele9_CaloIdL_TrackIdL_v"+string(buffer)).c_str()));
    if(HLT_DiMu9_Ele9_CaloIdL_TrackIdL_v<triggerBits->size()) HLT_DiMu9_Ele9_CaloIdL_TrackIdL = triggerBits->accept(HLT_DiMu9_Ele9_CaloIdL_TrackIdL_v);
    uint HLT_Mu8_DiEle12_CaloIdL_TrackIdL_v(trigNames.triggerIndex(("HLT_Mu8_DiEle12_CaloIdL_TrackIdL_v"+string(buffer)).c_str()));
    if(HLT_Mu8_DiEle12_CaloIdL_TrackIdL_v<triggerBits->size()) HLT_Mu8_DiEle12_CaloIdL_TrackIdL = triggerBits->accept(HLT_Mu8_DiEle12_CaloIdL_TrackIdL_v);
    uint HLT_TripleMu_12_10_5_v(trigNames.triggerIndex(("HLT_TripleMu_12_10_5_v"+string(buffer)).c_str()));
    if(HLT_TripleMu_12_10_5_v<triggerBits->size()) HLT_TripleMu_12_10_5 = triggerBits->accept(HLT_TripleMu_12_10_5_v);
    uint HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v(trigNames.triggerIndex(("HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v"+string(buffer)).c_str()));
    if(HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v<triggerBits->size()) HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL = triggerBits->accept(HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v);
    //Other
    uint HLT_DoubleEle33_CaloIdL_MW_v(trigNames.triggerIndex(("HLT_DoubleEle33_CaloIdL_MW_v"+string(buffer)).c_str()));
    if(HLT_DoubleEle33_CaloIdL_MW_v<triggerBits->size()) HLT_DoubleEle33_CaloIdL_MW = triggerBits->accept(HLT_DoubleEle33_CaloIdL_MW_v);
    uint HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v(trigNames.triggerIndex(("HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v"+string(buffer)).c_str()));
    if(HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v<triggerBits->size()) HLT_DoubleEle33_CaloIdL_GsfTrkIdVL = triggerBits->accept(HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v);

    //Analysis
    //uint TTHbb_SL_v(trigNames.triggerIndex(("TTHbb_SL_v"+string(buffer)).c_str()));
    //if(TTHbb_SL_v<triggerBits->size()) TTHbb_SL = triggerBits->accept(TTHbb_SL_v);
    //uint TTHbb_DL_v(trigNames.triggerIndex(("TTHbb_DL_v"+string(buffer)).c_str()));
    //if(TTHbb_DL_v<triggerBits->size()) TTHbb_DL = triggerBits->accept(TTHbb_DL_v);
    //uint TTHLep_2Mu_v(trigNames.triggerIndex(("TTHLep_2Mu_v"+string(buffer)).c_str()));
    //if(TTHLep_2Mu_v<triggerBits->size()) TTHLep_2Mu = triggerBits->accept(TTHLep_2Mu_v);
    //uint TTHLep_2Ele_v(trigNames.triggerIndex(("TTHLep_2Ele_v"+string(buffer)).c_str()));
    //if(TTHLep_2Ele_v<triggerBits->size()) TTHLep_2Ele = triggerBits->accept(TTHLep_2Ele_v);
    //uint TTHLep_MuEle_v(trigNames.triggerIndex(("TTHLep_MuEle_v"+string(buffer)).c_str()));
    //if(TTHLep_MuEle_v<triggerBits->size()) TTHLep_MuEle = triggerBits->accept(TTHLep_MuEle_v);
    //uint TTHLep_3L4L_v(trigNames.triggerIndex(("TTHLep_3L4L_v"+string(buffer)).c_str()));
    //if(TTHLep_3L4L_v<triggerBits->size()) TTHLep_3L4L = triggerBits->accept(TTHLep_3L4L_v);
    //if(HLT_IsoMu20_v<triggerBits->size()) cout<<"Trigger name is "<<HLT_IsoMu20_v<<" "<<path<<" "<<trigNames.triggerName(HLT_IsoMu20_v)<<" "<<triggerBits->accept(HLT_IsoMu20_v)<<endl;
  }
}

void TriggerSelector::SetBranches(){
  if(debug_)    std::cout<<"setting branches: calling AddBranch of baseTree"<<std::endl;
  //AddBranch(&Trigger_decision,   "Trigger_decision");
  //AddBranch(&Trigger_names,   "Trigger_names");
  //Common
  //Electron
  AddBranch(&HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ       ,"HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ");
  //Muon
  AddBranch(&HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ             ,"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ");
  AddBranch(&HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ           ,"HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ");
  AddBranch(&HLT_IsoMu20                                     ,"HLT_IsoMu20");
  AddBranch(&HLT_IsoMu22                                     ,"HLT_IsoMu22");
  AddBranch(&HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL ,"HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL");
  AddBranch(&HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL  ,"HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL");
  //TTHbb
  //Electron
  AddBranch(&HLT_Ele115_CaloIdVT_GsfTrkIdT                   ,"HLT_Ele115_CaloIdVT_GsfTrkIdT");
  AddBranch(&HLT_Ele27_eta2p1_WP75_Gsf                       ,"HLT_Ele27_eta2p1_WP75_Gsf");
  AddBranch(&HLT_Ele27_WP85_Gsf                              ,"HLT_Ele27_WP85_Gsf");
  AddBranch(&HLT_Ele27_eta2p1_WPLoose_Gsf                    ,"HLT_Ele27_eta2p1_WPLoose_Gsf");
  AddBranch(&HLT_Ele27_eta2p1_WPTight_Gsf                    ,"HLT_Ele27_eta2p1_WPTight_Gsf");
  //Muon
  AddBranch(&HLT_Mu45_eta2p1                                 ,"HLT_Mu45_eta2p1");
  AddBranch(&HLT_Mu50                                        ,"HLT_Mu50");
  AddBranch(&HLT_TkMu50                                      ,"HLT_TkMu50");
  AddBranch(&HLT_IsoMu17_eta2p1                              ,"HLT_IsoMu17_eta2p1");
  AddBranch(&HLT_IsoMu24_eta2p1                              ,"HLT_IsoMu24_eta2p1");
  AddBranch(&HLT_IsoMu18                                     ,"HLT_IsoMu18");
  //Cross
  AddBranch(&HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ       ,"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ");
  AddBranch(&HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL ,"HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL");
  AddBranch(&HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL  ,"HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL");
  //TTHLep
  //Electron
  AddBranch(&HLT_Ele23_WPLoose_Gsf                           ,"HLT_Ele23_WPLoose_Gsf"); //Data
  AddBranch(&HLT_Ele23_CaloIdL_TrackIdL_IsoVL                ,"HLT_Ele23_CaloIdL_TrackIdL_IsoVL"); //MC
  AddBranch(&HLT_Ele27_WPTight_Gsf,"HLT_Ele27_WPTight_Gsf");
  AddBranch(&HLT_Ele25_eta2p1_WPTight_Gsf,"HLT_Ele25_eta2p1_WPTight_Gsf");
  //Muon
  AddBranch(&HLT_IsoMu22_eta2p1,"HLT_IsoMu22_eta2p1");
  AddBranch(&HLT_IsoTkMu22_eta2p1,"HLT_IsoTkMu22_eta2p1");
  AddBranch(&HLT_IsoTkMu20                                   ,"HLT_IsoTkMu20");
  AddBranch(&HLT_IsoTkMu22                                   ,"HLT_IsoTkMu22");
  AddBranch(&HLT_IsoMu24,"HLT_IsoMu24");
  AddBranch(&HLT_IsoTkMu24,"HLT_IsoTkMu24");
  AddBranch(&HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ,"HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ");
  AddBranch(&HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ,"HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ");
  AddBranch(&HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL,"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL");
  AddBranch(&HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL,"HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL");
  AddBranch(&HLT_DoubleEle33_CaloIdL_GsfTrkIdVL,"HLT_DoubleEle33_CaloIdL_GsfTrkIdVL");
  //CrossEle-Mu
  AddBranch(&HLT_DiMu9_Ele9_CaloIdL_TrackIdL                 ,"HLT_DiMu9_Ele9_CaloIdL_TrackIdL");
  AddBranch(&HLT_Mu8_DiEle12_CaloIdL_TrackIdL                ,"HLT_Mu8_DiEle12_CaloIdL_TrackIdL");
  AddBranch(&HLT_TripleMu_12_10_5                            ,"HLT_TripleMu_12_10_5");
  AddBranch(&HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL           ,"HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL");
  AddBranch(&HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL,"HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL");
  AddBranch(&HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ,"HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ");
  //Other
  AddBranch(&HLT_DoubleEle33_CaloIdL_MW              ,"HLT_DoubleEle33_CaloIdL_MW");
  //Analysis
  AddBranch(&TTHbb_SL     ,"TTHbb_SL");
  AddBranch(&TTHbb_DL     ,"TTHbb_DL");
  AddBranch(&TTHLep_2Mu   ,"TTHLep_2Mu");
  AddBranch(&TTHLep_2Ele  ,"TTHLep_2Ele");
  AddBranch(&TTHLep_MuEle ,"TTHLep_MuEle");
  AddBranch(&TTHLep_3L4L  ,"TTHLep_3L4L");
  if(debug_)    std::cout<<"set branches"<<std::endl;
}

void TriggerSelector::Clear(){
  //Trigger_decision.clear();
  //Trigger_names.clear();
  //Common
  //Electron
  HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ = -9999;
  //Muon
  HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ             = -9999;
  HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ           = -9999;
  HLT_IsoMu20                                     = -9999;
  HLT_IsoMu22                                     = -9999;
  HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL = -9999;
  HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL  = -9999;
  //TTHbb
  //Electron
  HLT_Ele115_CaloIdVT_GsfTrkIdT = -9999;
  HLT_Ele27_eta2p1_WP75_Gsf     = -9999;
  HLT_Ele27_WP85_Gsf            = -9999;
  HLT_Ele27_eta2p1_WPLoose_Gsf  = -9999;
  HLT_Ele27_eta2p1_WPTight_Gsf  = -9999;
  //Muon
  HLT_Mu45_eta2p1    = -9999;
  HLT_Mu50           = -9999;
  HLT_TkMu50         = -9999;
  HLT_IsoMu17_eta2p1 = -9999;
  HLT_IsoMu24_eta2p1 = -9999;
  HLT_IsoMu18        = -9999;
  //Cross
  HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ       = -9999;
  HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL = -9999;
  HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL  = -9999;
  //TTHLep
  //Electron
  HLT_Ele23_WPLoose_Gsf            = -9999; //Data
  HLT_Ele23_CaloIdL_TrackIdL_IsoVL = -9999; //MC
  HLT_Ele27_WPTight_Gsf= -9999;
  HLT_Ele25_eta2p1_WPTight_Gsf= -9999;
  //Muon
  HLT_IsoMu22_eta2p1= -9999;
  HLT_IsoTkMu22_eta2p1= -9999;
  HLT_IsoTkMu20 = -9999;
  HLT_IsoTkMu22 = -9999;
  HLT_IsoMu24 = -9999;
  HLT_IsoTkMu24 = -9999;
  HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ = -9999;
  HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ = -9999;
  HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL = -9999;
  HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL = -9999;
  HLT_DoubleEle33_CaloIdL_GsfTrkIdVL = -9999;

  //CrossEle-Mu
  HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL= -9999;
  HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ= -9999;
  HLT_DiMu9_Ele9_CaloIdL_TrackIdL       = -9999;
  HLT_Mu8_DiEle12_CaloIdL_TrackIdL      = -9999;
  HLT_TripleMu_12_10_5                  = -9999;
  HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL = -9999;
  //Other
  HLT_DoubleEle33_CaloIdL_MW    = -9999;
  //Analysis
  TTHbb_SL     = -9999;
  TTHbb_DL     = -9999;
  TTHLep_2Mu   = -9999;
  TTHLep_2Ele  = -9999;
  TTHLep_MuEle = -9999;
  TTHLep_3L4L  = -9999;
}

void TriggerSelector::startTrigger(edm::EventSetup const& iSetup, edm::Run const & iRun){
  bool changed(true);
  if(_is_data) hltConfig_.init(iRun,iSetup,"HLT",changed);
  else         hltConfig_.init(iRun,iSetup,"HLT",changed);
}
