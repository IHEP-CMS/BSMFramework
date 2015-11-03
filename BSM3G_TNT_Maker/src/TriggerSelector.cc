#include "BSMFramework/BSM3G_TNT_Maker/interface/TriggerSelector.h"

TriggerSelector::TriggerSelector(std::string name, TTree* tree, bool debug, const pset& iConfig):baseTree(name,tree,debug){
  triggerResultsTag_  = iConfig.getParameter<edm::InputTag>("triggerResults");
  SetBranches();
}

TriggerSelector::~TriggerSelector(){
  delete tree_;
}

void TriggerSelector::Fill(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  
  if(debug_)    std::cout<<"getting met info"<<std::endl;
   
  Clear();

  //Trigget paths  
  edm::Handle<edm::TriggerResults> trigResults;
  edm::InputTag trigResultsTag("TriggerResults","","HLT");  
  iEvent.getByLabel(trigResultsTag,trigResults);
  const edm::TriggerNames& trigNames = iEvent.triggerNames(*trigResults);

  /*
  std::cout<<std::endl; std::cout<<std::endl;
  for (unsigned int i=0; i<trigNames.size(); i++) {
    std::cout<<trigNames.triggerName(i) << std::endl;
  }
  std::cout<<std::endl; std::cout<<std::endl;
  */

  unsigned int HLT_Ele105_CaloIdVT_GsfTrkIdT_v1(trigNames.triggerIndex("HLT_Ele105_CaloIdVT_GsfTrkIdT_v1"));
  unsigned int HLT_Ele105_CaloIdVT_GsfTrkIdT_v2(trigNames.triggerIndex("HLT_Ele105_CaloIdVT_GsfTrkIdT_v2"));
  unsigned int HLT_Ele105_CaloIdVT_GsfTrkIdT_v3(trigNames.triggerIndex("HLT_Ele105_CaloIdVT_GsfTrkIdT_v3"));
  unsigned int HLT_Ele105_CaloIdVT_GsfTrkIdT_v4(trigNames.triggerIndex("HLT_Ele105_CaloIdVT_GsfTrkIdT_v4"));
  unsigned int HLT_Ele105_CaloIdVT_GsfTrkIdT_v5(trigNames.triggerIndex("HLT_Ele105_CaloIdVT_GsfTrkIdT_v5"));
  if(HLT_Ele105_CaloIdVT_GsfTrkIdT_v1 < trigResults->size()) HLT_Ele105_CaloIdVT_GsfTrkIdT = trigResults->accept(HLT_Ele105_CaloIdVT_GsfTrkIdT_v1);
  if(HLT_Ele105_CaloIdVT_GsfTrkIdT_v2 < trigResults->size()) HLT_Ele105_CaloIdVT_GsfTrkIdT = trigResults->accept(HLT_Ele105_CaloIdVT_GsfTrkIdT_v2);
  if(HLT_Ele105_CaloIdVT_GsfTrkIdT_v3 < trigResults->size()) HLT_Ele105_CaloIdVT_GsfTrkIdT = trigResults->accept(HLT_Ele105_CaloIdVT_GsfTrkIdT_v3);
  if(HLT_Ele105_CaloIdVT_GsfTrkIdT_v4 < trigResults->size()) HLT_Ele105_CaloIdVT_GsfTrkIdT = trigResults->accept(HLT_Ele105_CaloIdVT_GsfTrkIdT_v4);
  if(HLT_Ele105_CaloIdVT_GsfTrkIdT_v5 < trigResults->size()) HLT_Ele105_CaloIdVT_GsfTrkIdT = trigResults->accept(HLT_Ele105_CaloIdVT_GsfTrkIdT_v5);
  unsigned int HLT_Ele27_eta2p1_WP75_Gsf_v1(trigNames.triggerIndex("HLT_Ele27_eta2p1_WP75_Gsf_v1"));
  unsigned int HLT_Ele27_eta2p1_WP75_Gsf_v2(trigNames.triggerIndex("HLT_Ele27_eta2p1_WP75_Gsf_v2"));
  unsigned int HLT_Ele27_eta2p1_WP75_Gsf_v3(trigNames.triggerIndex("HLT_Ele27_eta2p1_WP75_Gsf_v3"));
  unsigned int HLT_Ele27_eta2p1_WP75_Gsf_v4(trigNames.triggerIndex("HLT_Ele27_eta2p1_WP75_Gsf_v4"));
  unsigned int HLT_Ele27_eta2p1_WP75_Gsf_v5(trigNames.triggerIndex("HLT_Ele27_eta2p1_WP75_Gsf_v5"));
  if(HLT_Ele27_eta2p1_WP75_Gsf_v1 < trigResults->size()) HLT_Ele27_eta2p1_WP75_Gsf = trigResults->accept(HLT_Ele27_eta2p1_WP75_Gsf_v1);
  if(HLT_Ele27_eta2p1_WP75_Gsf_v2 < trigResults->size()) HLT_Ele27_eta2p1_WP75_Gsf = trigResults->accept(HLT_Ele27_eta2p1_WP75_Gsf_v2);
  if(HLT_Ele27_eta2p1_WP75_Gsf_v3 < trigResults->size()) HLT_Ele27_eta2p1_WP75_Gsf = trigResults->accept(HLT_Ele27_eta2p1_WP75_Gsf_v3);
  if(HLT_Ele27_eta2p1_WP75_Gsf_v4 < trigResults->size()) HLT_Ele27_eta2p1_WP75_Gsf = trigResults->accept(HLT_Ele27_eta2p1_WP75_Gsf_v4);
  if(HLT_Ele27_eta2p1_WP75_Gsf_v5 < trigResults->size()) HLT_Ele27_eta2p1_WP75_Gsf = trigResults->accept(HLT_Ele27_eta2p1_WP75_Gsf_v5);
  unsigned int HLT_Ele27_WP85_Gsf_v1(trigNames.triggerIndex("HLT_Ele27_WP85_Gsf_v1"));
  unsigned int HLT_Ele27_WP85_Gsf_v2(trigNames.triggerIndex("HLT_Ele27_WP85_Gsf_v2"));
  unsigned int HLT_Ele27_WP85_Gsf_v3(trigNames.triggerIndex("HLT_Ele27_WP85_Gsf_v3"));
  unsigned int HLT_Ele27_WP85_Gsf_v4(trigNames.triggerIndex("HLT_Ele27_WP85_Gsf_v4"));
  unsigned int HLT_Ele27_WP85_Gsf_v5(trigNames.triggerIndex("HLT_Ele27_WP85_Gsf_v5"));
  if(HLT_Ele27_WP85_Gsf_v1 < trigResults->size()) HLT_Ele27_WP85_Gsf = trigResults->accept(HLT_Ele27_WP85_Gsf_v1);
  if(HLT_Ele27_WP85_Gsf_v2 < trigResults->size()) HLT_Ele27_WP85_Gsf = trigResults->accept(HLT_Ele27_WP85_Gsf_v2);
  if(HLT_Ele27_WP85_Gsf_v3 < trigResults->size()) HLT_Ele27_WP85_Gsf = trigResults->accept(HLT_Ele27_WP85_Gsf_v3);
  if(HLT_Ele27_WP85_Gsf_v4 < trigResults->size()) HLT_Ele27_WP85_Gsf = trigResults->accept(HLT_Ele27_WP85_Gsf_v4);
  if(HLT_Ele27_WP85_Gsf_v5 < trigResults->size()) HLT_Ele27_WP85_Gsf = trigResults->accept(HLT_Ele27_WP85_Gsf_v5);
  unsigned int HLT_Ele27_eta2p1_WPLoose_Gsf_v1(trigNames.triggerIndex("HLT_Ele27_eta2p1_WPLoose_Gsf_v1"));
  unsigned int HLT_Ele27_eta2p1_WPLoose_Gsf_v2(trigNames.triggerIndex("HLT_Ele27_eta2p1_WPLoose_Gsf_v2"));
  unsigned int HLT_Ele27_eta2p1_WPLoose_Gsf_v3(trigNames.triggerIndex("HLT_Ele27_eta2p1_WPLoose_Gsf_v3"));
  unsigned int HLT_Ele27_eta2p1_WPLoose_Gsf_v4(trigNames.triggerIndex("HLT_Ele27_eta2p1_WPLoose_Gsf_v4"));
  unsigned int HLT_Ele27_eta2p1_WPLoose_Gsf_v5(trigNames.triggerIndex("HLT_Ele27_eta2p1_WPLoose_Gsf_v5"));
  if(HLT_Ele27_eta2p1_WPLoose_Gsf_v1 < trigResults->size()) HLT_Ele27_eta2p1_WPLoose_Gsf = trigResults->accept(HLT_Ele27_eta2p1_WPLoose_Gsf_v1);
  if(HLT_Ele27_eta2p1_WPLoose_Gsf_v2 < trigResults->size()) HLT_Ele27_eta2p1_WPLoose_Gsf = trigResults->accept(HLT_Ele27_eta2p1_WPLoose_Gsf_v2);
  if(HLT_Ele27_eta2p1_WPLoose_Gsf_v3 < trigResults->size()) HLT_Ele27_eta2p1_WPLoose_Gsf = trigResults->accept(HLT_Ele27_eta2p1_WPLoose_Gsf_v3);
  if(HLT_Ele27_eta2p1_WPLoose_Gsf_v4 < trigResults->size()) HLT_Ele27_eta2p1_WPLoose_Gsf = trigResults->accept(HLT_Ele27_eta2p1_WPLoose_Gsf_v4);
  if(HLT_Ele27_eta2p1_WPLoose_Gsf_v5 < trigResults->size()) HLT_Ele27_eta2p1_WPLoose_Gsf = trigResults->accept(HLT_Ele27_eta2p1_WPLoose_Gsf_v5);
  unsigned int HLT_Mu50_v1(trigNames.triggerIndex("HLT_Mu50_v1"));
  unsigned int HLT_Mu50_v2(trigNames.triggerIndex("HLT_Mu50_v2"));
  unsigned int HLT_Mu50_v3(trigNames.triggerIndex("HLT_Mu50_v3"));
  unsigned int HLT_Mu50_v4(trigNames.triggerIndex("HLT_Mu50_v4"));
  unsigned int HLT_Mu50_v5(trigNames.triggerIndex("HLT_Mu50_v5"));
  if(HLT_Mu50_v1 < trigResults->size()) HLT_Mu50 = trigResults->accept(HLT_Mu50_v1);
  if(HLT_Mu50_v2 < trigResults->size()) HLT_Mu50 = trigResults->accept(HLT_Mu50_v2);
  if(HLT_Mu50_v3 < trigResults->size()) HLT_Mu50 = trigResults->accept(HLT_Mu50_v3);
  if(HLT_Mu50_v4 < trigResults->size()) HLT_Mu50 = trigResults->accept(HLT_Mu50_v4);
  if(HLT_Mu50_v5 < trigResults->size()) HLT_Mu50 = trigResults->accept(HLT_Mu50_v5);
  unsigned int HLT_IsoMu20_v1(trigNames.triggerIndex("HLT_IsoMu20_v1"));
  unsigned int HLT_IsoMu20_v2(trigNames.triggerIndex("HLT_IsoMu20_v2"));
  unsigned int HLT_IsoMu20_v3(trigNames.triggerIndex("HLT_IsoMu20_v3"));
  unsigned int HLT_IsoMu20_v4(trigNames.triggerIndex("HLT_IsoMu20_v4"));
  unsigned int HLT_IsoMu20_v5(trigNames.triggerIndex("HLT_IsoMu20_v5"));
  if(HLT_IsoMu20_v1 < trigResults->size()) HLT_IsoMu20 = trigResults->accept(HLT_IsoMu20_v1);
  if(HLT_IsoMu20_v2 < trigResults->size()) HLT_IsoMu20 = trigResults->accept(HLT_IsoMu20_v2);
  if(HLT_IsoMu20_v3 < trigResults->size()) HLT_IsoMu20 = trigResults->accept(HLT_IsoMu20_v3);
  if(HLT_IsoMu20_v4 < trigResults->size()) HLT_IsoMu20 = trigResults->accept(HLT_IsoMu20_v4);
  if(HLT_IsoMu20_v5 < trigResults->size()) HLT_IsoMu20 = trigResults->accept(HLT_IsoMu20_v5);
  unsigned int HLT_IsoMu17_eta2p1_v1(trigNames.triggerIndex("HLT_IsoMu17_eta2p1_v1"));
  unsigned int HLT_IsoMu17_eta2p1_v2(trigNames.triggerIndex("HLT_IsoMu17_eta2p1_v2"));
  unsigned int HLT_IsoMu17_eta2p1_v3(trigNames.triggerIndex("HLT_IsoMu17_eta2p1_v3"));
  unsigned int HLT_IsoMu17_eta2p1_v4(trigNames.triggerIndex("HLT_IsoMu17_eta2p1_v4"));
  unsigned int HLT_IsoMu17_eta2p1_v5(trigNames.triggerIndex("HLT_IsoMu17_eta2p1_v5"));
  if(HLT_IsoMu17_eta2p1_v1 < trigResults->size()) HLT_IsoMu17_eta2p1 = trigResults->accept(HLT_IsoMu17_eta2p1_v1);
  if(HLT_IsoMu17_eta2p1_v2 < trigResults->size()) HLT_IsoMu17_eta2p1 = trigResults->accept(HLT_IsoMu17_eta2p1_v2);
  if(HLT_IsoMu17_eta2p1_v3 < trigResults->size()) HLT_IsoMu17_eta2p1 = trigResults->accept(HLT_IsoMu17_eta2p1_v3);
  if(HLT_IsoMu17_eta2p1_v4 < trigResults->size()) HLT_IsoMu17_eta2p1 = trigResults->accept(HLT_IsoMu17_eta2p1_v4);
  if(HLT_IsoMu17_eta2p1_v5 < trigResults->size()) HLT_IsoMu17_eta2p1 = trigResults->accept(HLT_IsoMu17_eta2p1_v5);
  unsigned int HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v1(trigNames.triggerIndex("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v1"));
  unsigned int HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v2(trigNames.triggerIndex("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v2"));
  unsigned int HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3(trigNames.triggerIndex("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3"));
  unsigned int HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v4(trigNames.triggerIndex("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v4"));
  unsigned int HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v5(trigNames.triggerIndex("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v5"));
  if(HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v1 < trigResults->size()) HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ = trigResults->accept(HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v1);
  if(HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v2 < trigResults->size()) HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ = trigResults->accept(HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v2);
  if(HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3 < trigResults->size()) HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ = trigResults->accept(HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3);
  if(HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v4 < trigResults->size()) HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ = trigResults->accept(HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v4);
  if(HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v5 < trigResults->size()) HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ = trigResults->accept(HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v5);
  unsigned int HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v1(trigNames.triggerIndex("HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v1"));
  unsigned int HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v2(trigNames.triggerIndex("HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v2"));
  unsigned int HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v3(trigNames.triggerIndex("HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v3"));
  unsigned int HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v4(trigNames.triggerIndex("HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v4"));
  unsigned int HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v5(trigNames.triggerIndex("HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v5"));
  if(HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v1 < trigResults->size()) HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL = trigResults->accept(HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v1);
  if(HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v2 < trigResults->size()) HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL = trigResults->accept(HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v2);
  if(HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v3 < trigResults->size()) HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL = trigResults->accept(HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v3);
  if(HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v4 < trigResults->size()) HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL = trigResults->accept(HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v4);
  if(HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v5 < trigResults->size()) HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL = trigResults->accept(HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v5);
  unsigned int HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v1(trigNames.triggerIndex("HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v1"));
  unsigned int HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v2(trigNames.triggerIndex("HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v2"));
  unsigned int HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v3(trigNames.triggerIndex("HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v3"));
  unsigned int HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v4(trigNames.triggerIndex("HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v4"));
  unsigned int HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v5(trigNames.triggerIndex("HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v5"));
  if(HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v1 < trigResults->size()) HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL = trigResults->accept(HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v1);
  if(HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v2 < trigResults->size()) HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL = trigResults->accept(HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v2);
  if(HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v3 < trigResults->size()) HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL = trigResults->accept(HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v3);
  if(HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v4 < trigResults->size()) HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL = trigResults->accept(HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v4);
  if(HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v5 < trigResults->size()) HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL = trigResults->accept(HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v5);
  unsigned int HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v1(trigNames.triggerIndex("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v1"));
  unsigned int HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v2(trigNames.triggerIndex("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v2"));
  unsigned int HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v3(trigNames.triggerIndex("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v3"));
  unsigned int HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v4(trigNames.triggerIndex("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v4"));
  unsigned int HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v5(trigNames.triggerIndex("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v5"));
  if(HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v1 < trigResults->size()) HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ = trigResults->accept(HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v1);
  if(HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v2 < trigResults->size()) HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ = trigResults->accept(HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v2);
  if(HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v3 < trigResults->size()) HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ = trigResults->accept(HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v3);
  if(HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v4 < trigResults->size()) HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ = trigResults->accept(HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v4);
  if(HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v5 < trigResults->size()) HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ = trigResults->accept(HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v5);
  unsigned int HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v1(trigNames.triggerIndex("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v1"));
  unsigned int HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v2(trigNames.triggerIndex("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v2"));
  unsigned int HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v3(trigNames.triggerIndex("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v3"));
  unsigned int HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v4(trigNames.triggerIndex("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v4"));
  unsigned int HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v5(trigNames.triggerIndex("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v5"));
  if(HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v1 < trigResults->size()) HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ = trigResults->accept(HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v1);
  if(HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v2 < trigResults->size()) HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ = trigResults->accept(HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v2);
  if(HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v3 < trigResults->size()) HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ = trigResults->accept(HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v3);
  if(HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v4 < trigResults->size()) HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ = trigResults->accept(HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v4);
  if(HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v5 < trigResults->size()) HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ = trigResults->accept(HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v5);
  unsigned int HLT_IsoMu24_eta2p1_v1(trigNames.triggerIndex("HLT_IsoMu24_eta2p1_v1"));
  unsigned int HLT_IsoMu24_eta2p1_v2(trigNames.triggerIndex("HLT_IsoMu24_eta2p1_v2"));
  unsigned int HLT_IsoMu24_eta2p1_v3(trigNames.triggerIndex("HLT_IsoMu24_eta2p1_v3"));
  unsigned int HLT_IsoMu24_eta2p1_v4(trigNames.triggerIndex("HLT_IsoMu24_eta2p1_v4"));
  unsigned int HLT_IsoMu24_eta2p1_v5(trigNames.triggerIndex("HLT_IsoMu24_eta2p1_v5"));
  if(HLT_IsoMu24_eta2p1_v1 < trigResults->size()) HLT_IsoMu24_eta2p1 = trigResults->accept(HLT_IsoMu24_eta2p1_v1);
  if(HLT_IsoMu24_eta2p1_v2 < trigResults->size()) HLT_IsoMu24_eta2p1 = trigResults->accept(HLT_IsoMu24_eta2p1_v2);
  if(HLT_IsoMu24_eta2p1_v3 < trigResults->size()) HLT_IsoMu24_eta2p1 = trigResults->accept(HLT_IsoMu24_eta2p1_v3);
  if(HLT_IsoMu24_eta2p1_v4 < trigResults->size()) HLT_IsoMu24_eta2p1 = trigResults->accept(HLT_IsoMu24_eta2p1_v4);
  if(HLT_IsoMu24_eta2p1_v5 < trigResults->size()) HLT_IsoMu24_eta2p1 = trigResults->accept(HLT_IsoMu24_eta2p1_v5);
  unsigned int HLT_IsoMu18_v1(trigNames.triggerIndex("HLT_IsoMu18_v1"));
  unsigned int HLT_IsoMu18_v2(trigNames.triggerIndex("HLT_IsoMu18_v2"));
  unsigned int HLT_IsoMu18_v3(trigNames.triggerIndex("HLT_IsoMu18_v3"));
  unsigned int HLT_IsoMu18_v4(trigNames.triggerIndex("HLT_IsoMu18_v4"));
  unsigned int HLT_IsoMu18_v5(trigNames.triggerIndex("HLT_IsoMu18_v5"));
  if(HLT_IsoMu18_v1 < trigResults->size()) HLT_IsoMu18 = trigResults->accept(HLT_IsoMu18_v1);
  if(HLT_IsoMu18_v2 < trigResults->size()) HLT_IsoMu18 = trigResults->accept(HLT_IsoMu18_v2);
  if(HLT_IsoMu18_v3 < trigResults->size()) HLT_IsoMu18 = trigResults->accept(HLT_IsoMu18_v3);
  if(HLT_IsoMu18_v4 < trigResults->size()) HLT_IsoMu18 = trigResults->accept(HLT_IsoMu18_v4);
  if(HLT_IsoMu18_v5 < trigResults->size()) HLT_IsoMu18 = trigResults->accept(HLT_IsoMu18_v5);

  //Single/double leptonic for ttHbb analysis
  if(HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ || HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL || HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL || 
     HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ || HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ) triggerDL = 1;
  else triggerDL = 0;
  if(HLT_Ele27_eta2p1_WPLoose_Gsf || HLT_IsoMu18 || HLT_Ele27_WP85_Gsf || HLT_IsoMu17_eta2p1) triggerSL = 1;
  else triggerSL = 0;
}

void TriggerSelector::SetBranches(){
  if(debug_)    std::cout<<"setting branches: calling AddBranch of baseTree"<<std::endl;
  AddBranch(&Trigger_decision,   "Trigger_decision");
  AddBranch(&Trigger_names,   "Trigger_names");
  AddBranch(&HLT_Ele105_CaloIdVT_GsfTrkIdT, "HLT_Ele105_CaloIdVT_GsfTrkIdT");
  AddBranch(&HLT_Ele27_eta2p1_WP75_Gsf, "HLT_Ele27_eta2p1_WP75_Gsf");
  AddBranch(&HLT_Ele27_WP85_Gsf, "HLT_Ele27_WP85_Gsf");
  AddBranch(&HLT_Ele27_eta2p1_WPLoose_Gsf, "HLT_Ele27_eta2p1_WPLoose_Gsf");
  AddBranch(&HLT_Mu50, "HLT_Mu50");
  AddBranch(&HLT_IsoMu20, "HLT_IsoMu20");
  AddBranch(&HLT_IsoMu17_eta2p1, "HLT_IsoMu17_eta2p1");
  AddBranch(&HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ, "HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ");
  AddBranch(&HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL, "HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL");
  AddBranch(&HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL, "HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL");
  AddBranch(&HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ, "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ");
  AddBranch(&HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ, "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ");
  AddBranch(&HLT_IsoMu24_eta2p1, "HLT_IsoMu24_eta2p1");
  AddBranch(&HLT_IsoMu24_eta2p1, "HLT_IsoMu18");
  AddBranch(&triggerSL, "triggerSL");
  AddBranch(&triggerDL, "triggerDL");
  if(debug_)    std::cout<<"set branches"<<std::endl;
}

void TriggerSelector::Clear(){
  Trigger_decision.clear();
  Trigger_names.clear(); 
  HLT_Ele105_CaloIdVT_GsfTrkIdT = -9999;
  HLT_Ele27_eta2p1_WP75_Gsf = -9999;
  HLT_Ele27_WP85_Gsf = -9999;
  HLT_Ele27_eta2p1_WPLoose_Gsf = -9999;
  HLT_Mu50 = -9999;
  HLT_IsoMu20 = -9999;
  HLT_IsoMu17_eta2p1 = -9999;
  HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ = -9999;
  HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL = -9999;
  HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL = -9999;
  HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ = -9999;
  HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ = -9999;
  HLT_IsoMu24_eta2p1 = -9999;
  HLT_IsoMu18 = -9999;
  triggerSL = -9999;
  triggerDL = -9999;
}

void TriggerSelector::startTrigger(edm::EventSetup const& iSetup, edm::Run const & iRun){
  bool changed(true);
  hltConfig_.init(iRun,iSetup,"HLT",changed);
}
