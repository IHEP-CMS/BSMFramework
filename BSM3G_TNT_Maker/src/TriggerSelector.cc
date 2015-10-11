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
  if(trigResults->accept(trigNames.triggerIndex("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v1"))) triggerDL = 1;
  else if(trigResults->accept(trigNames.triggerIndex("HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v1"))) triggerDL = 1;
  else if(trigResults->accept(trigNames.triggerIndex("HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v1"))) triggerDL = 1;
  else if(trigResults->accept(trigNames.triggerIndex("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v1"))) triggerDL = 1;
  else if(trigResults->accept(trigNames.triggerIndex("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v1"))) triggerDL = 1;
  else triggerDL = 0;
  if(trigResults->accept(trigNames.triggerIndex("HLT_Ele27_eta2p1_WP85_Gsf_HT200_v1"))) triggerSL = 1;
  else if(trigResults->accept(trigNames.triggerIndex("HLT_IsoMu24_eta2p1_v1"))) triggerSL = 1;
  else triggerSL = 0;
}

void TriggerSelector::SetBranches(){
  if(debug_)    std::cout<<"setting branches: calling AddBranch of baseTree"<<std::endl;
  
  AddBranch(&Trigger_decision,   "Trigger_decision");
  AddBranch(&Trigger_names,   "Trigger_names");
  AddBranch(&triggerSL, "triggerSL");
  AddBranch(&triggerDL, "triggerDL");

  if(debug_)    std::cout<<"set branches"<<std::endl;
}

void TriggerSelector::Clear(){
  Trigger_decision.clear();
  Trigger_names.clear();  
  triggerSL = -9999;
  triggerDL = -9999;

}

void TriggerSelector::startTrigger(edm::EventSetup const& iSetup, edm::Run const & iRun){
  bool changed(true);
  hltConfig_.init(iRun,iSetup,"HLT",changed);
}
