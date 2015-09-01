#include "NtupleMaker/BSM3G_TNT_Maker/interface/EventInfoSelector.h"
EventInfoSelector::EventInfoSelector(std::string name, TTree* tree, bool debug, const pset& iConfig):baseTree(name,tree,debug){
  if(debug) std::cout<<"in EventInfoSelector constructor"<<std::endl;
  _is_data = iConfig.getParameter<bool>("is_data");
  SetBranches();
}
EventInfoSelector::~EventInfoSelector(){
  delete tree_;
}
void EventInfoSelector::Fill(const edm::Event& iEvent){
  EVENT_event_     = iEvent.id().event();
  EVENT_run_       = iEvent.id().run();
  EVENT_lumiBlock_ = iEvent.id().luminosityBlock();
  EVENT_genWeight_ = 1;
  edm::Handle<GenEventInfoProduct> genEvtInfo;
  iEvent.getByLabel("generator",genEvtInfo);
  if(!_is_data){
    EVENT_genWeight_ = genEvtInfo.product()->weight();
  }
}
void EventInfoSelector::SetBranches(){
  if(debug_) std::cout<<"setting branches: calling AddBranch of baseTree"<<std::endl;
  AddBranch(&EVENT_event_,    "EVENT_event");
  AddBranch(&EVENT_run_,      "EVENT_run");
  AddBranch(&EVENT_lumiBlock_,"EVENT_lumiBlock");
  AddBranch(&EVENT_genWeight_,"EVENT_genWeight");
}
