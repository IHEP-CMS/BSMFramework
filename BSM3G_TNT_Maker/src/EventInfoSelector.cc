#include "BSMFramework/BSM3G_TNT_Maker/interface/EventInfoSelector.h"
EventInfoSelector::EventInfoSelector(std::string name, TTree* tree, bool debug, const pset& iConfig, edm::ConsumesCollector && ic):
  baseTree(name,tree,debug)
{
  genEvtInfo_   = ic.consumes<GenEventInfoProduct>(edm::InputTag("generator"));
  rhopogHandle_ = ic.consumes<double>(edm::InputTag("fixedGridRhoFastjetAll"));
  rhotthHandle_ = ic.consumes<double>(edm::InputTag("fixedGridRhoFastjetCentralNeutral"));
  _is_data = iConfig.getParameter<bool>("is_data");
  if(debug) std::cout<<"in EventInfoSelector constructor"<<std::endl;
  SetBranches();
}
EventInfoSelector::~EventInfoSelector(){
  delete tree_;
}
void EventInfoSelector::Fill(const edm::Event& iEvent){
  Initialise();
  EVENT_event_     = iEvent.id().event();
  EVENT_run_       = iEvent.id().run();
  EVENT_lumiBlock_ = iEvent.id().luminosityBlock();
  EVENT_genWeight_ = 1;
  edm::Handle<GenEventInfoProduct> genEvtInfo;
  iEvent.getByToken(genEvtInfo_,genEvtInfo);
  if(!_is_data){
    EVENT_genWeight_ = genEvtInfo->weight();
  }
  edm::Handle<double> rhopogHandle;
  iEvent.getByToken(rhopogHandle_,rhopogHandle);
  double rhopog = *rhopogHandle;
  edm::Handle<double> rhotthHandle;
  iEvent.getByToken(rhotthHandle_,rhotthHandle);
  double rhotth = *rhotthHandle;
  EVENT_rhopog_ = rhopog;
  EVENT_rhotth_ = rhotth;
}
void EventInfoSelector::SetBranches(){
  if(debug_) std::cout<<"setting branches: calling AddBranch of baseTree"<<std::endl;
  AddBranch(&EVENT_event_     ,"EVENT_event");
  AddBranch(&EVENT_run_       ,"EVENT_run");
  AddBranch(&EVENT_lumiBlock_ ,"EVENT_lumiBlock");
  AddBranch(&EVENT_genWeight_ ,"EVENT_genWeight");
  AddBranch(&EVENT_rhopog_    ,"EVENT_rhopog");
  AddBranch(&EVENT_rhotth_    ,"EVENT_rhotth");
}
void EventInfoSelector::Initialise(){
  EVENT_event_      = -9999;
  EVENT_run_        = -9999; 
  EVENT_lumiBlock_  = -9999;
  EVENT_genWeight_  = -9999;
  EVENT_rhopog_     = -9999;
  EVENT_rhotth_     = -9999; 
}
