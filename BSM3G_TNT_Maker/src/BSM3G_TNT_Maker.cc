#include "BSMFramework/BSM3G_TNT_Maker/interface/BSM3G_TNT_Maker.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
//Trigger
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include <memory>
#include <string>       // std::string
#include <iostream>     // std::cout
#include <sstream>      // std::stringstream, std::stringbuf
BSM3G_TNT_Maker::BSM3G_TNT_Maker(const edm::ParameterSet& iConfig):
  triggerBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"))),
  //now do what ever initialization is needed
  MaxN(200)
{
  debug_                 = iConfig.getParameter<bool>("debug_");
  _is_data               = iConfig.getParameter<bool>("is_data");
  _ifevtriggers          = iConfig.getParameter<bool>("ifevtriggers"); 
  _evtriggers            = iConfig.getParameter<vector<string> >("evtriggers");
  _fillBoostedJetinfo    = iConfig.getParameter<bool>("fillBoostedJetinfo");
  _fillelectronpatinfo   = iConfig.getParameter<bool>("fillelectronpatinfo"); 
  _filleventinfo         = iConfig.getParameter<bool>("filleventinfo"); 
  _fillgeninfo           = iConfig.getParameter<bool>("fillgeninfo"); 
  _filljetinfo           = iConfig.getParameter<bool>("filljetinfo"); 
  _fillMETinfo           = iConfig.getParameter<bool>("fillMETinfo");
  _fillmuoninfo          = iConfig.getParameter<bool>("fillmuoninfo");
  _fillPVinfo            = iConfig.getParameter<bool>("fillPVinfo");  
  _fillPileupReweight    = iConfig.getParameter<bool>("fillPileupReweight");
  _filltriggerinfo       = iConfig.getParameter<bool>("filltriggerinfo");
  _filltauinfo           = iConfig.getParameter<bool>("filltauinfo");
  _filltauboostedinfo    = iConfig.getParameter<bool>("filltauboostedinfo");
  _fillphotoninfo        = iConfig.getParameter<bool>("fillphotoninfo");
  lheEventProduct_ = consumes<LHERunInfoProduct,edm::InRun>(edm::InputTag("externalLHEProducer"));
  //lheEventProduct_(  consumes<LHERunInfoProduct,edm::InRun>(edm::InputTag("externalLHEProducer"))),

  edm::Service<TFileService> fs;
  evtree_ = fs->make<TTree>("evtree","evtree");
  evtree_->Branch("eventnum",&eventnum,"eventnum/I");
  evtree_->Branch("eventnumnegative",&eventnumnegative,"eventnumnegative/I");
  tree_   = fs->make<TTree>("BOOM","BOOM");
  if(_fillBoostedJetinfo)    BoostedJetselector = new BoostedJetSelector("miniAOD", tree_, debug_, iConfig, consumesCollector());
  if(_fillelectronpatinfo)   elpatselector      = new ElectronPatSelector("miniAOD", tree_, debug_, iConfig, consumesCollector());
  if(_filleventinfo)         eventinfoselector  = new EventInfoSelector("miniAOD", tree_, debug_, iConfig, consumesCollector());
  if(_fillgeninfo)           genselector        = new GenParticleSelector("miniAOD", tree_, debug_, iConfig, consumesCollector());
  if(_filljetinfo)           jetselector        = new JetSelector("miniAOD", tree_, debug_, iConfig, consumesCollector());
  if(_fillMETinfo)           metselector        = new METSelector("miniAOD", tree_, debug_, iConfig, consumesCollector());
  if(_fillmuoninfo)          muselector         = new MuonSelector("miniAOD", tree_, debug_, iConfig, consumesCollector());
  if(_fillPVinfo)            pvselector         = new PVSelector("miniAOD", tree_, debug_, iConfig, consumesCollector());
  if(_fillPileupReweight)    pileupreweight     = new PileupReweight("miniAOD", tree_, debug_, iConfig, consumesCollector());
  if(_filltauinfo)           tauselector        = new TauSelector("miniAOD", tree_, debug_, iConfig, consumesCollector());
  if(_filltauboostedinfo)    tauboostedselector = new TauboostedSelector("miniAOD", tree_, debug_, iConfig, consumesCollector());
  if(_filltriggerinfo)       trselector         = new TriggerSelector("miniAOD", tree_, debug_, iConfig);
  if(_fillphotoninfo)        photonselector     = new PhotonSelector("miniAOD", tree_, debug_, iConfig);
}
BSM3G_TNT_Maker::~BSM3G_TNT_Maker()
{
  //do anything here that needs to be done at desctruction time
  //(e.g. close files, deallocate resources etc.)
}
/////
//   Member functions
/////
// ------------ method called for each event  ------------
void BSM3G_TNT_Maker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //Namespace
  using namespace edm;
  using namespace pat;
  using namespace reco;
  //Event info for all the events you read
  eventnum = -1;
  eventnum = iEvent.id().event();
  eventnumnegative = 1;
  if(!_is_data){
    edm::Handle<GenEventInfoProduct> genEvtInfo;
    iEvent.getByLabel("generator",genEvtInfo);
    eventnumnegative = (genEvtInfo->weight())/abs(genEvtInfo->weight());
  }
  evtree_->Fill();
  //Require trigger on the event
  bool evtriggered = false;
  if(_ifevtriggers){
    edm::Handle<edm::TriggerResults> triggerBits;
    iEvent.getByToken(triggerBits_, triggerBits);
    const edm::TriggerNames &trigNames = iEvent.triggerNames(*triggerBits);
    for(uint tb = 0; tb<triggerBits->size(); tb++){
      for(uint tn = 0; tn<_evtriggers.size(); tn++){
        if(strstr(trigNames.triggerName(tb).c_str(),_evtriggers[tn].c_str()) && triggerBits->accept(tb)){
          evtriggered = true;
          break;
        } 
      }
    }
  }
  //Call classes
  if((_ifevtriggers && evtriggered) || !_ifevtriggers){
    if(_fillBoostedJetinfo)    BoostedJetselector->Fill(iEvent);
    if(_fillelectronpatinfo)   elpatselector->Fill(iEvent,iSetup); 
    if(_filleventinfo)         eventinfoselector->Fill(iEvent);
    if(_fillgeninfo)           genselector->Fill(iEvent); 
    if(_filljetinfo)           jetselector->Fill(iEvent);
    if(_fillMETinfo)           metselector->Fill(iEvent);
    if(_fillmuoninfo)          muselector->Fill(iEvent,iSetup);
    if(_fillPVinfo)            pvselector->Fill(iEvent); 
    if(_fillPileupReweight)    pileupreweight->Fill(iEvent);
    if(_filltauinfo)           tauselector->Fill(iEvent,iSetup); 
    if(_filltauboostedinfo)    tauboostedselector->Fill(iEvent,iSetup); 
    if(_filltriggerinfo)       trselector->Fill(iEvent,iSetup);
    if(_fillphotoninfo)        photonselector->Fill(iEvent);
    tree_->Fill();
  }
}
// ------------ method called once each job just before starting event loop  ------------
void BSM3G_TNT_Maker::beginJob()
{
}
// ------------ method called once each job just after ending the event loop  ------------
void BSM3G_TNT_Maker::endJob() 
{
}
// ------------ method called when starting to processes a run  ------------
void BSM3G_TNT_Maker::beginRun(edm::Run const & iRun, edm::EventSetup const& iSetup)
{
  if( _filltriggerinfo) trselector->startTrigger(iSetup, iRun);
}
// ------------ method called when ending the processing of a run  ------------
void BSM3G_TNT_Maker::endRun(edm::Run const & iRun, edm::EventSetup const&)
{
  //edm::Handle<LHERunInfoProduct> run;
  //typedef std::vector<LHERunInfoProduct::Header>::const_iterator headers_const_iterator;
  //iRun.getByToken(lheEventProduct_, run);
  //LHERunInfoProduct myLHERunInfoProduct = *(run.product());
  //for (headers_const_iterator iter=myLHERunInfoProduct.headers_begin(); iter!=myLHERunInfoProduct.headers_end(); iter++){
  //  std::cout << iter->tag() << std::endl;
  //  std::string weightTag = "initrwgt";
  //  if (iter->tag() != weightTag) continue;
  //  std::vector<std::string> lines = iter->lines();
  //  for (unsigned int iLine = 0; iLine<lines.size(); iLine++) {
  //    std::cout << lines.at(iLine);
  //  }
  //}
}
// ------------ method called when starting to processes a luminosity block  ------------
//void BSM3G_TNT_Maker::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
//{
//}
// ------------ method called when ending the processing of a luminosity block  ------------
//void BSM3G_TNT_Maker::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
//{
//}
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void BSM3G_TNT_Maker::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(BSM3G_TNT_Maker);
