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
  bjetnessselfilter      = iConfig.getParameter<bool>("bjetnessselfilter");
  _is_data               = iConfig.getParameter<bool>("is_data");
  _ifevtriggers          = iConfig.getParameter<bool>("ifevtriggers");
  _evtriggers            = iConfig.getParameter<vector<string> >("evtriggers");
  _fillgeninfo           = iConfig.getParameter<bool>("fillgeninfo");
  _fillgenHFCategoryinfo = iConfig.getParameter<bool>("fillgenHFCategoryinfo");
  _filleventinfo         = iConfig.getParameter<bool>("filleventinfo");
  _filltriggerinfo       = iConfig.getParameter<bool>("filltriggerinfo");
  _fillPVinfo            = iConfig.getParameter<bool>("fillPVinfo");
  _fillmuoninfo          = iConfig.getParameter<bool>("fillmuoninfo");
  _fillelectronpatinfo   = iConfig.getParameter<bool>("fillelectronpatinfo");
  _filltauinfo           = iConfig.getParameter<bool>("filltauinfo");
  _filljetinfo           = iConfig.getParameter<bool>("filljetinfo");
  _filltthjetinfo        = iConfig.getParameter<bool>("filltthjetinfo");
  _fillBoostedJetinfo    = iConfig.getParameter<bool>("fillBoostedJetinfo");
  _fillTopSubJetinfo     = iConfig.getParameter<bool>("fillTopSubJetinfo");
  _fillTauJetnessinfo    = iConfig.getParameter<bool>("fillTauJetnessinfo");
  _fillBJetnessinfo      = iConfig.getParameter<bool>("fillBJetnessinfo");
  _fillBJetnessFVinfo    = iConfig.getParameter<bool>("fillBJetnessFVinfo");
  _fillBTagReweight      = iConfig.getParameter<bool>("fillBTagReweight");
  _fillPileupReweight    = iConfig.getParameter<bool>("fillPileupReweight");
  _fillMETinfo           = iConfig.getParameter<bool>("fillMETinfo");
  _fillphotoninfo        = iConfig.getParameter<bool>("fillphotoninfo");
  analysisFilter         = iConfig.getParameter<bool>("analysisFilter");

  edm::Service<TFileService> fs;
  evtree_ = fs->make<TTree>("evtree","evtree");
  evtree_->Branch("eventnum",&eventnum,"eventnum/I");
  evtree_->Branch("eventnumnegative",&eventnumnegative,"eventnumnegative/I");
  tree_   = fs->make<TTree>("BOOM","BOOM");
  cout << "BSM3G_TNT_Maker:BSM3G_TNT_Maker >> Number of entries in TTree = " << tree_->GetEntries() << endl;

  //> Register helper classes to the framework.
  if(_fillgeninfo)           genselector        = new GenParticleSelector("miniAOD", tree_, debug_, iConfig, consumesCollector());
  if(_fillgenHFCategoryinfo) genhfselector      = new GenHFHadrMatchSelector("miniAOD", tree_, debug_, iConfig, consumesCollector());
  if(_filleventinfo)         eventinfoselector  = new EventInfoSelector("miniAOD", tree_, debug_, iConfig, consumesCollector());
  if(_filltriggerinfo)       trselector         = new TriggerSelector("miniAOD", tree_, debug_, iConfig);
  if(_fillPVinfo)            pvselector         = new PVSelector("miniAOD", tree_, debug_, iConfig, consumesCollector());
  if(_fillmuoninfo)          muselector         = new MuonSelector("miniAOD", tree_, debug_, iConfig, consumesCollector());
  if(_fillelectronpatinfo)   elpatselector      = new ElectronPatSelector("miniAOD", tree_, debug_, iConfig, consumesCollector());
  if(_filltauinfo)           tauselector        = new TauSelector("miniAOD", tree_, debug_, iConfig, consumesCollector());
  if(_filljetinfo)           jetselector        = new JetSelector("miniAOD", tree_, debug_, iConfig, consumesCollector());
  if(_filltthjetinfo)        tthjetselector     = new TTHJetSelector("miniAOD", tree_, debug_, iConfig, consumesCollector());
  if(_fillBoostedJetinfo)    BoostedJetselector = new BoostedJetSelector("miniAOD", tree_, debug_, iConfig, consumesCollector());
  if(_fillTopSubJetinfo)     TopSubJetselector  = new TopSubJetSelector("miniAOD", tree_, debug_, iConfig);
  if(_fillTauJetnessinfo)    TauJetnessselector = new TauJetnessSelector("miniAOD", tree_, debug_, iConfig, consumesCollector());
  if(_fillBJetnessinfo)      BJetnessselector   = new BJetnessSelector("miniAOD", tree_, debug_, iConfig, consumesCollector());
  if(_fillBJetnessFVinfo)    BJetnessFVselector = new BJetnessFVSelector("miniAOD", tree_, debug_, iConfig, consumesCollector());
  if(_fillBTagReweight)      btagreweight       = new BTagReweight("miniAOD", tree_, debug_, iConfig, consumesCollector());
  if(_fillPileupReweight)    pileupreweight     = new PileupReweight("miniAOD", tree_, debug_, iConfig, consumesCollector());
  if(_fillMETinfo)           metselector        = new METSelector("miniAOD", tree_, debug_, iConfig, consumesCollector());
  if(_fillphotoninfo)        photonselector     = new PhotonSelector("miniAOD", tree_, debug_, iConfig);
  if(analysisFilter)         {TTHbb_eventselector = new TTHbb_eventSelector("miniAOD", tree_, debug_, iConfig, consumesCollector()); TTHbb_eventselector->SetBranches();}
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
  //if (eventnum == 315650){
  //  if(debug_) cout << ">>>>>>>>>>>>>>>>>>>>>>>>>> " << eventnum << " <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< I'M HERE!" << endl;
  //}
  eventnumnegative = 1;
  if(!_is_data){
    edm::Handle<GenEventInfoProduct> genEvtInfo;
    iEvent.getByLabel("generator",genEvtInfo);
    eventnumnegative = (genEvtInfo->weight())/abs(genEvtInfo->weight());
  }
  evtree_->Fill();
  bool evtriggered = false;
  if(_ifevtriggers){//Flag to reqeust a trigger was fired in order to save event
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
  // Any if statment here works as the ONLY form of event selection.
  // Otherwise there should be no reason an event is not saved.
  if((_ifevtriggers && evtriggered) || !_ifevtriggers){//Require trigger on the event if flag, otherwise save without requirement.
    bjetnesssel_filter = 0;
    if(_fillBJetnessinfo)      BJetnessselector->Fill(iEvent, iSetup, bjetnesssel_filter);
    if((bjetnessselfilter && bjetnesssel_filter==1) || !bjetnessselfilter){
      if(analysisFilter)         TTHbb_eventselector->Fill(iEvent,iSetup);
      if(_fillBJetnessFVinfo)    BJetnessFVselector->Fill(iEvent, iSetup);
      if(_fillgeninfo)           genselector->Fill(iEvent);
      if(_fillgenHFCategoryinfo) genhfselector->Fill(iEvent);
      if(_filleventinfo)         eventinfoselector->Fill(iEvent);
      if(_filltriggerinfo)       trselector->Fill(iEvent,iSetup);
      if(_fillPVinfo)            pvselector->Fill(iEvent);
      if(_fillmuoninfo)          muselector->Fill(iEvent,iSetup);
      if(_fillelectronpatinfo)   elpatselector->Fill(iEvent,iSetup);
      if(_filltauinfo)           tauselector->Fill(iEvent,iSetup);
      if(_filljetinfo)           jetselector->Fill(iEvent);
      if(_filltthjetinfo)        tthjetselector->Fill(iEvent,iSetup);
      if(_fillBoostedJetinfo)    BoostedJetselector->Fill(iEvent);
      if(_fillTopSubJetinfo)     TopSubJetselector->Fill(iEvent);
      if(_fillTauJetnessinfo)    TauJetnessselector->Fill(iEvent, iSetup);
      if(_fillBTagReweight)      btagreweight->Fill(iEvent);
      if(_fillPileupReweight)    pileupreweight->Fill(iEvent);
      if(_fillMETinfo)           metselector->Fill(iEvent);
      if(_fillphotoninfo)        photonselector->Fill(iEvent);
      tree_->Fill();
    }
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
//void BSM3G_TNT_Maker::endRun(edm::Run const&, edm::EventSetup const&)
//{
//}
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
