#include "BSMFramework/BSM3G_TNT_Maker/interface/BSM3G_TNT_Maker.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include <memory>
#include <string>
BSM3G_TNT_Maker::BSM3G_TNT_Maker(const edm::ParameterSet& iConfig):
  //now do what ever initialization is needed
  MaxN(200)
{
  debug_                 = iConfig.getParameter<bool>("debug_");
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
  _fillBJetnessinfo      = iConfig.getParameter<bool>("fillBJetnessinfo"); 
  _fillBTagReweight      = iConfig.getParameter<bool>("fillBTagReweight");
  _fillMETinfo           = iConfig.getParameter<bool>("fillMETinfo");
  _fillphotoninfo        = iConfig.getParameter<bool>("fillphotoninfo");

  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("BOOM","BOOM");

  if(_fillgeninfo)           genselector        = new GenParticleSelector("miniAOD", tree_, debug_, iConfig);
  if(_fillgenHFCategoryinfo) genhfselector      = new GenHFHadrMatchSelector("miniAOD", tree_, debug_, iConfig, consumesCollector());
  if(_filleventinfo)         eventinfoselector  = new EventInfoSelector("miniAOD", tree_, debug_, iConfig);
  if(_filltriggerinfo)       trselector         = new TriggerSelector("miniAOD", tree_, debug_, iConfig);
  if(_fillPVinfo)            pvselector         = new PVSelector("miniAOD", tree_, debug_, iConfig);
  if(_fillmuoninfo)          muselector         = new MuonSelector("miniAOD", tree_, debug_, iConfig, consumesCollector());
  if(_fillelectronpatinfo)   elpatselector      = new ElectronPatSelector("miniAOD", tree_, debug_, iConfig, consumesCollector());
  if(_filltauinfo)           tauselector        = new TauSelector("miniAOD", tree_, debug_, iConfig);
  if(_filljetinfo)           jetselector        = new JetSelector("miniAOD", tree_, debug_, iConfig);
  if(_filltthjetinfo)        tthjetselector     = new TTHJetSelector("miniAOD", tree_, debug_, iConfig, consumesCollector());
  if(_fillBoostedJetinfo)    BoostedJetselector = new BoostedJetSelector("miniAOD", tree_, debug_, iConfig);
  if(_fillTopSubJetinfo)     TopSubJetselector  = new TopSubJetSelector("miniAOD", tree_, debug_, iConfig);
  if(_fillBJetnessinfo)      BJetnessselector   = new BJetnessSelector("miniAOD", tree_, debug_, iConfig, consumesCollector());
  if(_fillBTagReweight)      btagreweight       = new BTagReweight("miniAOD", tree_, debug_, iConfig);
  if(_fillMETinfo)           metselector        = new METSelector("miniAOD", tree_, debug_, iConfig);
  if(_fillphotoninfo)        photonselector     = new PhotonSelector("miniAOD", tree_, debug_, iConfig);
}
BSM3G_TNT_Maker::~BSM3G_TNT_Maker()
{
  //do anything here that needs to be done at desctruction time
  //(e.g. close files, deallocate resources etc.)
  /*
  delete genselector; 
  delete genhfselector;
  delete eventinfoselector; 
  delete trselector;  
  delete pvselector; 
  delete muselector; 
  delete elpatselector; 
  delete tauselector; 
  delete jetselector; 
  delete tthjetselector; 
  delete BoostedJetselector; 
  delete TopSubJetselector; 
  delete BJetnessselector; 
  delete btagreweight; 
  delete metselector; 
  delete photonselector; 
  */
}
/////
//   Member functions
/////
// ------------ method called for each event  ------------
void BSM3G_TNT_Maker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace pat;
  using namespace reco;

  int savebjetnessevt = 0;
  if(_fillBJetnessinfo) BJetnessselector->Fill(iEvent, iSetup, savebjetnessevt);
  if((_fillBJetnessinfo && savebjetnessevt==1) || !_fillBJetnessinfo){
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
    if(_fillBTagReweight)      btagreweight->Fill(iEvent);
    if(_fillMETinfo)           metselector->Fill(iEvent);
    if(_fillphotoninfo)        photonselector->Fill(iEvent);
    tree_->Fill();
  }
#ifdef THIS_IS_AN_EVENT_EXAMPLE
  Handle<ExampleData> pIn;
  iEvent.getByLabel("example",pIn);
#endif
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  ESHandle<SetupData> pSetup;
  iSetup.get<SetupRecord>().get(pSetup);
#endif
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
