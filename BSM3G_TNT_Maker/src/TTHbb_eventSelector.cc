#include "BSMFramework/BSM3G_TNT_Maker/interface/TTHbb_eventSelector.h"
#include <iostream>
#include <fstream>
#include <ios>
using namespace edm;
using namespace std;


// Code uses member initialisation e.g. Foo(int num): bar(num) {};
// instead of member assignment e.g. Foo(int num){bar = num;} =>> FASTER!!
TTHbb_eventSelector::TTHbb_eventSelector(string treeName, TTree* tree, bool debug, const pset& iConfig, edm::ConsumesCollector && ic):
  baseTree(treeName,tree,debug),
  // Must register that module/helper will be making data access requests. Can be done
  // in constructor or 'beginJob' method:
  //
  //        e.g. edm::EDGetTokenT<T> consumes<T>(edm::InputTag const&)
  //
  // 'T' = c++ class used to get data.
  // Above e.g. tells framework the data will be used when the module is called.
  // ConsumesCollector gathers all consumes information for EDConsumerBase.
  // Can be passed around to register data on behalf of module.
  electronVetoIdMapToken_(ic.consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electronVetoIdMap"))),
  electronLooseIdMapToken_(ic.consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electronLooseIdMap"))),
  electronMediumIdMapToken_(ic.consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electronMediumIdMap"))),
  electronTightIdMapToken_(ic.consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electronTightIdMap"))),
  eleMVATrigIdMapToken_(ic.consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleMVATrigIdMap"))),
  eleMVAnonTrigIdMap_(ic.consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleMVAnonTrigIdMap"))),
  eleMVATrigwp90IdMap_(ic.consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleMVATrigwp90IdMap"))),
  eleMVAnonTrigwp90IdMap_(ic.consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleMVAnonTrigwp90IdMap"))),
  eleHEEPIdMapToken_(ic.consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleHEEPIdMap"))),
  elemvaValuesMapToken_nonTrig_(ic.consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("elemvaValuesMap_nonTrig"))),
  elemvaCategoriesMapToken_nonTrig_(ic.consumes<edm::ValueMap<int> >(iConfig.getParameter<edm::InputTag>("elemvaCategoriesMap_nonTrig"))),
  elemvaValuesMapToken_Trig_(ic.consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("elemvaValuesMap_Trig"))),
  elemvaCategoriesMapToken_Trig_(ic.consumes<edm::ValueMap<int> >(iConfig.getParameter<edm::InputTag>("elemvaCategoriesMap_Trig"))),
  electron_pat_(ic.consumes<edm::View<pat::Electron> >(iConfig.getParameter<edm::InputTag>("patElectrons"))),
  muon_h_(ic.consumes<edm::View<pat::Muon> >(iConfig.getParameter<edm::InputTag>("muons"))),
  jets_(ic.consumes<pat::JetCollection >(iConfig.getParameter<edm::InputTag>("jets"))),
  rhopogHandle_(ic.consumes<double>(edm::InputTag("fixedGridRhoFastjetAll"))),
  rhoJERHandle_(ic.consumes<double>(edm::InputTag("fixedGridRhoAll"))),
  jecPayloadNamesAK4PFchsMC1_(iConfig.getParameter<edm::FileInPath>("jecPayloadNamesAK4PFchsMC1")),
  jecPayloadNamesAK4PFchsMC2_(iConfig.getParameter<edm::FileInPath>("jecPayloadNamesAK4PFchsMC2")),
  jecPayloadNamesAK4PFchsMC3_(iConfig.getParameter<edm::FileInPath>("jecPayloadNamesAK4PFchsMC3")),
  jecPayloadNamesAK4PFchsMCUnc_(iConfig.getParameter<edm::FileInPath>("jecPayloadNamesAK4PFchsMCUnc")),
  jecPayloadNamesAK4PFchsDATA1_(iConfig.getParameter<edm::FileInPath>("jecPayloadNamesAK4PFchsDATA1")),
  jecPayloadNamesAK4PFchsDATA2_(iConfig.getParameter<edm::FileInPath>("jecPayloadNamesAK4PFchsDATA2")),
  jecPayloadNamesAK4PFchsDATA3_(iConfig.getParameter<edm::FileInPath>("jecPayloadNamesAK4PFchsDATA3")),
  jecPayloadNamesAK4PFchsDATA4_(iConfig.getParameter<edm::FileInPath>("jecPayloadNamesAK4PFchsDATA4")),
  jecPayloadNamesAK4PFchsDATAUnc_(iConfig.getParameter<edm::FileInPath>("jecPayloadNamesAK4PFchsDATAUnc")),
  jerAK4PFchs_(iConfig.getParameter<edm::FileInPath>("jerAK4PFchs").fullPath()),
  jerAK4PFchsSF_(iConfig.getParameter<edm::FileInPath>("jerAK4PFchsSF").fullPath()),
  vtx_h_(ic.consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
  _vtx_ndof_min(iConfig.getParameter<int>("vtx_ndof_min")),
  _vtx_rho_max(iConfig.getParameter<int>("vtx_rho_max")),
  _vtx_position_z_max(iConfig.getParameter<double>("vtx_position_z_max")),
  _is_data(iConfig.getParameter<bool>("is_data"))
{
  if(!_is_data) prunedGenToken_ = ic.consumes<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("pruned"));
  JECInitialisation();
  //SetBranches();
}
TTHbb_eventSelector::~TTHbb_eventSelector(){
  //delete TTree created by baseTree constructor:
  delete tree_;
}

void TTHbb_eventSelector::JECInitialisation(){
  //AK4chs - MC: Get the factorized jet corrector parameters.
  std::vector<std::string> jecPayloadNamesAK4PFchsMC_;
  // Get path to JEC .txt from config and add into vector.
  jecPayloadNamesAK4PFchsMC_.push_back(jecPayloadNamesAK4PFchsMC1_.fullPath());
  jecPayloadNamesAK4PFchsMC_.push_back(jecPayloadNamesAK4PFchsMC2_.fullPath());
  jecPayloadNamesAK4PFchsMC_.push_back(jecPayloadNamesAK4PFchsMC3_.fullPath());
  std::vector<JetCorrectorParameters> vParAK4PFchsMC;
  // Loop over vector of JEC files and create JetCorrectorParameters objects for each.
  // const_iterator used: we really dont want to be able to edit the original .txt file!!
  for ( std::vector<std::string>::const_iterator payloadBegin = jecPayloadNamesAK4PFchsMC_.begin(),
          payloadEnd = jecPayloadNamesAK4PFchsMC_.end(), ipayload = payloadBegin; ipayload != payloadEnd; ++ipayload ) {
    JetCorrectorParameters pars(*ipayload);
    vParAK4PFchsMC.push_back(pars);
  }
  jecAK4PFchsMC_    = boost::shared_ptr<FactorizedJetCorrector>  ( new FactorizedJetCorrector(vParAK4PFchsMC) );
  // Get JEC uncertainty
  jecAK4PFchsMCUnc_ = boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(jecPayloadNamesAK4PFchsMCUnc_.fullPath()) );
  //AK4chs - DATA: Get the factorized jet corrector parameters.
  std::vector<std::string> jecPayloadNamesAK4PFchsDATA_;
  jecPayloadNamesAK4PFchsDATA_.push_back(jecPayloadNamesAK4PFchsDATA1_.fullPath());
  jecPayloadNamesAK4PFchsDATA_.push_back(jecPayloadNamesAK4PFchsDATA2_.fullPath());
  jecPayloadNamesAK4PFchsDATA_.push_back(jecPayloadNamesAK4PFchsDATA3_.fullPath());
  jecPayloadNamesAK4PFchsDATA_.push_back(jecPayloadNamesAK4PFchsDATA4_.fullPath());
  std::vector<JetCorrectorParameters> vParAK4PFchsDATA;
  for ( std::vector<std::string>::const_iterator payloadBegin = jecPayloadNamesAK4PFchsDATA_.begin(),
          payloadEnd = jecPayloadNamesAK4PFchsDATA_.end(), ipayload = payloadBegin; ipayload != payloadEnd; ++ipayload ) {
    JetCorrectorParameters pars(*ipayload);
    vParAK4PFchsDATA.push_back(pars);
  }
  jecAK4PFchsDATA_    = boost::shared_ptr<FactorizedJetCorrector>  ( new FactorizedJetCorrector(vParAK4PFchsDATA) );
  jecAK4PFchsDATAUnc_ = boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(jecPayloadNamesAK4PFchsDATAUnc_.fullPath()) );
}

void TTHbb_eventSelector::SetBranches(){
  if(debug_) std::cout<<"TTHbb_eventSelector::SetBranches >>> Setting branches: calling AddBranch of baseTree"<<std::endl;
  AddBranch(&is_singlelep,"is_singlelep");
  AddBranch(&is_dilepton,"is_dilepton");
}

// If 'good' vertex requirements are not met, returns false.
bool TTHbb_eventSelector::isGoodVertex(const reco::Vertex& vtx){
  // Fake primary vertex is the same as the beam spot.
  // Often found if not enough quality tracks or no sensible vertex is reconstructed.
  if(vtx.isFake())                                   return false;
  // # Degrees of freedom: reject vertices effectively consisting of one track.
  // Indicator of fit quality.
  if(vtx.ndof()<_vtx_ndof_min)                       return false;

  if(vtx.position().Rho()>_vtx_rho_max)              return false;
  if(fabs(vtx.position().Z()) > _vtx_position_z_max) return false;
  return true;
}

double TTHbb_eventSelector::rel_iso_dbc_mu(const pat::Muon& lepton){
  return((lepton.pfIsolationR04().sumChargedHadronPt + max(lepton.pfIsolationR04().sumNeutralHadronEt + lepton.pfIsolationR04().sumPhotonEt - 0.5 * lepton.pfIsolationR04().sumPUPt,0.0))/lepton.pt()
        );
}

double TTHbb_eventSelector::rel_iso_dbc_ele(const pat::Electron& el, double rhopog){
  double SumChHadPt       = el.pfIsolationVariables().sumChargedHadronPt;
  double SumNeuHadEt      = el.pfIsolationVariables().sumNeutralHadronEt;
  double SumPhotonEt      = el.pfIsolationVariables().sumPhotonEt;
  double EffArea          = get_effarea(el.superCluster()->position().eta());
  double SumNeutralCorrEt = std::max( 0.0, SumNeuHadEt+SumPhotonEt - rhopog*EffArea );
  double relIsoRhoEA = (SumChHadPt + SumNeutralCorrEt)/el.pt();
  return relIsoRhoEA;
}

double TTHbb_eventSelector::get_effarea(double eta){
  double effarea = -1;
  if(abs(eta) < 1.0)        effarea = 0.1752;
  else if(abs(eta) < 1.479) effarea = 0.1862;
  else if(abs(eta) < 2.0)   effarea = 0.1411;
  else if(abs(eta) < 2.2)   effarea = 0.1534;
  else if(abs(eta) < 2.3)   effarea = 0.1903;
  else if(abs(eta) < 2.4)   effarea = 0.2243;
  else                      effarea = 0.2687;
  return effarea;
}

// If 'loose' muon requirements are not met, return false.
bool TTHbb_eventSelector::is_loose_muon(const pat::Muon& mu, const reco::Vertex& vtx){
  bool isloosemu = false;
  if(mu.pt()>15 && TMath::Abs(mu.eta()) < 2.4 && mu.isTightMuon(vtx) && rel_iso_dbc_mu(mu) < 0.25) isloosemu = true;
  return isloosemu;
}

// If 'tight' muon requirements are not met, return false.
bool TTHbb_eventSelector::is_tight_muon(const pat::Muon& mu, const reco::Vertex& vtx){
  bool isloosemu = false;
  if(mu.pt()>26 && TMath::Abs(mu.eta()) < 2.1 && mu.isTightMuon(vtx) && rel_iso_dbc_mu(mu) < 0.15) isloosemu = true;
  return isloosemu;
}

bool TTHbb_eventSelector::is_loose_electron(const pat::Electron& ele, double rhopog){
  bool isele = false;
  if(ele.pt()>15 && TMath::Abs(ele.eta())<2.4 &&
     !(fabs(ele.superCluster()->position().eta()) > 1.4442 && fabs(ele.superCluster()->position().eta()) < 1.5660)){//Check if in crack.
     isele = true;
  }
  return isele;
}
bool TTHbb_eventSelector::is_tight_electron(const pat::Electron& ele, double rhopog){
  bool isele = false;
  if(ele.pt()>30 && TMath::Abs(ele.eta())<2.1 &&
     !(fabs(ele.superCluster()->position().eta()) > 1.4442 && fabs(ele.superCluster()->position().eta()) < 1.5660)){//check if in crack.
     isele = true;
  }
  return isele;
}

//Look for loose electron
//Require good jets
//This function has to be updated in order to select good jet using the TTHbb definition
bool TTHbb_eventSelector::is_good_jet(const pat::Jet &j,double rho, double rhoJER, int vtxsize, double channel_jpt_cut){
  bool isgoodjet = true;
  //Jet Energy Corrections and Uncertainties
  double corrAK4PFchs     = 1;
  reco::Candidate::LorentzVector uncorrJetAK4PFchs = j.correctedP4(0);
  if(!_is_data){
    jecAK4PFchsMC_->setJetEta( uncorrJetAK4PFchs.eta()    );
    jecAK4PFchsMC_->setJetPt ( uncorrJetAK4PFchs.pt()     );
    jecAK4PFchsMC_->setJetE  ( uncorrJetAK4PFchs.energy() );
    jecAK4PFchsMC_->setRho	( rho  );
    jecAK4PFchsMC_->setNPV	( vtxsize  );
    jecAK4PFchsMC_->setJetA  ( j.jetArea()	     );
    corrAK4PFchs = jecAK4PFchsMC_->getCorrection();
  } else {
    jecAK4PFchsDATA_->setJetEta( uncorrJetAK4PFchs.eta()    );
    jecAK4PFchsDATA_->setJetPt ( uncorrJetAK4PFchs.pt()     );
    jecAK4PFchsDATA_->setJetE  ( uncorrJetAK4PFchs.energy() );
    jecAK4PFchsDATA_->setRho	( rho  );
    jecAK4PFchsDATA_->setNPV	( vtxsize  );
    jecAK4PFchsDATA_->setJetA  ( j.jetArea()	     );
    corrAK4PFchs = jecAK4PFchsDATA_->getCorrection();
  }
  float JERScaleFactor     = 1;
  float JERScaleFactorUP   = 1;
  float JERScaleFactorDOWN = 1;
  if(!_is_data) GetJER(j, corrAK4PFchs, rhoJER, true, JERScaleFactor, JERScaleFactorUP, JERScaleFactorDOWN);
  //Acceptance
  double jetpt = (j.correctedJet("Uncorrected").pt()*corrAK4PFchs*JERScaleFactor);
  if(jetpt < channel_jpt_cut)        isgoodjet = false; //Please note that this requirement is for the SL channel, while for DL channel we require pT > 20!
  if(fabs(j.eta())>2.4) isgoodjet = false;
  //ID requirements
  if(j.neutralHadronEnergyFraction() >= 0.99) isgoodjet = false;
  if(j.chargedEmEnergyFraction()     >= 0.99) isgoodjet = false;
  if(j.neutralEmEnergyFraction()     >= 0.99) isgoodjet = false;
  if(j.numberOfDaughters()           <= 1)    isgoodjet = false;
  if(j.chargedHadronEnergyFraction() <= 0.0)  isgoodjet = false;
  if(j.chargedMultiplicity()         <= 0.0)  isgoodjet = false;
  //cout<<setw(20)<<"Jet pt,eta,phi"<<setw(20)<<jetpt<<setw(20)<<j.eta()<<setw(20)<<j.phi()<<endl;
  return isgoodjet;
}

void TTHbb_eventSelector::GetJER(pat::Jet jet, float JesSF, float rhoJER, bool AK4PFchs, float &JERScaleFactor, float &JERScaleFactorUP, float &JERScaleFactorDOWN){
  if(!jet.genJet()) return;
  double jetEta=fabs(jet.eta());
  double cFactorJER = 1.0;
  double cFactorJERdown = 1.0;
  double cFactorJERup = 1.0;
  //https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution#JER_Scaling_factors_and_Uncertai
  if( jetEta<0.5 ){
    cFactorJER = 1.109;
    cFactorJERdown = 1.109-0.008;
    cFactorJERup   = 1.109+0.008;
  } else if( jetEta<0.8 ){
    cFactorJER = 1.138;
    cFactorJERdown = 1.138-0.013;
    cFactorJERup   = 1.138+0.013;
  } else if( jetEta<1.1 ){
    cFactorJER = 1.114;
    cFactorJERdown = 1.114-0.013;
    cFactorJERup   = 1.114+0.013;
  } else if( jetEta<1.3 ){
    cFactorJER = 1.123;
    cFactorJERdown = 1.123-0.024;
    cFactorJERup   = 1.123+0.024;
  } else if( jetEta<1.7 ){
    cFactorJER = 1.084;
    cFactorJERdown = 1.084-0.011;
    cFactorJERup   = 1.084+0.011;
  } else if( jetEta<1.9 ){
    cFactorJER = 1.082;
    cFactorJERdown = 1.082-0.035;
    cFactorJERup   = 1.082+0.035;
  } else if( jetEta<2.1 ){
    cFactorJER = 1.140;
    cFactorJERdown = 1.140-0.047;
    cFactorJERup   = 1.140+0.047;
  } else if( jetEta<2.3 ){
    cFactorJER = 1.067;
    cFactorJERdown = 1.067-0.053;
    cFactorJERup   = 1.067+0.053;
  } else if( jetEta<2.5 ){
    cFactorJER = 1.177;
    cFactorJERdown = 1.177-0.041;
    cFactorJERup   = 1.177+0.041;
  } else if( jetEta<2.8 ){
    cFactorJER = 1.364;
    cFactorJERdown = 1.364-0.039;
    cFactorJERup   = 1.364+0.039;
  } else if( jetEta<3.0 ){
    cFactorJER = 1.857;
    cFactorJERdown = 1.857-0.071;
    cFactorJERup   = 1.857+0.071;
  } else if( jetEta<3.2 ){
    cFactorJER = 0.328;
    cFactorJERdown = 0.328-0.022;
    cFactorJERup   = 0.328+0.022;
  } else if( jetEta<5.0 ){
    cFactorJER = 1.160;
    cFactorJERdown = 1.160-0.029;
    cFactorJERup   = 1.160+0.029;
  }
  //double recoJetPt = jet.pt();//(jet.correctedJet("Uncorrected").pt())*JesSF;
  double recoJetPt = (jet.correctedJet("Uncorrected").pt())*JesSF;
  double genJetPt  = jet.genJet()->pt();
  double diffPt    = recoJetPt - genJetPt;
  JME::JetResolution resolution;
  JME::JetResolutionScaleFactor res_sf;
  if(AK4PFchs){
    resolution = JME::JetResolution(jerAK4PFchs_);
    res_sf = JME::JetResolutionScaleFactor(jerAK4PFchsSF_);
  } else {
    //resolution = JME::JetResolution(jerAK4PFPuppi_);
    //res_sf = JME::JetResolutionScaleFactor(jerAK4PFPuppiSF_);
  }
  JME::JetParameters parameters;
  parameters.setJetPt(jet.pt());
  parameters.setJetEta(jet.eta());
  parameters.setRho(rhoJER);
  float relpterr = resolution.getResolution(parameters);
  if(genJetPt>0. && deltaR(jet.eta(),jet.phi(),jet.genJet()->eta(),jet.genJet()->phi())<0.2
     && (abs(jet.pt()-jet.genJet()->pt())<3*relpterr*jet.pt())) {
    JERScaleFactor     = (std::max(0., genJetPt + cFactorJER*diffPt))/recoJetPt;
    JERScaleFactorUP   = (std::max(0., genJetPt + cFactorJERup*diffPt))/recoJetPt;
    JERScaleFactorDOWN = (std::max(0., genJetPt + cFactorJERdown*diffPt))/recoJetPt;
  } else {
    JERScaleFactor     = 1.;
    JERScaleFactorUP   = 1.;
    JERScaleFactorDOWN = 1.;
  }
}

void TTHbb_eventSelector::FillJetVectors(const edm::Event& iEvent, std::vector<const reco::Candidate*> looseleps, std::vector<const pat::Jet*>& leading_jets, std::vector<const pat::Jet*>& subleading_jets, std::vector<const pat::Jet*>& leading_btags, std::vector<const pat::Jet*>& subleading_btags, vector<pair<double,int> >& leading_jet_csv_pos, vector<pair<double,int> >& subleading_jet_csv_pos){

  edm::Handle<reco::VertexCollection> vtx_h;
  iEvent.getByToken(vtx_h_, vtx_h);
  edm::Handle<double> rhopogHandle;
  iEvent.getByToken(rhopogHandle_,rhopogHandle);
  double rhopog = *rhopogHandle;
  edm::Handle<double> rhoJERHandle;
  iEvent.getByToken(rhoJERHandle_,rhoJERHandle);
  double rhoJER = *rhoJERHandle;
  edm::Handle<pat::JetCollection> jets;
  iEvent.getByToken(jets_, jets);

  // OR performed using 'subleading' jets definition.
  // Leading jets are subset of subleading jets.

  int jet_pos = 0; //This counter helps to order jets
  for(const pat::Jet &j : *jets){
    int vtxsize = vtx_h->size();

    //=== Good subleading jets ===
    if(!is_good_jet(j,rhopog,rhoJER,vtxsize,20)){jet_pos++; continue;}
    bool jetmatchedlepts = false;
    for(uint gl=0; gl<looseleps.size(); gl++) if(deltaR(looseleps[gl]->p4(),j.p4())<0.4) jetmatchedlepts = true;
    if(jetmatchedlepts){jet_pos++; continue;}
    double csvcurrjet = j.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
    subleading_jets.push_back((const pat::Jet*)&j);
    subleading_jet_csv_pos.push_back(make_pair(csvcurrjet,jet_pos));
    if(csvcurrjet>0.8484) {
      subleading_btags.push_back((const pat::Jet*)&j);
    }

    //=== Good leading jets ===
    if(!is_good_jet(j,rhopog,rhoJER,vtxsize,30)){jet_pos++; continue;}
    leading_jets.push_back((const pat::Jet*)&j);
    leading_jet_csv_pos.push_back(make_pair(csvcurrjet,jet_pos));
    if(csvcurrjet>0.8484) {
      leading_btags.push_back((const pat::Jet*)&j);
    }
    jet_pos++;
  }

  return;
}

void TTHbb_eventSelector::FillLeptonVectors(const edm::Event& iEvent, std::vector<const reco::Candidate*>& looseleps, std::vector<const reco::Candidate*>& loose_electrons, std::vector<const reco::Candidate*>& loose_muons, std::vector<const reco::Candidate*>& tightleps, std::vector<const reco::Candidate*>& tight_electrons, std::vector<const reco::Candidate*>& tight_muons){
  edm::Handle<reco::VertexCollection> vtx_h;
  iEvent.getByToken(vtx_h_, vtx_h);
  edm::Handle<double> rhopogHandle;
  iEvent.getByToken(rhopogHandle_,rhopogHandle);
  double rhopog = *rhopogHandle;
  edm::Handle<edm::ValueMap<bool>  > mvanontrig_id_decisions;
  iEvent.getByToken(eleMVAnonTrigIdMap_,mvanontrig_id_decisions);

  edm::Handle<edm::ValueMap<bool>  > medium_id_decisions;
  iEvent.getByToken(electronMediumIdMapToken_,medium_id_decisions);

  edm::Handle<edm::View<pat::Electron> > electron_pat;
  iEvent.getByToken(electron_pat_, electron_pat);
  edm::Handle<edm::View<pat::Muon> > muon_h;
  iEvent.getByToken(muon_h_, muon_h);
  const reco::Vertex &PV = vtx_h->front();

  //=== Muons ===
  //is_XXXX_muon selections have no trigger requirements atm.
  for(const pat::Muon &mu : *muon_h){
    if(!is_loose_muon(mu,PV)) continue;
    looseleps.push_back((const reco::Candidate*)&mu);
    loose_muons.push_back((const reco::Candidate*)&mu);
    if(!is_tight_muon(mu,PV)) continue;
    tightleps.push_back((const reco::Candidate*)&mu);
    tight_muons.push_back((const reco::Candidate*)&mu);
  }

  //=== Electrons ===
  for(edm::View<pat::Electron>::const_iterator ele = electron_pat->begin(); ele != electron_pat->end(); ele++){
    const Ptr<pat::Electron> elPtr(electron_pat, ele - electron_pat->begin() );
    if(!(is_loose_electron(*ele,rhopog))) continue;
    //bool isPassMvanontrig = (*mvanontrig_id_decisions)[ elPtr ];
    //if(!(isPassMvanontrig)) continue;
    bool isPassEletrig = (*medium_id_decisions)[ elPtr ];
    if(!(isPassEletrig)) continue;
    if(!(rel_iso_dbc_ele(*ele,rhopog)<0.15)) continue;
    const pat::Electron &lele = *ele;
    looseleps.push_back((const reco::Candidate*)&lele);
    loose_electrons.push_back((const reco::Candidate*)&lele);
    if(!(is_tight_electron(*ele,rhopog))) continue;
    tightleps.push_back((const reco::Candidate*)&lele);
    tight_electrons.push_back((const reco::Candidate*)&lele);
  }
  return;
}

int TTHbb_eventSelector::dilepton_channel(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  // Data access requests.
  // Using 'consumes' for these objects in constructor which returns edm::EDGetToken so use '.getByToken'.

  is_DL = 0;
  edm::Handle<reco::VertexCollection> vtx_h;
  iEvent.getByToken(vtx_h_, vtx_h);

  // Access electron candidate information.
  edm::Handle<edm::ValueMap<bool>  > veto_id_decisions;
  edm::Handle<edm::ValueMap<bool>  > loose_id_decisions;
  edm::Handle<edm::ValueMap<bool>  > medium_id_decisions;
  edm::Handle<edm::ValueMap<bool>  > tight_id_decisions;
  edm::Handle<edm::ValueMap<bool>  > heep_id_decisions;
  edm::Handle<edm::ValueMap<bool>  > mvatrig_id_decisions;
  edm::Handle<edm::ValueMap<bool>  > mvanontrig_id_decisions;
  edm::Handle<edm::ValueMap<bool>  > mvatrigwp90_id_decisions;
  edm::Handle<edm::ValueMap<bool>  > mvanontrigwp90_id_decisions;
  edm::Handle<edm::ValueMap<float> > elemvaValues_nonTrig;
  edm::Handle<edm::ValueMap<int> >   elemvaCategories_nonTrig;
  edm::Handle<edm::ValueMap<float> > elemvaValues_Trig;
  edm::Handle<edm::ValueMap<int> >   elemvaCategories_Trig;
  iEvent.getByToken(electronVetoIdMapToken_,   veto_id_decisions);
  iEvent.getByToken(electronLooseIdMapToken_,  loose_id_decisions);
  iEvent.getByToken(electronMediumIdMapToken_, medium_id_decisions);
  iEvent.getByToken(electronTightIdMapToken_,  tight_id_decisions);
  iEvent.getByToken(eleMVATrigIdMapToken_,     mvatrig_id_decisions);
  iEvent.getByToken(eleMVAnonTrigIdMap_,       mvanontrig_id_decisions);
  iEvent.getByToken(eleMVATrigwp90IdMap_,      mvatrigwp90_id_decisions);
  iEvent.getByToken(eleMVAnonTrigwp90IdMap_,   mvanontrigwp90_id_decisions);
  iEvent.getByToken(eleHEEPIdMapToken_,        heep_id_decisions);
  iEvent.getByToken(elemvaValuesMapToken_nonTrig_,     elemvaValues_nonTrig);
  iEvent.getByToken(elemvaCategoriesMapToken_nonTrig_, elemvaCategories_nonTrig);
  iEvent.getByToken(elemvaValuesMapToken_Trig_,        elemvaValues_Trig);
  iEvent.getByToken(elemvaCategoriesMapToken_Trig_,    elemvaCategories_Trig);

  // Access data from EventSetup Record (passed to method as argument):
  edm::ESHandle<TransientTrackBuilder> ttrkbuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",ttrkbuilder);

  //=======================
  // TTH OBJECT DEFINITION
  //=======================

  // === Primary Vertex ===
  // If no PV is found, skip event.
  if(vtx_h->empty()) return is_DL;
  reco::VertexCollection::const_iterator firstgoodVertex = vtx_h->end();
  for(reco::VertexCollection::const_iterator it = vtx_h->begin(); it != firstgoodVertex; it++){
    if(isGoodVertex(*it)){
      firstgoodVertex = it;
      break;
    }
  }
  if(firstgoodVertex == vtx_h->end()) return is_DL;
  //const reco::Vertex &PV = vtx_h->front();


  //=== Leptons ===
  std::vector<const reco::Candidate*> looseleps;
  std::vector<const reco::Candidate*> loose_electrons;
  std::vector<const reco::Candidate*> loose_muons;

  std::vector<const reco::Candidate*> tightleps;
  std::vector<const reco::Candidate*> tight_electrons;
  std::vector<const reco::Candidate*> tight_muons;

  FillLeptonVectors(iEvent, looseleps, loose_electrons, loose_muons, tightleps, tight_electrons, tight_muons);

  //=== Jets ===
  std::vector<const pat::Jet*> subleading_jets;
  std::vector<const pat::Jet*> leading_jets;
  std::vector<const pat::Jet*> subleading_btags;
  std::vector<const pat::Jet*> leading_btags;

  int jet_num = 0; // # Analysis jets
  int jetb_num = 0; // # Analysis b-tag.
  vector<pair<double,int> > leading_jet_csv_pos;
  vector<pair<double,int> > subleading_jet_csv_pos;

  FillJetVectors(iEvent, looseleps, leading_jets, subleading_jets, leading_btags, subleading_btags, leading_jet_csv_pos, subleading_jet_csv_pos);

  jet_num = leading_jets.size();
  jetb_num = subleading_btags.size();

  //====================================
  //            TTH SELECTION
  //====================================
  // Presuming all the above definitions
  // match those of the current ttH
  // analysis, you can now seperate
  // events into SL or DL channels.
  //====================================
  //=== DL ===
  //Ele pT for lead is from tightleps (pT>30 GeV)
  //Jet pT > 30 GeV (as in single lepton channel)

  //if(debug_) std::cout << ("TTHbb_eventSelector:dilepton") << ": DL selection" << std::endl;
  //if(debug_) std::cout << ("TTHbb_eventSelector:dilepton") << ": # tightleps: " << tightleps.size() << std::endl;
  //if(debug_) std::cout << ("TTHbb_eventSelector:dilepton") << ": # looseleps: " << looseleps.size() << std::endl;
  //if(debug_) std::cout << ("TTHbb_eventSelector:dilepton") << ": # jet_num: " << jet_num << std::endl;
  //if(debug_) std::cout << ("TTHbb_eventSelector:dilepton") << ": # jetb_num : " << jetb_num << std::endl;

  if(tightleps.size()>=1 && looseleps.size()==2 && jet_num>=2 && jetb_num>=1) {
      if (loose_muons.size()==1 && loose_electrons.size()==1){
        is_DL = 1;
      }
      if (loose_muons.size() == 2){
        is_DL = 2;
      }
      if (loose_electrons.size()==2){
        is_DL = 3;
      }
  }

  // Return types:
  //              0 = failed selection
  //              1 = el+mu
  //              2 = mu+mu
  //              3 = el+el
  return is_DL;

}

int TTHbb_eventSelector::singlelepton_channel(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  is_SL=0;

  edm::Handle<reco::VertexCollection> vtx_h;
  iEvent.getByToken(vtx_h_, vtx_h);

  // Access electron candidate information.
  edm::Handle<edm::ValueMap<bool>  > veto_id_decisions;
  edm::Handle<edm::ValueMap<bool>  > loose_id_decisions;
  edm::Handle<edm::ValueMap<bool>  > medium_id_decisions;
  edm::Handle<edm::ValueMap<bool>  > tight_id_decisions;
  edm::Handle<edm::ValueMap<bool>  > heep_id_decisions;
  edm::Handle<edm::ValueMap<bool>  > mvatrig_id_decisions;
  edm::Handle<edm::ValueMap<bool>  > mvanontrig_id_decisions;
  edm::Handle<edm::ValueMap<bool>  > mvatrigwp90_id_decisions;
  edm::Handle<edm::ValueMap<bool>  > mvanontrigwp90_id_decisions;
  edm::Handle<edm::ValueMap<float> > elemvaValues_nonTrig;
  edm::Handle<edm::ValueMap<int> >   elemvaCategories_nonTrig;
  edm::Handle<edm::ValueMap<float> > elemvaValues_Trig;
  edm::Handle<edm::ValueMap<int> >   elemvaCategories_Trig;
  iEvent.getByToken(electronVetoIdMapToken_,   veto_id_decisions);
  iEvent.getByToken(electronLooseIdMapToken_,  loose_id_decisions);
  iEvent.getByToken(electronMediumIdMapToken_, medium_id_decisions);
  iEvent.getByToken(electronTightIdMapToken_,  tight_id_decisions);
  iEvent.getByToken(eleMVATrigIdMapToken_,     mvatrig_id_decisions);
  iEvent.getByToken(eleMVAnonTrigIdMap_,       mvanontrig_id_decisions);
  iEvent.getByToken(eleMVATrigwp90IdMap_,      mvatrigwp90_id_decisions);
  iEvent.getByToken(eleMVAnonTrigwp90IdMap_,   mvanontrigwp90_id_decisions);
  iEvent.getByToken(eleHEEPIdMapToken_,        heep_id_decisions);
  iEvent.getByToken(elemvaValuesMapToken_nonTrig_,     elemvaValues_nonTrig);
  iEvent.getByToken(elemvaCategoriesMapToken_nonTrig_, elemvaCategories_nonTrig);
  iEvent.getByToken(elemvaValuesMapToken_Trig_,        elemvaValues_Trig);
  iEvent.getByToken(elemvaCategoriesMapToken_Trig_,    elemvaCategories_Trig);

  // Access data from EventSetup Record (passed to method as argument):
  edm::ESHandle<TransientTrackBuilder> ttrkbuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",ttrkbuilder);

  //=======================
  // TTH OBJECT DEFINITION
  //=======================

  // === Primary Vertex ===
  // If no PV is found, skip event.
  if(vtx_h->empty()) return is_SL;
  reco::VertexCollection::const_iterator firstgoodVertex = vtx_h->end();
  for(reco::VertexCollection::const_iterator it = vtx_h->begin(); it != firstgoodVertex; it++){
    if(isGoodVertex(*it)){
      firstgoodVertex = it;
      break;
    }
  }
  if(firstgoodVertex == vtx_h->end()) return is_SL;
  //const reco::Vertex &PV = vtx_h->front();


  //=== Leptons ===
  std::vector<const reco::Candidate*> looseleps;
  std::vector<const reco::Candidate*> loose_electrons;
  std::vector<const reco::Candidate*> loose_muons;

  std::vector<const reco::Candidate*> tightleps;
  std::vector<const reco::Candidate*> tight_electrons;
  std::vector<const reco::Candidate*> tight_muons;

  FillLeptonVectors(iEvent, looseleps, loose_electrons, loose_muons, tightleps, tight_electrons, tight_muons);

  //=== Jets ===
  std::vector<const pat::Jet*> subleading_jets;
  std::vector<const pat::Jet*> leading_jets;
  std::vector<const pat::Jet*> subleading_btags;
  std::vector<const pat::Jet*> leading_btags;

  int jet_num = 0; // # Analysis jets
  int jetb_num = 0; // # Analysis b-tag.
  vector<pair<double,int> > leading_jet_csv_pos;
  vector<pair<double,int> > subleading_jet_csv_pos;

  FillJetVectors(iEvent, looseleps, leading_jets, subleading_jets, leading_btags, subleading_btags, leading_jet_csv_pos, subleading_jet_csv_pos);

  jet_num = leading_jets.size();
  jetb_num = leading_btags.size();


  //====================================
  //            TTH SELECTION
  //====================================
  // Presuming all the above definitions
  // match those of the current ttH
  // analysis, you can now seperate
  // events into SL or DL channels.
  //====================================

  if(debug_) std::cout << ("TTHbb_eventSelector:singlelepton_channel") << ": SL selection" << std::endl;
  if(debug_) std::cout << ("TTHbb_eventSelector:singlelepton_channel") << ": tightleps.size(): " << tightleps.size() << std::endl;
  if(debug_) std::cout << ("TTHbb_eventSelector:singlelepton_channel") << ": looseleps.size() : " << looseleps.size() << std::endl;
  if(debug_) std::cout << ("TTHbb_eventSelector:singlelepton_channel") << ": jet_num: " << jet_num << std::endl;
  if(debug_) std::cout << ("TTHbb_eventSelector:singlelepton_channel") << ": jetb_num : " << jetb_num << std::endl;

  //=== SL ===
  if(tightleps.size()==1 && looseleps.size()==1 && jet_num>=4 && jetb_num>=2){
    if(tight_muons.size()==1 && loose_electrons.size()==0){
      is_SL = 1;
    }
    if(tight_electrons.size()==1 && loose_muons.size()==0){
      is_SL = 2;
    }
  }

  // Return types:
  //              0 = failed selection
  //              1 = mu
  //              2 = el
  return is_SL;

}

void TTHbb_eventSelector::Fill(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  tmp_evtNum     = iEvent.id().event();
  if(debug_) std::cout << ("TTHbb_eventSelector:Fill") << ": Running Fill method"<<endl;
  if(debug_) std::cout << ("TTHbb_eventSelector:Fill") << ": Event Number = " << tmp_evtNum << std::endl;

  is_singlelep = singlelepton_channel(iEvent,iSetup);
  is_dilepton = dilepton_channel(iEvent,iSetup);
  if(debug_) std::cout << ("TTHbb_eventSelector:Fill") << ": is_singlelep = " << is_singlelep << std::endl;
  if(debug_) std::cout << ("TTHbb_eventSelector:Fill") << ": is_dilepton = " << is_dilepton << std::endl;
  return;
}

void TTHbb_eventSelector::Clean(){

  //delete *jecAK4PFchsMC_;
  //delete *jecAK4PFchsMCUnc_;
  //delete *jecAK4PFchsDATA_;
  //delete *jecAK4PFchsDATAUnc_;

  return;
}
