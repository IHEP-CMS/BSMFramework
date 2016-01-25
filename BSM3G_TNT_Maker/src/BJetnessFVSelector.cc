#include "BSMFramework/BSM3G_TNT_Maker/interface/BJetnessFVSelector.h"
#include "BSMFramework/BSM3G_TNT_Maker/interface/ElectronPatSelector.h"
//Track builder infos
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexUpdator.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexTrackCompatibilityEstimator.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexSmoother.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/TrackExtra.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "RecoBTag/BTagTools/interface/SignedTransverseImpactParameter.h"
#include "TMath.h"
KalmanVertexFitter vtxFitterFV(true);
using namespace std;
using namespace pat;
using namespace edm;
using namespace reco;
BJetnessFVSelector::BJetnessFVSelector(std::string name, TTree* tree, bool debug, const pset& iConfig, edm::ConsumesCollector && ic):
  baseTree(name,tree,debug),
  eleMVATrigIdMapToken_(ic.consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleMVATrigIdMap")))
  {
  _muonToken        = iConfig.getParameter<edm::InputTag>("muons");
  _patElectronToken = iConfig.getParameter<edm::InputTag>("patElectrons");
  jetToken_         = iConfig.getParameter<edm::InputTag>("jets");
  _vertexInputTag   = iConfig.getParameter<edm::InputTag>("vertices");
  SetBranches();
}
BJetnessFVSelector::~BJetnessFVSelector(){
  delete tree_;
}
void BJetnessFVSelector::Fill(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  Clear();
  /////
  //   Recall collections
  ///// 
  edm::Handle<reco::VertexCollection> vtx_h;
  iEvent.getByLabel(_vertexInputTag, vtx_h);
  edm::Handle<edm::View<pat::Muon> > muon_h;
  iEvent.getByLabel(_muonToken, muon_h);
  edm::Handle<edm::View<pat::Electron> > electron_pat;
  iEvent.getByLabel(_patElectronToken, electron_pat);
  edm::Handle<edm::ValueMap<bool>  > mvatrig_id_decisions;
  iEvent.getByToken(eleMVATrigIdMapToken_, mvatrig_id_decisions);
  edm::Handle<pat::JetCollection> jets;                                       
  iEvent.getByLabel(jetToken_, jets);                                         
  edm::Handle<double> rhopogHandle;
  iEvent.getByLabel("fixedGridRhoFastjetAll",rhopogHandle);
  double rhopog = *rhopogHandle;
  edm::ESHandle<TransientTrackBuilder> ttrkbuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",ttrkbuilder);
  /////
  //   First clean the jet according the TTHbb selection
  /////
  //Require a good vertex 
  if(vtx_h->empty()) return; // skip the event if no PV found
  const reco::Vertex &PV = vtx_h->front();
  //Look for loose muons, electrons to clean jets
  vector<const reco::Candidate*> looseleps;
  //Muons
  for(const pat::Muon &mu : *muon_h){
    if(!is_loose_muon(mu,PV)) continue;
    looseleps.push_back((const reco::Candidate*)&mu);
  }
  //Electrons
  for(edm::View<pat::Electron>::const_iterator ele = electron_pat->begin(); ele != electron_pat->end(); ele++){
    const Ptr<pat::Electron> elPtr(electron_pat, ele - electron_pat->begin() );
    bool isPassMvatrig = (*mvatrig_id_decisions)[ elPtr ];
    if(!(is_loose_electron(*ele,rhopog) && isPassMvatrig)) continue;
    const pat::Electron &lele = *ele;
    looseleps.push_back((const reco::Candidate*)&lele);
  }
  //Get the good jets of the event
  //Iterate to access jet by decreasing b-tagging value
  int jet_pos = 0; //This counter helps to order jets
  int jet_num = 0; //This counter accounts for the number of good jets in the events
                   //The definition of good jet in the event must be the same of the TTHbb analysis
                   //so that jet_num corresponds to the number of jets that define the categories in the TTHbb search
  vector<pair<double,int> > jet_csv_pos;
  for(const pat::Jet &j : *jets){ 
    if(!is_good_jet(j)){jet_pos++; continue;}
    bool jetmatchedlepts = false;
    for(uint gl=0; gl<looseleps.size(); gl++) if(deltaR(looseleps[gl]->p4(),j.p4())<0.4) jetmatchedlepts = true;
    if(jetmatchedlepts){jet_pos++; continue;}
    double csvcurrjet = j.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
    jet_csv_pos.push_back(make_pair(csvcurrjet,jet_pos));
    jet_pos++;
    jet_num++;
  }
  sort(jet_csv_pos.rbegin(), jet_csv_pos.rend());//Order by descreasing csv value
  if(jet_num!=0){
    /////
    //   From here on it starts to define the BJetness variables
    /////
    //You need to provide as input the jets selected in the event (selection according to the TTHbb analysis),
    //excluding the jet with the highest CSV (start from jn=1 below)   
    vector<pat::Jet> evtjets; evtjets.clear();
    int maxjetnum = 6; //This value has been chosen after optimisation 
    if(jet_num<maxjetnum)  maxjetnum = jet_num;
    for(int jn=1; jn<maxjetnum; jn++) evtjets.push_back((*jets)[jet_csv_pos[jn].second]);
    //Define the variables you want to access 
    double bjetnessFV_num_loosenoipnoiso_leps = -1;
    double bjetnessFV_numjettrksnopv          = -1;
    double bjetnessFV_pvTrkOVcollTrk          = -1;
    double bjetnessFV_avip3d_val              = -1;
    double bjetnessFV_avip3d_sig              = -1;
    double bjetnessFV_avsip3d_sig             = -1;
    double bjetnessFV_avip1d_sig              = -1;   
    //This is the method to access the BJetness variables
    get_bjetness_vars(
                      //Inputs:
                      evtjets,      //Jets selected in the event according to the TTHbb selection 
                      PV,           //Prinary vertex of the event 
                      *ttrkbuilder, //Transient tracker builder to measure impact parameters
                      electron_pat, muon_h, //Leptons collections to count the number of electrons and muons 
                      //BJetness variables  
                      bjetnessFV_num_loosenoipnoiso_leps,bjetnessFV_numjettrksnopv,bjetnessFV_pvTrkOVcollTrk,bjetnessFV_avip3d_val,bjetnessFV_avip3d_sig,bjetnessFV_avsip3d_sig,bjetnessFV_avip1d_sig
                     );
    //Fill the quantities for the event
    //Num_of_trks
    BJetnessFV_num_loosenoipnoiso_leps.push_back(bjetnessFV_num_loosenoipnoiso_leps);
    BJetnessFV_numjettrksnopv.push_back(bjetnessFV_numjettrksnopv);
    BJetnessFV_pvTrkOVcollTrk.push_back(bjetnessFV_pvTrkOVcollTrk);
    //ImpactParameter  
    BJetnessFV_avip3d_val.push_back(bjetnessFV_avip3d_val);
    BJetnessFV_avip3d_sig.push_back(bjetnessFV_avip3d_sig);
    BJetnessFV_avsip3d_sig.push_back(bjetnessFV_avsip3d_sig);
    BJetnessFV_avip1d_sig.push_back(bjetnessFV_avip1d_sig); 
 }else{//if(jet_num!=0)
   //Num_of_trks
   BJetnessFV_num_loosenoipnoiso_leps.push_back(-999);
   BJetnessFV_numjettrksnopv.push_back(-999);
   BJetnessFV_pvTrkOVcollTrk.push_back(-999);
   //ImpactParameter  
   BJetnessFV_avip3d_val.push_back(-999);
   BJetnessFV_avip3d_sig.push_back(-999);
   BJetnessFV_avsip3d_sig.push_back(-999);
   BJetnessFV_avip1d_sig.push_back(-999);
 }
}
/////
//   Variables
/////
void BJetnessFVSelector::SetBranches(){
  if(debug_) std::cout<<"setting branches: calling AddBranch of baseTree"<<std::endl;
  //Num_of_trks
  AddBranch(&BJetnessFV_num_loosenoipnoiso_leps  ,"BJetnessFV_num_loosenoipnoiso_leps");
  AddBranch(&BJetnessFV_numjettrksnopv           ,"BJetnessFV_numjettrksnopv");
  AddBranch(&BJetnessFV_pvTrkOVcollTrk           ,"BJetnessFV_pvTrkOVcollTrk");
  //ImpactParameter
  AddBranch(&BJetnessFV_avip3d_val               ,"BJetnessFV_avip3d_val");
  AddBranch(&BJetnessFV_avip3d_sig               ,"BJetnessFV_avip3d_sig");
  AddBranch(&BJetnessFV_avsip3d_sig              ,"BJetnessFV_avsip3d_sig");
  AddBranch(&BJetnessFV_avip1d_sig               ,"BJetnessFV_avip1d_sig");
  if(debug_) std::cout<<"set branches"<<std::endl;
}
void BJetnessFVSelector::Clear(){
  //Num_of_trks
  BJetnessFV_num_loosenoipnoiso_leps.clear();
  BJetnessFV_numjettrksnopv.clear();
  BJetnessFV_pvTrkOVcollTrk.clear();
  //ImpactParameter
  BJetnessFV_avip3d_val.clear();
  BJetnessFV_avip3d_sig.clear();
  BJetnessFV_avsip3d_sig.clear();
  BJetnessFV_avip1d_sig.clear();
}
/////
//   Methods to be aligned to the TTHbb selection
/////
//Look for loose muon (definition for the jet cleaning)
bool BJetnessFVSelector::is_loose_muon(const pat::Muon& mu, const reco::Vertex& vtx){
  bool isloosemu = false;
  if(mu.pt()>15 &&
    TMath::Abs(mu.eta()) < 2.4 &&
    mu.isTightMuon(vtx) &&  
    rel_iso_dbc_mu(mu) < 0.25
    ) isloosemu = true;
  return isloosemu;
}
double BJetnessFVSelector::rel_iso_dbc_mu(const pat::Muon& lepton){
  return((lepton.pfIsolationR04().sumChargedHadronPt + max(lepton.pfIsolationR04().sumNeutralHadronEt + lepton.pfIsolationR04().sumPhotonEt - 0.5 * lepton.pfIsolationR04().sumPUPt,0.0))/lepton.pt()
        );
}
//Look for loose electron (definition for the jet cleaning)
bool BJetnessFVSelector::is_loose_electron(const pat::Electron& ele, double rhopog){
  bool isele = false;
  if(ele.pt()>15 && TMath::Abs(ele.eta())<2.4 &&
     !(fabs(ele.superCluster()->position().eta()) > 1.4442 && fabs(ele.superCluster()->position().eta()) < 1.5660)){
    if(fabs(ele.superCluster()->position().eta())<1.4442 
      && (ele.full5x5_sigmaIetaIeta()<0.012)                 
      && (ele.hcalOverEcal()<0.09)                           
      && (ele.ecalPFClusterIso()/ele.pt()<0.37)            
      && (ele.hcalPFClusterIso()/ele.pt()<0.25)            
      && (ele.dr03TkSumPt()/ele.pt()<0.18)                 
      && (fabs(ele.deltaEtaSuperClusterTrackAtVtx())<0.0095)
      && (fabs(ele.deltaPhiSuperClusterTrackAtVtx())<0.065)){  
      isele = true;
    }
    if(fabs(ele.superCluster()->position().eta())>1.5660
      && ele.full5x5_sigmaIetaIeta()<0.033
      && ele.hcalOverEcal()<0.09
      && (ele.ecalPFClusterIso()/ele.pt())<0.45
      && (ele.hcalPFClusterIso()/ele.pt())<0.28
      && (ele.dr03TkSumPt()/ele.pt())<0.18){
      isele = true;
    }
  }
  //ele.passConversionVeto()
  return isele; 
}
double BJetnessFVSelector::rel_iso_dbc_ele(const pat::Electron& lepton, double rhopog){
  double effarea          = get_effarea(lepton.superCluster()->position().eta());  
  double SumChHadPt       = lepton.pfIsolationVariables().sumChargedHadronPt;
  double SumNeuHadEt      = lepton.pfIsolationVariables().sumNeutralHadronEt;
  double SumPhotonEt      = lepton.pfIsolationVariables().sumPhotonEt;
  double SumNeutralCorrEt = std::max( 0.0, SumNeuHadEt+SumPhotonEt - rhopog*effarea );
  return ((SumChHadPt + SumNeutralCorrEt)/lepton.pt());
}
double BJetnessFVSelector::get_effarea(double eta){
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
//Require good jets (according to TTHbb analysis)
bool BJetnessFVSelector::is_good_jet(const pat::Jet &j){
  bool isgoodjet = true;
  //Acceptance
  if(j.pt() < 30)       isgoodjet = false; //Please note that this requirement is for the SL channel, while for DL channel we require pT > 20! 
  if(fabs(j.eta())>2.4) isgoodjet = false; 
  //ID requirements
  if(j.neutralHadronEnergyFraction() >= 0.99) isgoodjet = false;
  if(j.chargedEmEnergyFraction()     >= 0.99) isgoodjet = false;
  if(j.neutralEmEnergyFraction()     >= 0.99) isgoodjet = false;
  if(j.numberOfDaughters()           <= 1)    isgoodjet = false;
  if(j.chargedHadronEnergyFraction() <= 0.0)  isgoodjet = false;
  if(j.chargedMultiplicity()         <= 0.0)  isgoodjet = false;
  return isgoodjet;
}
/////
//   Methods for the BJetness variables
/////
void BJetnessFVSelector::get_bjetness_vars(
                                           vector<pat::Jet> evtjets, const reco::Vertex& vtx, const TransientTrackBuilder& ttrkbuilder, edm::Handle<edm::View<pat::Electron> > electron_pat, edm::Handle<edm::View<pat::Muon> > muon_h,
                                           double& bjetnessFV_num_loosenoipnoiso_leps, double& bjetnessFV_numjettrksnopv, double& bjetnessFV_pvTrkOVcollTrk, double& bjetnessFV_avip3d_val, double& bjetnessFV_avip3d_sig, double& bjetnessFV_avsip3d_sig, double& bjetnessFV_avip1d_sig
                                          ){
  //Get BJetness trk info
  vector<Track> jetschtrks; jetschtrks.clear(); 
  double num_pvtrks              = 0;
  double num_npvtrks             = 0;
  double num_loosenoipnoiso_eles = 0;  
  double num_loose_mus           = 0;           
  vector<tuple<double, double, double> > jetsdir; jetsdir.clear(); 
  get_bjetness_trkinfos(evtjets, vtx, jetschtrks, num_pvtrks, num_npvtrks, electron_pat, muon_h, num_loosenoipnoiso_eles, num_loose_mus, jetsdir);
  bjetnessFV_num_loosenoipnoiso_leps = num_loosenoipnoiso_eles+num_loose_mus;
  bjetnessFV_numjettrksnopv          = num_npvtrks;
  if(jetschtrks.size()!=0){
    bjetnessFV_pvTrkOVcollTrk        = num_pvtrks/double(jetschtrks.size()); 
    //Get BJetness Impact Parameters
    double ip_valtemp = 0;
    //3D
    double jetchtrks_avip3d_val  = 0;
    double jetchtrks_avip3d_sig  = 0;
    double jetchtrks_avsip3d_sig = 0;
    get_avip3d(jetschtrks, ttrkbuilder, vtx, jetsdir, jetchtrks_avip3d_val,jetchtrks_avip3d_sig,jetchtrks_avsip3d_sig);
    ip_valtemp = jetchtrks_avip3d_val/jetschtrks.size();
    if(ip_valtemp==ip_valtemp) bjetnessFV_avip3d_val = ip_valtemp;
    else                       bjetnessFV_avip3d_val = -996;
    ip_valtemp = jetchtrks_avip3d_sig/jetschtrks.size();
    if(ip_valtemp==ip_valtemp) bjetnessFV_avip3d_sig = ip_valtemp;
    else                       bjetnessFV_avip3d_sig = -996; 
    ip_valtemp = jetchtrks_avsip3d_sig/jetschtrks.size();
    if(ip_valtemp==ip_valtemp) bjetnessFV_avsip3d_sig = ip_valtemp;
    else                       bjetnessFV_avsip3d_sig = -996;
    //1D
    double jetchtrks_avip1d_sig  = 0;
    get_avip1d(jetschtrks, ttrkbuilder, vtx, jetsdir, jetchtrks_avip1d_sig);
    ip_valtemp = jetchtrks_avip1d_sig/jetschtrks.size();
    if(ip_valtemp==ip_valtemp) bjetnessFV_avip1d_sig = ip_valtemp;
    else                       bjetnessFV_avip1d_sig = -996;    
  }else{
    bjetnessFV_pvTrkOVcollTrk        = -998;
    bjetnessFV_avip3d_val            = -998;
    bjetnessFV_avip3d_sig            = -998;
    bjetnessFV_avsip3d_sig           = -998;
    bjetnessFV_avip1d_sig            = -998;
  }
}
//Get the BJetness trk info 
void BJetnessFVSelector::get_bjetness_trkinfos(vector<pat::Jet> evtjets, const reco::Vertex& vtx, vector<Track>& jetchtrks, double& bjetness_num_pvtrks, double& bjetness_num_npvtrks, edm::Handle<edm::View<pat::Electron> > electron_pat, edm::Handle<edm::View<pat::Muon> > muon_h, double& bjetness_num_loosenoipnoiso_eles, double& bjetness_num_loose_mus, vector<tuple<double, double, double> >& jetsdir){
  //Loop over evt jet
  for(uint j=0; j<evtjets.size(); j++){
    pat::Jet jet = evtjets[j];
    //Access jet daughters
    vector<CandidatePtr> jdaus(jet.daughterPtrVector());
    sort(jdaus.begin(), jdaus.end(), [](const reco::CandidatePtr &p1, const reco::CandidatePtr &p2) {return p1->pt() > p2->pt();});
    for(uint jd=0; jd<jdaus.size(); jd++){
      const pat::PackedCandidate &jcand = dynamic_cast<const pat::PackedCandidate &>(*jdaus[jd]);
      //dR requirement
      if(deltaR(jcand.p4(),jet.p4())>0.4) continue;
      Track trk = Track(jcand.pseudoTrack());
      bool isgoodtrk = is_goodtrk(trk,vtx);
      //Minimal conditions for a BJetness jet constituent 
      if(isgoodtrk && jcand.charge()!=0 && jcand.fromPV()>1){
        jetchtrks.push_back(trk);
        if(jcand.fromPV()==3) bjetness_num_pvtrks++;
        if(jcand.fromPV()==2) bjetness_num_npvtrks++;
        jetsdir.push_back(make_tuple(jet.px(),jet.py(),jet.pz()));
        if(fabs(jcand.pdgId())==13 && is_loosePOG_jetmuon(jcand,muon_h)) bjetness_num_loose_mus++;
        if(fabs(jcand.pdgId())==11 && is_loosePOGNoIPNoIso_jetelectron(jcand,electron_pat,vtx)) bjetness_num_loosenoipnoiso_eles++;
      }//Ch trks 
    }//Loop on jet daus 
  }//Loop on evt jet
}
//Check that the track is a good track
bool BJetnessFVSelector::is_goodtrk(Track trk,const reco::Vertex& vtx){
 bool isgoodtrk = false;
 if(trk.pt()>1 &&
   trk.hitPattern().numberOfValidHits()>=8 &&
   trk.hitPattern().numberOfValidPixelHits()>=2 &&
   trk.normalizedChi2()<5 &&
   std::abs(trk.dxy(vtx.position()))<0.2 &&
   std::abs(trk.dz(vtx.position()))<17
   ) isgoodtrk = true;
 return isgoodtrk;
}
//Look for loose muon (definition to look for candidates among jet daughters)
bool BJetnessFVSelector::is_loosePOG_jetmuon(const pat::PackedCandidate &jcand, edm::Handle<edm::View<pat::Muon> > muon_h){
  bool ismu = false;
  for(const pat::Muon &mu : *muon_h){
    if(deltaR(jcand.p4(),mu.p4())<0.1 && fabs(jcand.pt()-mu.pt())/mu.pt()<0.05){
     if(mu.isLooseMuon()) ismu = true;
     if(ismu) break;
    }
  }  
  return ismu;
}
//Look for loose electron ((definition to look for candidates among jet daughters)
bool BJetnessFVSelector::is_loosePOGNoIPNoIso_jetelectron(const pat::PackedCandidate &jcand, edm::Handle<edm::View<pat::Electron> > electron_pat, const reco::Vertex& vtx){
  bool isele = false;
  for(edm::View<pat::Electron>::const_iterator ele = electron_pat->begin(); ele != electron_pat->end(); ele++){
    const pat::Electron &lele = *ele; 
    if(deltaR(jcand.p4(),lele.p4())<0.1 && fabs(jcand.pt()-lele.pt())/lele.pt()<0.05){
      double ooEmooP = 999;
      if(lele.ecalEnergy()==0)                   ooEmooP = 1e30;
      else if(!std::isfinite(lele.ecalEnergy())) ooEmooP = 1e30;
      else                                       ooEmooP = fabs(1.0/lele.ecalEnergy() - lele.eSuperClusterOverP()/lele.ecalEnergy());
      if(!(fabs(lele.superCluster()->position().eta()) > 1.4442 && fabs(lele.superCluster()->position().eta()) < 1.5660)){
        if(fabs(lele.superCluster()->position().eta())<1.4442
          && (lele.full5x5_sigmaIetaIeta()<0.0103)
          && (lele.hcalOverEcal()<0.104)
          && (fabs(lele.deltaEtaSuperClusterTrackAtVtx())<0.0105)
          && (fabs(lele.deltaPhiSuperClusterTrackAtVtx())<0.115)
          && ooEmooP<0.102
          && lele.gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS)<=2
          && lele.passConversionVeto()   
          ){
            isele = true;
        }
        if(fabs(lele.superCluster()->position().eta())>1.5660
          && (lele.full5x5_sigmaIetaIeta()<0.0301)
          && (lele.hcalOverEcal()<0.0897)
          && (fabs(lele.deltaEtaSuperClusterTrackAtVtx())<0.00814)
          && (fabs(lele.deltaPhiSuperClusterTrackAtVtx())<0.182)
          && ooEmooP<0.126
          && lele.gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS)<=1
          && lele.passConversionVeto()   
          ){
            isele = true;
        }
      } 
      if(isele) break;
    } 
  } 
  return isele; 
}
//Methods related to IP
void BJetnessFVSelector::get_avip3d(vector<Track> trks, const TransientTrackBuilder& ttrkbuilder, reco::Vertex vtx, vector<tuple<double, double, double> >& jetsdir, double& jetchtrks_avip3d_val, double& jetchtrks_avip3d_sig, double& jetchtrks_avsip3d_sig){
  double valtemp = 0;
  vector<TransientTrack> ttrks = get_ttrks(trks,ttrkbuilder);
  for(uint t=0; t<ttrks.size(); t++){
    valtemp = IPTools::absoluteImpactParameter3D(ttrks[t],vtx).second.value();
    if(valtemp==valtemp) jetchtrks_avip3d_val  += valtemp;
    valtemp = IPTools::absoluteImpactParameter3D(ttrks[t],vtx).second.significance();
    if(valtemp==valtemp) jetchtrks_avip3d_sig  += valtemp;
    GlobalVector jetsdirgv(get<0>(jetsdir[t]),get<1>(jetsdir[t]),get<2>(jetsdir[t]));
    valtemp = IPTools::signedImpactParameter3D(ttrks[t],jetsdirgv,vtx).second.significance();
    if(valtemp==valtemp) jetchtrks_avsip3d_sig += valtemp;
  }
}
void BJetnessFVSelector::get_avip1d(vector<Track> trks, const TransientTrackBuilder& ttrkbuilder, reco::Vertex vtx, vector<tuple<double, double, double> >& jetsdir, double& jetchtrks_avip1d_sig){
  double valtemp = 0;
  vector<TransientTrack> ttrks = get_ttrks(trks,ttrkbuilder);
  SignedTransverseImpactParameter stip;
  for(uint t=0; t<ttrks.size(); t++){
    GlobalVector jetsdirgv(get<0>(jetsdir[t]),get<1>(jetsdir[t]),get<2>(jetsdir[t]));
    valtemp = fabs(stip.zImpactParameter(ttrks[t],jetsdirgv,vtx).second.significance());
    if(valtemp==valtemp) jetchtrks_avip1d_sig  += valtemp;
  }
}
vector<TransientTrack> BJetnessFVSelector::get_ttrks(vector<Track> trks, const TransientTrackBuilder& ttrkbuilder){
 vector<TransientTrack> ttrks;
 for(uint tr=0; tr<trks.size(); tr++){
  TransientTrack ttrk = ttrkbuilder.build(&trks[tr]);
  ttrks.push_back(ttrk);
 }
 return ttrks;
}
