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
  eleMVATrigIdMapToken_(ic.consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleMVATrigIdMap"))),
  eleMVAnonTrigIdMapToken_(ic.consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleMVAnonTrigIdMap")))
{
  vtx_h_              = ic.consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"));
  electron_pat_       = ic.consumes<edm::View<pat::Electron> >(iConfig.getParameter<edm::InputTag>("patElectrons"));
  muon_h_             = ic.consumes<edm::View<pat::Muon> >(iConfig.getParameter<edm::InputTag>("muons"));
  jets_               = ic.consumes<pat::JetCollection >(iConfig.getParameter<edm::InputTag>("jets"));
  rhopogHandle_       = ic.consumes<double>(edm::InputTag("fixedGridRhoFastjetAll"));
  rhoJERHandle_       = ic.consumes<double>(edm::InputTag("fixedGridRhoAll"));
  _vtx_ndof_min       = iConfig.getParameter<int>("vtx_ndof_min");
  _vtx_rho_max        = iConfig.getParameter<int>("vtx_rho_max");
  _vtx_position_z_max = iConfig.getParameter<double>("vtx_position_z_max");
  _is_data = iConfig.getParameter<bool>("is_data");
  jecPayloadNamesAK4PFchsMC1_   = iConfig.getParameter<edm::FileInPath>("jecPayloadNamesAK4PFchsMC1");
  jecPayloadNamesAK4PFchsMC2_   = iConfig.getParameter<edm::FileInPath>("jecPayloadNamesAK4PFchsMC2");
  jecPayloadNamesAK4PFchsMC3_   = iConfig.getParameter<edm::FileInPath>("jecPayloadNamesAK4PFchsMC3");
  jecPayloadNamesAK4PFchsMCUnc_ = iConfig.getParameter<edm::FileInPath>("jecPayloadNamesAK4PFchsMCUnc");
  jecPayloadNamesAK4PFchsDATA1_   = iConfig.getParameter<edm::FileInPath>("jecPayloadNamesAK4PFchsDATA1");
  jecPayloadNamesAK4PFchsDATA2_   = iConfig.getParameter<edm::FileInPath>("jecPayloadNamesAK4PFchsDATA2");
  jecPayloadNamesAK4PFchsDATA3_   = iConfig.getParameter<edm::FileInPath>("jecPayloadNamesAK4PFchsDATA3");
  jecPayloadNamesAK4PFchsDATA4_   = iConfig.getParameter<edm::FileInPath>("jecPayloadNamesAK4PFchsDATA4");
  jecPayloadNamesAK4PFchsDATAUnc_ = iConfig.getParameter<edm::FileInPath>("jecPayloadNamesAK4PFchsDATAUnc");
  jerAK4PFchs_     = iConfig.getParameter<edm::FileInPath>("jerAK4PFchs").fullPath();
  jerAK4PFchsSF_   = iConfig.getParameter<edm::FileInPath>("jerAK4PFchsSF").fullPath();
  JECInitialization();
  SetBranches();
  evSel = new TTHbb_eventSelector("miniAOD", tree, debug, iConfig, std::move(ic));
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
  iEvent.getByToken(vtx_h_, vtx_h);
  edm::Handle<edm::View<pat::Electron> > electron_pat;
  iEvent.getByToken(electron_pat_, electron_pat);
  edm::Handle<edm::View<pat::Muon> > muon_h;
  iEvent.getByToken(muon_h_, muon_h);
  edm::Handle<pat::JetCollection> jets;
  iEvent.getByToken(jets_, jets);
  edm::Handle<double> rhopogHandle;
  iEvent.getByToken(rhopogHandle_,rhopogHandle);
  double rhopog = *rhopogHandle;
  edm::Handle<double> rhoJERHandle;
  iEvent.getByToken(rhoJERHandle_,rhoJERHandle);
  double rhoJER = *rhoJERHandle;
  edm::Handle<edm::ValueMap<bool>  > mvatrig_id_decisions;
  iEvent.getByToken(eleMVATrigIdMapToken_, mvatrig_id_decisions);
  edm::Handle<edm::ValueMap<bool>  > mvanontrig_id_decisions;
  iEvent.getByToken(eleMVAnonTrigIdMapToken_, mvanontrig_id_decisions);
  edm::ESHandle<TransientTrackBuilder> ttrkbuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",ttrkbuilder);

  bool isSL_ = evSel->singlelepton_channel(iEvent,iSetup);
  bool isDL_ = evSel->dilepton_channel(iEvent,iSetup);



  /////
  //   First clean the jet according the TTHbb selection
  /////
  //Require a good vertex (This part has to be clarified (do we want the first PV to be a good one?))
  if(vtx_h->empty()) return; // skip the event if no PV found
  reco::VertexCollection::const_iterator firstgoodVertex = vtx_h->end();
  for(reco::VertexCollection::const_iterator it = vtx_h->begin(); it != firstgoodVertex; it++){
    if(isGoodVertex(*it)){
      firstgoodVertex = it;
      break;
    }
  }
  if(firstgoodVertex == vtx_h->end()) return;
  const reco::Vertex &PV = vtx_h->front(); //It still takes the first vertex



  //Look for muons, electrons
  vector<const reco::Candidate*> looseleps;
  vector<const reco::Candidate*> tightleps;


  //Muons
  for(const pat::Muon &mu : *muon_h){
    if(!is_loose_muon(mu,PV)) continue;
    looseleps.push_back((const reco::Candidate*)&mu);
    if(!is_tight_muon(mu,PV)) continue;
    tightleps.push_back((const reco::Candidate*)&mu);
    //cout<<setw(20)<<"MuonFV pt,eta,phi"<<setw(20)<<mu.pt()<<setw(20)<<mu.eta()<<setw(20)<<mu.phi()<<endl;
  }


  //Electrons
  for(edm::View<pat::Electron>::const_iterator ele = electron_pat->begin(); ele != electron_pat->end(); ele++){
    const Ptr<pat::Electron> elPtr(electron_pat, ele - electron_pat->begin() );
    //bool isPassMvatrig = (*mvatrig_id_decisions)[ elPtr ];
    //if(!(is_loose_electron(*ele,rhopog) && isPassMvatrig)) continue;
    bool isPassMvanontrig = (*mvanontrig_id_decisions)[ elPtr ];
    if(!(isPassMvanontrig)) continue;
    if(!(rel_iso_dbc_ele(*ele,rhopog)<0.15)) continue;
    const pat::Electron &lele = *ele;
    looseleps.push_back((const reco::Candidate*)&lele);
    if(!(is_tight_electron(*ele,rhopog) && isPassMvanontrig)) continue;
    tightleps.push_back((const reco::Candidate*)&lele);
    //cout<<setw(20)<<"ElectronFV pt,eta,phi"<<setw(20)<<ele->pt()<<setw(20)<<ele->eta()<<setw(20)<<ele->phi()<<endl;
  }



  //Get the good jets of the event
  //Iterate to access jet by decreasing b-tagging value
  int jet_pos = 0; //This counter helps to order jets
  int jet_num = 0; //This counter accounts for the number of good jets in the events
                   //The definition of good jet in the event must be the same of the TTHbb analysis
                   //so that jet_num corresponds to the number of jets that define the categories in the TTHbb search
  int jetb_num = 0;
  vector<pair<double,int> > jet_csv_pos;
  vector<pair<double,int> > jet_cmva_pos;
  for(const pat::Jet &j : *jets){
    int vtxsize = vtx_h->size();
    if(!is_good_jet(j,rhopog,rhoJER,vtxsize)){jet_pos++; continue;}
    bool jetmatchedlepts = false;
    for(uint gl=0; gl<looseleps.size(); gl++) if(deltaR(looseleps[gl]->p4(),j.p4())<0.4) jetmatchedlepts = true;
    if(jetmatchedlepts){jet_pos++; continue;}
    double csvcurrjet = j.bDiscriminator("newpfCombinedInclusiveSecondaryVertexV2BJetTags");
    jet_csv_pos.push_back(make_pair(csvcurrjet,jet_pos));
    BJetnessFV_jetcsv.push_back(csvcurrjet);
    double cmvacurrjet = j.bDiscriminator("newpfCombinedMVAV2BJetTags");
    jet_cmva_pos.push_back(make_pair(cmvacurrjet,jet_pos));
    BJetnessFV_pfCombinedMVAV2BJetTags.push_back(cmvacurrjet);
    double jetprobjet = j.bDiscriminator("newpfJetProbabilityBJetTags");
    BJetnessFV_pfJetProbabilityBJetTags.push_back(jetprobjet);
    if(csvcurrjet>0.8) jetb_num++;
    jet_pos++;
    jet_num++;
  }





  /////
  //   Select only TTHbb events (mainly lep sel)
  /////
  if(tightleps.size()==1 && looseleps.size()==1 && jet_num>=4 && jetb_num>=2) BJetnessFV_isSingleLepton = 1;
  //This selection is not really the dilepton one!
  //Ele pT for lead is from tightleps (pT>30 GeV)
  //Jet pT > 30 GeV (as in single lepton channel)
  //if(tightleps.size()==1 && looseleps.size()==2 && jet_num>=3 && jetb_num>=2) BJetnessFV_isDoubleLepton = 1;
  if(!(BJetnessFV_isSingleLepton==1 || BJetnessFV_isDoubleLepton==1)) return;






  /////
  //   You need to provide as input the jets selected in the event (selection according to the TTHbb analysis),
  //   which have to be ordered by decreasing b-tagging value
  /////
  sort(jet_csv_pos.rbegin(), jet_csv_pos.rend());//Order by descreasing csv value
  if(jet_num!=0){
    //cout<<"Num of jet is"<<setw(20)<<jet_num<<" "<<jetb_num<<endl;
    /////
    //   From here on it starts to define the BJetness variables
    //   Note that we exclude the jet with the highest CSV (start from jn=1 below) and consider up to maximum 6 jets in an event
    /////
    vector<pat::Jet> evtjets; evtjets.clear();
    int maxjetnum = 6; //This value has been chosen after optimisation
    if(jet_num<maxjetnum) maxjetnum = jet_num;
    for(int jn=1; jn<maxjetnum; jn++) evtjets.push_back((*jets)[jet_csv_pos[jn].second]);
    //Define the variables you want to access
    double bjetnessFV_num_leps        = -1;
    double bjetnessFV_npvTrkOVcollTrk = -1;
    double bjetnessFV_avip3d_val      = -1;
    double bjetnessFV_avip3d_sig      = -1;
    double bjetnessFV_avsip3d_sig     = -1;
    double bjetnessFV_avip1d_sig      = -1;
    //This is the method to access the BJetness variables
    get_bjetness_vars(
                      //Inputs:
                      evtjets,      //Jets used to build the BJetness jets
                      PV,           //Prinary vertex of the event
                      *ttrkbuilder, //Transient tracker builder to measure impact parameters
                      electron_pat, muon_h, //Leptons collections to count the number of electrons and muons
                      //BJetness variables
                      bjetnessFV_num_leps,bjetnessFV_npvTrkOVcollTrk,bjetnessFV_avip3d_val,bjetnessFV_avip3d_sig,bjetnessFV_avsip3d_sig,bjetnessFV_avip1d_sig
                     );
    //Fill the quantities for the event
    //Num_of_trks
    BJetnessFV_num_leps.push_back(bjetnessFV_num_leps);
    BJetnessFV_npvTrkOVcollTrk.push_back(bjetnessFV_npvTrkOVcollTrk);
    //ImpactParameter
    BJetnessFV_avip3d_val.push_back(bjetnessFV_avip3d_val);
    BJetnessFV_avip3d_sig.push_back(bjetnessFV_avip3d_sig);
    BJetnessFV_avsip3d_sig.push_back(bjetnessFV_avsip3d_sig);
    BJetnessFV_avip1d_sig.push_back(bjetnessFV_avip1d_sig);
  }else{//if(jet_num!=0)
    //Num_of_trks
    BJetnessFV_num_leps.push_back(-999);
    BJetnessFV_npvTrkOVcollTrk.push_back(-999);
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
  //Evt selection
  AddBranch(&BJetnessFV_isSingleLepton               ,"BJetnessFV_isSingleLepton");
  AddBranch(&BJetnessFV_isDoubleLepton               ,"BJetnessFV_isDoubleLepton");
  //BTag discriminators
  AddBranch(&BJetnessFV_jetcsv                   ,"BJetnessFV_jetcsv");
  AddBranch(&BJetnessFV_pfJetProbabilityBJetTags ,"BJetnessFV_pfJetProbabilityBJetTags");
  AddBranch(&BJetnessFV_pfCombinedMVAV2BJetTags  ,"BJetnessFV_pfCombinedMVAV2BJetTags");
  //Num_of_trks
  AddBranch(&BJetnessFV_num_leps  ,"BJetnessFV_num_leps");
  AddBranch(&BJetnessFV_npvTrkOVcollTrk           ,"BJetnessFV_npvTrkOVcollTrk");
  //ImpactParameter
  AddBranch(&BJetnessFV_avip3d_val               ,"BJetnessFV_avip3d_val");
  AddBranch(&BJetnessFV_avip3d_sig               ,"BJetnessFV_avip3d_sig");
  AddBranch(&BJetnessFV_avsip3d_sig              ,"BJetnessFV_avsip3d_sig");
  AddBranch(&BJetnessFV_avip1d_sig               ,"BJetnessFV_avip1d_sig");
  if(debug_) std::cout<<"set branches"<<std::endl;
}
void BJetnessFVSelector::Clear(){
  //Evt selection
  BJetnessFV_isSingleLepton = 0;
  BJetnessFV_isDoubleLepton = 0;
  //BTag discriminators
  BJetnessFV_jetcsv.clear();
  BJetnessFV_pfJetProbabilityBJetTags.clear();
  BJetnessFV_pfCombinedMVAV2BJetTags.clear();
  //Num_of_trks
  BJetnessFV_num_leps.clear();
  BJetnessFV_npvTrkOVcollTrk.clear();
  //ImpactParameter
  BJetnessFV_avip3d_val.clear();
  BJetnessFV_avip3d_sig.clear();
  BJetnessFV_avsip3d_sig.clear();
  BJetnessFV_avip1d_sig.clear();
}
/////
//   Methods to be aligned to the TTHbb selection
/////
//Ask for good vertices
bool BJetnessFVSelector::isGoodVertex(const reco::Vertex& vtx){
  if(vtx.isFake())                                   return false;
  if(vtx.ndof()<_vtx_ndof_min)                       return false;
  if(vtx.position().Rho()>_vtx_rho_max)              return false;
  if(fabs(vtx.position().Z()) > _vtx_position_z_max) return false;
  return true;
}
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
bool BJetnessFVSelector::is_tight_muon(const pat::Muon& mu, const reco::Vertex& vtx){
  bool isloosemu = false;
  if(mu.pt()>25 &&
    TMath::Abs(mu.eta()) < 2.1 &&
    mu.isTightMuon(vtx) &&
    rel_iso_dbc_mu(mu) < 0.15
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
bool BJetnessFVSelector::is_tight_electron(const pat::Electron& ele, double rhopog){
  bool isele = false;
  if(ele.pt()>30 && TMath::Abs(ele.eta())<2.1 &&
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
double BJetnessFVSelector::rel_iso_dbc_ele(const pat::Electron& el, double rhopog){
  double SumChHadPt       = el.pfIsolationVariables().sumChargedHadronPt;
  double SumNeuHadEt      = el.pfIsolationVariables().sumNeutralHadronEt;
  double SumPhotonEt      = el.pfIsolationVariables().sumPhotonEt;
  double EffArea          = get_effarea(el.superCluster()->position().eta());
  double SumNeutralCorrEt = std::max( 0.0, SumNeuHadEt+SumPhotonEt - rhopog*EffArea );
  double relIsoRhoEA = (SumChHadPt + SumNeutralCorrEt)/el.pt();
  return relIsoRhoEA;
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
//This function has to be updated in order to select good jet using the TTHbb definition
bool BJetnessFVSelector::is_good_jet(const pat::Jet &j,double rho, double rhoJER, int vtxsize){
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
  if(jetpt < 30)       isgoodjet = false; //Please note that this requirement is for the SL channel, while for DL channel we require pT > 20!
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
void BJetnessFVSelector::JECInitialization(){
  //AK4chs - MC: Get the factorized jet corrector parameters.
  std::vector<std::string> jecPayloadNamesAK4PFchsMC_;
  jecPayloadNamesAK4PFchsMC_.push_back(jecPayloadNamesAK4PFchsMC1_.fullPath());
  jecPayloadNamesAK4PFchsMC_.push_back(jecPayloadNamesAK4PFchsMC2_.fullPath());
  jecPayloadNamesAK4PFchsMC_.push_back(jecPayloadNamesAK4PFchsMC3_.fullPath());
  std::vector<JetCorrectorParameters> vParAK4PFchsMC;
  for ( std::vector<std::string>::const_iterator payloadBegin = jecPayloadNamesAK4PFchsMC_.begin(),
          payloadEnd = jecPayloadNamesAK4PFchsMC_.end(), ipayload = payloadBegin; ipayload != payloadEnd; ++ipayload ) {
    JetCorrectorParameters pars(*ipayload);
    vParAK4PFchsMC.push_back(pars);
  }
  jecAK4PFchsMC_    = boost::shared_ptr<FactorizedJetCorrector>  ( new FactorizedJetCorrector(vParAK4PFchsMC) );
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
void BJetnessFVSelector::GetJER(pat::Jet jet, float JesSF, float rhoJER, bool AK4PFchs, float &JERScaleFactor, float &JERScaleFactorUP, float &JERScaleFactorDOWN){
  if(!jet.genJet()) return;
  double jetEta=fabs(jet.eta());
  double cFactorJER = 1.0;
  double cFactorJERdown = 1.0;
  double cFactorJERup = 1.0;
  //https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution#JER_Scaling_factors_and_Uncertai
  //NEEDS UPDATING
  if( jetEta<0.5 ){
    cFactorJER = 1.122;
    cFactorJERdown = 1.122-0.026;
    cFactorJERup   = 1.122+0.026;
  } else if( jetEta<0.8 ){
    cFactorJER = 1.167;
    cFactorJERdown = 1.167-0.048;
    cFactorJERup   = 1.167+0.048;
  } else if( jetEta<1.1 ){
    cFactorJER = 1.168;
    cFactorJERdown = 1.168-0.046;
    cFactorJERup   = 1.168+0.046;
  } else if( jetEta<1.3 ){
    cFactorJER = 1.029;
    cFactorJERdown = 1.029-0.066;
    cFactorJERup   = 1.029+0.066;
  } else if( jetEta<1.7 ){
    cFactorJER = 1.115;
    cFactorJERdown = 1.115-0.030;
    cFactorJERup   = 1.115+0.030;
  } else if( jetEta<1.9 ){
    cFactorJER = 1.041;
    cFactorJERdown = 1.041-0.062;
    cFactorJERup   = 1.041+0.062;
  } else if( jetEta<2.1 ){
    cFactorJER = 1.167;
    cFactorJERdown = 1.167-0.086;
    cFactorJERup   = 1.167+0.086;
  } else if( jetEta<2.3 ){
    cFactorJER = 1.094;
    cFactorJERdown = 1.094-0.093;
    cFactorJERup   = 1.094+0.093;
  } else if( jetEta<2.5 ){
    cFactorJER = 1.168;
    cFactorJERdown = 1.168-0.120;
    cFactorJERup   = 1.168+0.120;
  } else if( jetEta<2.8 ){
    cFactorJER = 1.266;
    cFactorJERdown = 1.266-0.132;
    cFactorJERup   = 1.266+0.132;
  } else if( jetEta<3.0 ){
    cFactorJER = 1.595;
    cFactorJERdown = 1.595-0.175;
    cFactorJERup   = 1.595+0.175;
  } else if( jetEta<3.2 ){
    cFactorJER = 0.998;
    cFactorJERdown = 0.998-0.066;
    cFactorJERup   = 0.998+0.066;
  } else if( jetEta<5.0 ){
    cFactorJER = 1.226;
    cFactorJERdown = 1.226-0.145;
    cFactorJERup   = 1.226+0.145;
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
/////
//   Methods for the BJetness variables
/////
void BJetnessFVSelector::get_bjetness_vars(
                                           vector<pat::Jet> evtjets, const reco::Vertex& vtx, const TransientTrackBuilder& ttrkbuilder, edm::Handle<edm::View<pat::Electron> > electron_pat, edm::Handle<edm::View<pat::Muon> > muon_h,
                                           double& bjetnessFV_num_leps, double& bjetnessFV_npvTrkOVcollTrk, double& bjetnessFV_avip3d_val, double& bjetnessFV_avip3d_sig, double& bjetnessFV_avsip3d_sig, double& bjetnessFV_avip1d_sig
                                          ){
  //Get BJetness trk info
  vector<Track> jetschtrks; jetschtrks.clear();
  double num_pvtrks  = 0;
  double num_npvtrks = 0;
  double num_eles    = 0;
  double num_mus     = 0;
  vector<tuple<double, double, double> > jetsdir; jetsdir.clear();
  get_bjetness_trkinfos(evtjets, vtx, jetschtrks, num_pvtrks, num_npvtrks, electron_pat, muon_h, num_eles, num_mus, jetsdir);
  bjetnessFV_num_leps = num_eles+num_mus;
  if(jetschtrks.size()!=0){
    bjetnessFV_npvTrkOVcollTrk       = num_npvtrks/double(jetschtrks.size());
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
    bjetnessFV_npvTrkOVcollTrk       = -998;
    bjetnessFV_avip3d_val            = -998;
    bjetnessFV_avip3d_sig            = -998;
    bjetnessFV_avsip3d_sig           = -998;
    bjetnessFV_avip1d_sig            = -998;
  }
}
//Get the BJetness trk info
void BJetnessFVSelector::get_bjetness_trkinfos(vector<pat::Jet> evtjets, const reco::Vertex& vtx, vector<Track>& jetchtrks, double& bjetness_num_pvtrks, double& bjetness_num_npvtrks, edm::Handle<edm::View<pat::Electron> > electron_pat, edm::Handle<edm::View<pat::Muon> > muon_h, double& bjetness_num_eles, double& bjetness_num_mus, vector<tuple<double, double, double> >& jetsdir){
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
        if(fabs(jcand.pdgId())==13 && is_loosePOG_jetmuon(jcand,muon_h)) bjetness_num_mus++;
        if(fabs(jcand.pdgId())==11 && is_softLep_jetelectron(jcand,electron_pat,vtx)) bjetness_num_eles++;
        //if(fabs(jcand.pdgId())==11 && is_loosePOGNoIPNoIso_jetelectron(jcand,electron_pat,vtx)) bjetness_num_eles++;
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
bool BJetnessFVSelector::is_softLep_jetelectron(const pat::PackedCandidate &jcand, edm::Handle<edm::View<pat::Electron> > electron_pat, const reco::Vertex& vtx){
  bool isele = false;
  for(edm::View<pat::Electron>::const_iterator ele = electron_pat->begin(); ele != electron_pat->end(); ele++){
    const pat::Electron &lele = *ele;
    if(deltaR(jcand.p4(),lele.p4())<0.1 && fabs(jcand.pt()-lele.pt())/lele.pt()<0.05 ){
      const HitPattern &hitPattern = lele.gsfTrack().get()->hitPattern();
      uint32_t hit = hitPattern.getHitPattern(HitPattern::TRACK_HITS, 0);
      bool hitCondition = !(HitPattern::validHitFilter(hit) && ((HitPattern::pixelBarrelHitFilter(hit) && HitPattern::getLayer(hit) < 3) || HitPattern::pixelEndcapHitFilter(hit)));
      if(!hitCondition && lele.passConversionVeto()) isele = true;
      if(isele) break;
    }
  }
  return isele;
}
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
