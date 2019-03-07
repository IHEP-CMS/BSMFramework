#include "BSMFramework/BSM3G_TNT_Maker/interface/MuonSelector.h"
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
//#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "RecoBTag/BTagTools/interface/SignedTransverseImpactParameter.h"
#include "TMath.h"
KalmanVertexFitter vertexfittermu(true);
using namespace std;
using namespace pat;
using namespace edm;
using namespace reco;
namespace{
  struct ByEta{
    bool operator()(const pat::PackedCandidate *c1, const pat::PackedCandidate *c2) const{
      return c1->eta()<c2->eta();
    }
    bool operator()(double c1eta, const pat::PackedCandidate *c2) const{
      return c1eta<c2->eta();
    }
    bool operator()(const pat::PackedCandidate *c1, double c2eta) const{
      return c1->eta()<c2eta;
    }
  };
}

MuonSelector::MuonSelector(std::string name, TTree* tree, bool debug, const pset& iConfig, edm::ConsumesCollector && ic):
  baseTree(name,tree,debug),
  triggerBits_(ic.consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"))),
  triggerObjects_(ic.consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("objects")))
{
  vtx_h_              = ic.consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"));
  beamSpot_           = ic.consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpot"));
  muon_h_             = ic.consumes<edm::View<pat::Muon> >(iConfig.getParameter<edm::InputTag>("muons"));
  pfToken_            = ic.consumes<pat::PackedCandidateCollection>(edm::InputTag("packedPFCandidates"));
  _Muon_pt_min        = iConfig.getParameter<double>("Muon_pt_min");
  _Muon_eta_max       = iConfig.getParameter<double>("Muon_eta_max");
  _vtx_ndof_min       = iConfig.getParameter<int>("vtx_ndof_min");
  _vtx_rho_max        = iConfig.getParameter<int>("vtx_rho_max");
  _vtx_position_z_max = iConfig.getParameter<double>("vtx_position_z_max");
  _is_data            = iConfig.getParameter<bool>("is_data");
  _reduced            = iConfig.getParameter<bool>("reduced");
  SetBranches();
}
MuonSelector::~MuonSelector(){
  delete tree_;
}
void MuonSelector::Fill(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  Clear();
  /////
  //   Recall collections
  ///// 
  edm::Handle<reco::VertexCollection> vtx_h;
  iEvent.getByToken(vtx_h_, vtx_h);
  reco::BeamSpot beamSpot;
  edm::Handle<reco::BeamSpot> beamSpotHandle;
  iEvent.getByToken(beamSpot_, beamSpotHandle);
  edm::Handle<edm::View<pat::Muon> > muon_h;
  iEvent.getByToken(muon_h_, muon_h);
  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
  iEvent.getByToken(triggerObjects_, triggerObjects);
  edm::Handle<edm::TriggerResults> triggerBits;
  iEvent.getByToken(triggerBits_, triggerBits);
  edm::Handle<pat::PackedCandidateCollection> pcc;
  iEvent.getByToken(pfToken_, pcc);
  edm::ESHandle<TransientTrackBuilder> ttrkbuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",ttrkbuilder);
  /////
  //   Require a good vertex 
  ///// 
  //reco::VertexCollection::const_iterator firstgoodVertex = vtx_h->end();
  //for(reco::VertexCollection::const_iterator it = vtx_h->begin(); it != firstgoodVertex; it++){
  //  if(isGoodVertex(*it)){
  //    firstgoodVertex = it;
  //    break;
  //  }
  //}
  //if(firstgoodVertex == vtx_h->end()) return;
  //const reco::Vertex &firstGoodVertex = *firstgoodVertex;
  //if(vtx_h->empty()) return; // skip the event if no PV found
  const reco::Vertex &firstGoodVertex = vtx_h->front();  
  //bool isgoodvtx = isGoodVertex(firstGoodVertex);
  //if(!isgoodvtx) return;
  /////
  //   Get muon information
  /////
  for(edm::View<pat::Muon>::const_iterator mu = muon_h->begin(); mu != muon_h->end(); mu++){
    /////
    //   BSM variables
    /////
    //Acceptance 
    if(mu->pt() < _Muon_pt_min)         continue;
    if(fabs(mu->eta()) > _Muon_eta_max) continue;  
    //Kinematics
    Muon_pt.push_back(mu->pt());
    Muon_eta.push_back(mu->eta());
    Muon_phi.push_back(mu->phi());
    Muon_energy.push_back(mu->energy());
    Muon_px.push_back(mu->px());
    Muon_py.push_back(mu->py());
    Muon_pz.push_back(mu->pz());
    Muon_p.push_back(mu->p());
    Muon_dB.push_back(mu->dB());
    //Charge
    Muon_charge.push_back(mu->charge());
    //ID
    Muon_soft.push_back(mu->isSoftMuon(firstGoodVertex));
    Muon_loose.push_back(mu->isLooseMuon());
    Muon_medium.push_back(mu->isMediumMuon());
    Muon_tight.push_back(mu->isTightMuon(firstGoodVertex));
    Muon_isHighPt.push_back(mu->isHighPtMuon(firstGoodVertex));
    Muon_POGisGood.push_back(muon::isGoodMuon(*mu, muon::TMOneStationTight));
    Muon_pdgId.push_back(mu->pdgId());
    Muon_pf.push_back(mu->isPFMuon());   
    Muon_isGlobal.push_back(mu->isGlobalMuon());   
    Muon_isTrackerMuon.push_back(mu->isTrackerMuon());
    reco::TrackRef tunePBestTrack = mu->tunePMuonBestTrack();
    if(tunePBestTrack.isNonnull()){
      reco::Muon::MuonTrackType tunePBestTrackType = mu->tunePMuonBestTrackType();
      Muon_tunePBestTrackType.push_back(tunePBestTrackType);
    }else{
      Muon_tunePBestTrackType.push_back(-999);
    }
    //Isolation
    double SumChHadPt  = mu->pfIsolationR04().sumChargedHadronPt;
    double SumNeuHadEt = mu->pfIsolationR04().sumNeutralHadronEt;
    double SumPhotonEt = mu->pfIsolationR04().sumPhotonEt;
    double SumPU       = mu->pfIsolationR04().sumPUPt;
    Muon_isoR04Charged.push_back(SumChHadPt);
    Muon_isoR04NeutralHadron.push_back(SumNeuHadEt);
    Muon_isoR04Photon.push_back(SumPhotonEt);
    Muon_isoR04PU.push_back(SumPU);
    double SumNeutralCorrEt = std::max( 0.0, SumNeuHadEt+SumPhotonEt - 0.5*SumPU );
    double relIsoDeltaBeta = (SumChHadPt + SumNeutralCorrEt)/mu->pt();
    Muon_relIsoDeltaBetaR04.push_back(relIsoDeltaBeta);
    Muon_isoR04CharParPt.push_back((mu->pfIsolationR04().sumChargedParticlePt));
    SumChHadPt  = mu->pfIsolationR03().sumChargedHadronPt;
    SumNeuHadEt = mu->pfIsolationR03().sumNeutralHadronEt;
    SumPhotonEt = mu->pfIsolationR03().sumPhotonEt;
    SumPU       = mu->pfIsolationR03().sumPUPt;
    Muon_isoR03Charged.push_back(SumChHadPt);
    Muon_isoR03NeutralHadron.push_back(SumNeuHadEt);
    Muon_isoR03Photon.push_back(SumPhotonEt);
    Muon_isoR03PU.push_back(SumPU);
    SumNeutralCorrEt = std::max( 0.0, SumNeuHadEt+SumPhotonEt - 0.5*SumPU );
    relIsoDeltaBeta = (SumChHadPt + SumNeutralCorrEt)/mu->pt();
    Muon_relIsoDeltaBetaR03.push_back(relIsoDeltaBeta);
    Muon_isoR03CharParPt.push_back((mu->pfIsolationR03().sumChargedParticlePt));
    Muon_trackIso.push_back(mu->trackIso());
    Muon_TrackerIso.push_back(mu->isolationR03().sumPt/mu->pt());
    Muon_ecalIso.push_back(mu->ecalIso());
    Muon_hcalIso.push_back(mu->hcalIso());
    Muon_caloIso.push_back(mu->caloIso());
    Muon_isoSum.push_back((mu->trackIso() + mu->ecalIso() + mu->hcalIso()));
    Muon_pfEcalEnergy.push_back(mu->pfEcalEnergy());
    //Track related variables 
    reco::TrackRef gtk = mu->globalTrack();
    if(gtk.isNonnull()) Muon_chi2.push_back(gtk->normalizedChi2());
    else                Muon_chi2.push_back(-999);
    Muon_chi2LocalPosition.push_back(mu->combinedQuality().chi2LocalPosition); 
    Muon_matchedStat.push_back(mu->numberOfMatchedStations());
    if(gtk.isNonnull()) Muon_validHits.push_back(gtk->hitPattern().numberOfValidMuonHits()); 
    else                Muon_validHits.push_back(-999);  
    if(mu->innerTrack().isNonnull()){
      Muon_validHitsInner.push_back(mu->innerTrack()->hitPattern().numberOfValidPixelHits());
      Muon_TLayers.push_back(mu->innerTrack()->hitPattern().trackerLayersWithMeasurement());
      Muon_ndof.push_back(mu->innerTrack()->ndof());
      Muon_validFraction.push_back(mu->innerTrack()->validFraction());
      Muon_pixelLayersWithMeasurement.push_back(mu->innerTrack()->hitPattern().pixelLayersWithMeasurement());
      Muon_qualityhighPurity.push_back(mu->innerTrack()->quality(reco::TrackBase::highPurity));
    }else{
      Muon_validHitsInner.push_back(-999);
      Muon_TLayers.push_back(-999);
      Muon_ndof.push_back(-999);
      Muon_validFraction.push_back(-999);
      Muon_pixelLayersWithMeasurement.push_back(-999);
      Muon_qualityhighPurity.push_back(-999);
    }
    Muon_trkKink.push_back(mu->combinedQuality().trkKink);
    Muon_segmentCompatibility.push_back(mu->segmentCompatibility());
    //IP
    if(mu->innerTrack().isNonnull()){
      Muon_dz_pv.push_back(mu->innerTrack()->dz(firstGoodVertex.position()));
      Muon_dxy_pv.push_back(mu->innerTrack()->dxy(firstGoodVertex.position()));
      if(beamSpotHandle.isValid()){
        beamSpot = *beamSpotHandle;
        math::XYZPoint point(beamSpot.x0(),beamSpot.y0(), beamSpot.z0());
        Muon_dz_bs.push_back(mu->innerTrack()->dz(point));
        Muon_dxy_bs.push_back(-1.*(mu->innerTrack()->dxy(point)));
      }else{
        Muon_dz_bs.push_back(-998);
        Muon_dxy_bs.push_back(-998);
        edm::LogInfo("MyAnalyzer") << "No beam spot available from EventSetup \n";
      }
      Muon_dzError.push_back(mu->innerTrack()->dzError());
      Muon_dxyError.push_back(mu->innerTrack()->d0Error());
      Muon_vtx.push_back(mu->innerTrack()->vx());
      Muon_vty.push_back(mu->innerTrack()->vy());
      Muon_vtz.push_back(mu->innerTrack()->vz());
    }else{
      Muon_dz_pv.push_back(-999);
      Muon_dxy_pv.push_back(-999);
      Muon_dz_bs.push_back(-999);
      Muon_dxy_bs.push_back(-999);
      Muon_dzError.push_back(-999);
      Muon_dxyError.push_back(-999);
      Muon_vtx.push_back(-999);
      Muon_vty.push_back(-999);
      Muon_vtz.push_back(-999);
    }
    /////
    //   MC info
    /////
    if(!_is_data){
      const reco::GenParticle * genpart = mu->genParticle(); 
      if(genpart){
        Muon_gen_pt.push_back(genpart->pt());
        Muon_gen_eta.push_back(genpart->eta());
        Muon_gen_phi.push_back(genpart->phi());
        Muon_gen_en.push_back(genpart->energy());
        Muon_gen_pdgId.push_back(genpart->pdgId());
        Muon_gen_isPromptFinalState.push_back(genpart->isPromptFinalState());
        Muon_gen_isDirectPromptTauDecayProductFinalState.push_back(genpart->isDirectPromptTauDecayProductFinalState());
      }else{
        Muon_gen_pt.push_back(-999);
        Muon_gen_eta.push_back(-999);
        Muon_gen_phi.push_back(-999);
        Muon_gen_en.push_back(-999);
        Muon_gen_pdgId.push_back(-999);
        Muon_gen_isPromptFinalState.push_back(-999);
        Muon_gen_isDirectPromptTauDecayProductFinalState.push_back(-999);
      }
    }
    //TRIGGER MATCHING - for now only * is implemented
    int isMatchedToTrigger = 0;
    //if(!_is_data) isMatchedToTrigger = MatchingToTrigger(iEvent, triggerObjects, triggerBits, mu->eta(), mu->phi());
    Muon_isMatchedToTrigger.push_back(isMatchedToTrigger);
  }
}
void MuonSelector::SetBranches(){
  if(debug_)    std::cout<<"setting branches: calling AddBranch of baseTree"<<std::endl;
  //Kinematics
  AddBranch(&Muon_pt                ,"Muon_pt");
  AddBranch(&Muon_eta               ,"Muon_eta");
  AddBranch(&Muon_phi               ,"Muon_phi");
  AddBranch(&Muon_energy            ,"Muon_energy");
  if(!_reduced){
    AddBranch(&Muon_px                ,"Muon_px");
    AddBranch(&Muon_py                ,"Muon_py");
    AddBranch(&Muon_pz                ,"Muon_pz");
    AddBranch(&Muon_p                 ,"Muon_p");
    AddBranch(&Muon_dB                ,"Muon_dB");
  }
  //Charge
  AddBranch(&Muon_charge             ,"Muon_charge");
  //ID
  AddBranch(&Muon_soft               ,"Muon_soft");
  AddBranch(&Muon_loose              ,"Muon_loose");
  AddBranch(&Muon_medium             ,"Muon_medium");
  AddBranch(&Muon_tight              ,"Muon_tight");
  AddBranch(&Muon_isHighPt           ,"Muon_isHighPt");
  if(!_reduced){
    AddBranch(&Muon_POGisGood          ,"Muon_POGisGood");
    AddBranch(&Muon_pdgId              ,"Muon_pdgId");
    AddBranch(&Muon_pf                 ,"Muon_pf");
    AddBranch(&Muon_isGlobal           ,"Muon_isGlobal");
    AddBranch(&Muon_isTrackerMuon      ,"Muon_isTrackerMuon");
    AddBranch(&Muon_tunePBestTrackType ,"Muon_tunePBestTrackType");
    //Isolation
    AddBranch(&Muon_isoR04Charged       ,"Muon_isoR04Charged");
    AddBranch(&Muon_isoR04NeutralHadron ,"Muon_isoR04NeutralHadron");
    AddBranch(&Muon_isoR04Photon        ,"Muon_isoR04Photon");
    AddBranch(&Muon_isoR04PU            ,"Muon_isoR04PU");
  }
  AddBranch(&Muon_relIsoDeltaBetaR04  ,"Muon_relIsoDeltaBetaR04");
  AddBranch(&Muon_TrackerIso          ,"Muon_TrackerIso");
  if(!_reduced){
    AddBranch(&Muon_isoR04CharParPt     ,"Muon_isoR04CharParPt");
    AddBranch(&Muon_isoR03Charged       ,"Muon_isoR03Charged");
    AddBranch(&Muon_isoR03NeutralHadron ,"Muon_isoR03NeutralHadron");
    AddBranch(&Muon_isoR03Photon        ,"Muon_isoR03Photon");
    AddBranch(&Muon_isoR03PU            ,"Muon_isoR03PU");
    AddBranch(&Muon_relIsoDeltaBetaR03  ,"Muon_relIsoDeltaBetaR03");
    AddBranch(&Muon_isoR03CharParPt     ,"Muon_isoR03CharParPt");
    AddBranch(&Muon_trackIso            ,"Muon_trackIso");
    AddBranch(&Muon_ecalIso             ,"Muon_ecalIso");
    AddBranch(&Muon_hcalIso             ,"Muon_hcalIso");
    AddBranch(&Muon_caloIso             ,"Muon_caloIso");
    AddBranch(&Muon_isoSum              ,"Muon_isoSum");
    AddBranch(&Muon_pfEcalEnergy        ,"Muon_pfEcalEnergy");
    //Track related variables and neutral part isolation
    AddBranch(&Muon_chi2                       ,"Muon_chi2");
    AddBranch(&Muon_chi2LocalPosition          ,"Muon_chi2LocalPosition");
    AddBranch(&Muon_matchedStat                ,"Muon_matchedStat");
    AddBranch(&Muon_validHits                  ,"Muon_validHits");
    AddBranch(&Muon_validHitsInner             ,"Muon_validHitsInner");
    AddBranch(&Muon_TLayers                    ,"Muon_TLayers");
    AddBranch(&Muon_ndof                       ,"Muon_ndof");
    AddBranch(&Muon_validFraction              ,"Muon_validFraction");
    AddBranch(&Muon_pixelLayersWithMeasurement ,"Muon_pixelLayersWithMeasurement");
    AddBranch(&Muon_qualityhighPurity          ,"Muon_qualityhighPurity");
    AddBranch(&Muon_trkKink                    ,"Muon_trkKink");
    AddBranch(&Muon_segmentCompatibility       ,"Muon_segmentCompatibility");
    //IP
    AddBranch(&Muon_dxy_pv                 ,"Muon_dxy_pv");
    AddBranch(&Muon_dz_pv                  ,"Muon_dz_pv");
    AddBranch(&Muon_dz_bs                  ,"Muon_dz_bs");
    AddBranch(&Muon_dxy_bs                 ,"Muon_dxy_bs");
    AddBranch(&Muon_dzError                ,"Muon_dzError");
    AddBranch(&Muon_dxyError               ,"Muon_dxyError");
    AddBranch(&Muon_vtx                    ,"Muon_vtx");
    AddBranch(&Muon_vty                    ,"Muon_vty");
    AddBranch(&Muon_vtz                    ,"Muon_vtz");
    //MC info
    if(!_is_data){
      AddBranch(&Muon_gen_pt                                      ,"Muon_gen_pt");
      AddBranch(&Muon_gen_eta                                     ,"Muon_gen_eta");
      AddBranch(&Muon_gen_phi                                     ,"Muon_gen_phi");
      AddBranch(&Muon_gen_en                                      ,"Muon_gen_en");
      AddBranch(&Muon_gen_pdgId                                   ,"Muon_gen_pdgId");
      AddBranch(&Muon_gen_isPromptFinalState                      ,"Muon_gen_isPromptFinalState");
      AddBranch(&Muon_gen_isDirectPromptTauDecayProductFinalState ,"Muon_gen_isDirectPromptTauDecayProductFinalState");
    }
  }
  AddBranch(&Muon_isMatchedToTrigger,"Muon_isMatchedToTrigger");
  if(debug_) std::cout<<"set branches"<<std::endl;
}
void MuonSelector::Clear(){
  //Kinematics  
  Muon_pt.clear();
  Muon_eta.clear();
  Muon_phi.clear();
  Muon_energy.clear();
  Muon_px.clear();
  Muon_py.clear();
  Muon_pz.clear();
  Muon_p.clear();
  Muon_dB.clear();
  //Charge
  Muon_charge.clear(); 
  //ID
  Muon_soft.clear();
  Muon_loose.clear();
  Muon_medium.clear();
  Muon_tight.clear();
  Muon_isHighPt.clear();
  Muon_POGisGood.clear();
  Muon_pdgId.clear();
  Muon_pf.clear();   
  Muon_isGlobal.clear();   
  Muon_isTrackerMuon.clear();
  Muon_tunePBestTrackType.clear();
  //Isolation
  Muon_isoR04Charged.clear();
  Muon_isoR04NeutralHadron.clear();
  Muon_isoR04Photon.clear();
  Muon_isoR04PU.clear();
  Muon_relIsoDeltaBetaR04.clear();
  Muon_isoR04CharParPt.clear();
  Muon_isoR03Charged.clear();
  Muon_isoR03NeutralHadron.clear();
  Muon_isoR03Photon.clear();
  Muon_isoR03PU.clear();
  Muon_relIsoDeltaBetaR03.clear();
  Muon_isoR03CharParPt.clear();
  Muon_trackIso.clear();
  Muon_TrackerIso.clear();
  Muon_ecalIso.clear();
  Muon_hcalIso.clear(); 
  Muon_caloIso.clear();
  Muon_isoSum.clear();
  Muon_pfEcalEnergy.clear();
  //Track related variables
  Muon_chi2.clear(); 
  Muon_chi2LocalPosition.clear();
  Muon_matchedStat.clear(); 
  Muon_validHits.clear();
  Muon_validHitsInner.clear(); 
  Muon_TLayers.clear(); 
  Muon_ndof.clear();
  Muon_validFraction.clear();
  Muon_pixelLayersWithMeasurement.clear();
  Muon_qualityhighPurity.clear();
  Muon_trkKink.clear();
  Muon_segmentCompatibility.clear();
  //IP
  Muon_dz_pv.clear();
  Muon_dxy_pv.clear(); 
  Muon_dz_bs.clear();
  Muon_dxy_bs.clear();
  Muon_dzError.clear();
  Muon_dxyError.clear();
  Muon_vtx.clear();
  Muon_vty.clear();
  Muon_vtz.clear();
  //MC info
  if(!_is_data){
    Muon_gen_pt.clear();
    Muon_gen_eta.clear();
    Muon_gen_phi.clear();
    Muon_gen_en.clear();
    Muon_gen_pdgId.clear();
    Muon_gen_isPromptFinalState.clear();
    Muon_gen_isDirectPromptTauDecayProductFinalState.clear();
  }
  Muon_isMatchedToTrigger.clear();
}
bool MuonSelector::isGoodVertex(const reco::Vertex& vtx){
  if(vtx.isFake())                                   return false;
  if(vtx.ndof()<_vtx_ndof_min)                       return false;
  if(vtx.position().Rho()>_vtx_rho_max)              return false;
  if(fabs(vtx.position().Z()) > _vtx_position_z_max) return false;
  return true;
}
double MuonSelector::get_effarea(double eta){
  double effarea = -1;
  if(abs(eta) < 0.8)      effarea = 0.0566;
  else if(abs(eta) < 1.3) effarea = 0.0562;
  else if(abs(eta) < 2.0) effarea = 0.0363;
  else if(abs(eta) < 2.2) effarea = 0.0119;
  else                    effarea = 0.0064;
  return effarea;
}
int MuonSelector::MatchingToTrigger(const edm::Event& iEvent, edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects, edm::Handle<edm::TriggerResults> triggerBits, float eta, float phi){
  const edm::TriggerNames & triggerNames = iEvent.triggerNames(*triggerBits);
  float deltaRMin = 99.;
  for (pat::TriggerObjectStandAlone obj : *triggerObjects) { 
    obj.unpackPathNames(triggerNames);
    if (obj.hasFilterLabel("hltL3fL1sMu22L1f0Tkf24QL3trkIsoFiltered0p09") || obj.hasFilterLabel("hltL3crIsoL1sMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p09")) {
      //for (unsigned h = 0; h < obj.filterLabels().size(); ++h) std::cout << " " << obj.filterLabels()[h];  
      //cout<<obj.pt()<<endl;
      float deltaPhi = TMath::Abs(phi-obj.phi());
      float deltaEta = eta-obj.eta();
      if(deltaPhi > TMath::Pi()) deltaPhi = TMath::TwoPi() - deltaPhi;
      float deltaR = TMath::Sqrt(deltaEta*deltaEta + deltaPhi*deltaPhi);
      if(deltaR<deltaRMin) deltaRMin=deltaR;
    }
  }
  int result = 0;
  if(deltaRMin<0.3) result = 1;
  return result;
}
