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
#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"
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
  baseTree(name,tree,debug)
{
  vtx_h_              = ic.consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"));
  beamSpot_           = ic.consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpot"));
  muon_h_             = ic.consumes<edm::View<pat::Muon> >(iConfig.getParameter<edm::InputTag>("muons"));
  jets_               = ic.consumes<pat::JetCollection >(iConfig.getParameter<edm::InputTag>("jets"));
  jetsToken           = ic.consumes<edm::View<pat::Jet>>(iConfig.getParameter<edm::InputTag>("jets"));
  pfToken_            = ic.consumes<pat::PackedCandidateCollection>(edm::InputTag("packedPFCandidates"));
  rhoHandle_          = ic.consumes<double>(edm::InputTag("fixedGridRhoFastjetCentralNeutral"));
  qgToken_            = ic.consumes<edm::ValueMap<float>>(edm::InputTag("QGTagger", "qgLikelihood"));
  _Muon_pt_min        = iConfig.getParameter<double>("Muon_pt_min");
  _Muon_eta_max       = iConfig.getParameter<double>("Muon_eta_max");
  _vtx_ndof_min       = iConfig.getParameter<int>("vtx_ndof_min");
  _vtx_rho_max        = iConfig.getParameter<int>("vtx_rho_max");
  _vtx_position_z_max = iConfig.getParameter<double>("vtx_position_z_max");
  _super_TNT          = iConfig.getParameter<bool>("super_TNT");
  _is_data            = iConfig.getParameter<bool>("is_data");
  _AJVar              = iConfig.getParameter<bool>("AJVar");
  _tthlepVar          = iConfig.getParameter<bool>("tthlepVar");
  _qglVar             = iConfig.getParameter<bool>("qglVar");
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
  edm::Handle<pat::JetCollection> jets;
  iEvent.getByToken(jets_, jets);
  edm::Handle<double> rhoHandle;
  iEvent.getByToken(rhoHandle_,rhoHandle);
  double rho = *rhoHandle;
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
  if(vtx_h->empty()) return; // skip the event if no PV found
  const reco::Vertex &firstGoodVertex = vtx_h->front();  
  bool isgoodvtx = isGoodVertex(firstGoodVertex);
  if(!isgoodvtx) return;
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
    if(mu->innerTrack().isNonnull()){
      Muon_pt_it.push_back(mu->innerTrack()->pt());
      Muon_ptErr_it.push_back(mu->innerTrack()->ptError());
      if(mu->innerTrack()->pt()!=0) Muon_pTErrOVpT_it.push_back(mu->innerTrack()->ptError()/mu->innerTrack()->pt());
      else                          Muon_pTErrOVpT_it.push_back(-998);   
    }else{
      Muon_pt_it.push_back(-999);
      Muon_ptErr_it.push_back(-999);
      Muon_pTErrOVpT_it.push_back(-999);
    }
    if(mu->muonBestTrack().isNonnull()){
      Muon_pt_bt.push_back(mu->muonBestTrack()->pt());
      Muon_ptErr_bt.push_back(mu->muonBestTrack()->ptError());
      if(mu->muonBestTrack()->pt()!=0) Muon_pTErrOVpT_bt.push_back(mu->muonBestTrack()->ptError()/mu->muonBestTrack()->pt());
      else                             Muon_pTErrOVpT_bt.push_back(-998); 
    }else{
      Muon_pt_bt.push_back(-999);
      Muon_ptErr_bt.push_back(-999);
      Muon_pTErrOVpT_bt.push_back(-999);
    }
    reco::TrackRef tunePBestTrack = mu->tunePMuonBestTrack();
    if(tunePBestTrack.isNonnull()) Muon_pt_tunePbt.push_back(tunePBestTrack->pt());
    else                           Muon_pt_tunePbt.push_back(-999);
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
    if(_AJVar){
      if(beamSpotHandle.isValid() && mu->innerTrack().isNonnull()){//AJ vars (both pv and bs are in this if condition, tought for pv is not mandatory)
	beamSpot = *beamSpotHandle;
	GlobalPoint thebs(beamSpot.x0(),beamSpot.y0(),beamSpot.z0());
	GlobalPoint thepv(firstGoodVertex.position().x(),firstGoodVertex.position().y(),firstGoodVertex.position().z()); 
	TrackRef muit = mu->innerTrack();
	TransientTrack muonTransTkPtr = ttrkbuilder->build(muit);
	GlobalPoint mu_pca_bs = muonTransTkPtr.trajectoryStateClosestToPoint(thebs).position();
	GlobalPoint mu_pca_pv = muonTransTkPtr.trajectoryStateClosestToPoint(thepv).position();
	Muon_track_PCAx_pv.push_back(mu_pca_pv.x());
	Muon_track_PCAy_pv.push_back(mu_pca_pv.y());
	Muon_track_PCAz_pv.push_back(mu_pca_pv.z());
	Muon_track_PCAx_bs.push_back(mu_pca_bs.x());
	Muon_track_PCAy_bs.push_back(mu_pca_bs.y());
	Muon_track_PCAz_bs.push_back(mu_pca_bs.z());
	const float muonMass = 0.1056583715;
	float muonSigma      = muonMass*1e-6;
	float chi2 = 0.0;
	float ndf  = 0.0;
	KinematicParticleFactoryFromTransientTrack pFactory;
	RefCountedKinematicParticle muonParticle = pFactory.particle(muonTransTkPtr, muonMass, chi2, ndf, muonSigma);
	Muon_trackFitErrorMatrix_00.push_back(muonParticle->stateAtPoint(mu_pca_bs).kinematicParametersError().matrix()(0,0));
	Muon_trackFitErrorMatrix_01.push_back(muonParticle->stateAtPoint(mu_pca_bs).kinematicParametersError().matrix()(0,1));
	Muon_trackFitErrorMatrix_02.push_back(muonParticle->stateAtPoint(mu_pca_bs).kinematicParametersError().matrix()(0,2));
	Muon_trackFitErrorMatrix_11.push_back(muonParticle->stateAtPoint(mu_pca_bs).kinematicParametersError().matrix()(1,1));
	Muon_trackFitErrorMatrix_12.push_back(muonParticle->stateAtPoint(mu_pca_bs).kinematicParametersError().matrix()(1,2));
	Muon_trackFitErrorMatrix_22.push_back(muonParticle->stateAtPoint(mu_pca_bs).kinematicParametersError().matrix()(2,2));
      }else{
	Muon_track_PCAx_bs.push_back(-999);
	Muon_track_PCAy_bs.push_back(-999);
	Muon_track_PCAz_bs.push_back(-999);
	Muon_track_PCAx_pv.push_back(-999);
	Muon_track_PCAy_pv.push_back(-999);
	Muon_track_PCAz_pv.push_back(-999);
	Muon_trackFitErrorMatrix_00.push_back(-999);
	Muon_trackFitErrorMatrix_01.push_back(-999);
	Muon_trackFitErrorMatrix_02.push_back(-999);
	Muon_trackFitErrorMatrix_11.push_back(-999);
	Muon_trackFitErrorMatrix_12.push_back(-999);
	Muon_trackFitErrorMatrix_22.push_back(-999);
      }
      if(mu->muonBestTrack().isNonnull()){
	Muon_dz_bt.push_back(fabs(mu->muonBestTrack()->dz(firstGoodVertex.position())));
	Muon_dxy_bt.push_back(fabs(mu->muonBestTrack()->dxy(firstGoodVertex.position())));
      }else{
	Muon_dz_bt.push_back(-999);
	Muon_dxy_bt.push_back(-999);
      } 
    }
    //////
    ///   TTH variables
    ////// 
    if(_tthlepVar){
      //Isolation
      double miniIso      = 999;
      double miniIsoCh    = 999;
      double miniIsoNeu   = 999;
      double miniIsoPUsub = 999;
      get_muminiIso_info(*pcc,rho,*mu,miniIso,miniIsoCh,miniIsoNeu,miniIsoPUsub);
      Muon_miniIsoRel.push_back(miniIso/mu->pt());
      Muon_miniIsoCh.push_back(miniIsoCh);
      Muon_miniIsoNeu.push_back(miniIsoNeu);
      Muon_miniIsoPUsub.push_back(miniIsoPUsub);
      //Other
      double mujet_mindr    = 999;
      double mujet_pt       = -1;
      double muptOVmujetpt  = -1;
      double mujet_pfCombinedInclusiveSecondaryVertexV2BJetTags = -1;
      double mujet_pfJetProbabilityBJetTag = -1;
      double mujet_pfCombinedMVABJetTags = -1;
      double mujet_qgl = -3;
      double mujetx  = -999;
      double mujety  = -999;
      double mujetz  = -999;
      double muptrel = -999;
      double muparel = -999;
      double muparelsub = -999;
      int lepjetidx = -1;
      get_mujet_info(*mu,iEvent,iSetup,
                     mujet_mindr,mujet_pt,muptOVmujetpt,
                     mujet_pfCombinedInclusiveSecondaryVertexV2BJetTags,mujet_pfJetProbabilityBJetTag,mujet_pfCombinedMVABJetTags,mujet_qgl,
                     mujetx,mujety,mujetz,muptrel,muparel,muparelsub,lepjetidx);
      Muon_jetdr.push_back(mujet_mindr);
      Muon_jetpt.push_back(mujet_pt);
      Muon_jetptratio.push_back(muptOVmujetpt);
      Muon_jetcsv.push_back(mujet_pfCombinedInclusiveSecondaryVertexV2BJetTags);
      Muon_ptrel.push_back(muptrel);
      Muon_parel.push_back(muparel);
      Muon_parelsub.push_back(muparelsub);
      Muon_mujet_pfJetProbabilityBJetTag.push_back(mujet_pfJetProbabilityBJetTag);
      Muon_mujet_pfCombinedMVABJetTags.push_back(mujet_pfCombinedMVABJetTags);
      Muon_mujet_qgl.push_back(mujet_qgl);      
      Muon_IP3Dsig_it.push_back(fabs(mu->dB(pat::Muon::PV3D))/mu->edB(pat::Muon::PV3D));
      int pvass = pvassociation(*mu,*pcc);
      Muon_pvass.push_back(pvass);
      const math::XYZVector& lepton_momentum = mu->momentum(); 
      const math::XYZVector axis(mujetx,mujety,mujetz);
      double etarel = relEta(axis,lepton_momentum);
      Muon_etarel.push_back(etarel);
      double etarelsub = relEta(axis-lepton_momentum,lepton_momentum);
      Muon_etarelsub.push_back(etarelsub);
      double etapar = parEta(axis,lepton_momentum);
      Muon_etapar.push_back(etapar);
      double etaparsub = parEta(axis-lepton_momentum,lepton_momentum);
      Muon_etaparsub.push_back(etaparsub);
      Muon_ptOVen.push_back(mu->pt()/mu->energy());
      //Mass
      Muon_mumass.push_back(mu->p4().M());
      double muwmass    = 0;//get_lepWmass(*mu,iEvent,lepjetidx);
      double mutopmass  = 0;//get_lepTopmass(*mu,iEvent,lepjetidx);
      double muwtopmass = 0;//get_lepWTopmass(*mu,iEvent,lepjetidx);
      if(lepjetidx!=-1){
        const pat::Jet & lepjet = (*jets)[lepjetidx];
        Muon_mujet_mass.push_back(lepjet.p4().M());
      }else{
        Muon_mujet_mass.push_back(0);
      }
      Muon_mujet_Wmass.push_back(muwmass);
      Muon_mujet_Topmass.push_back(mutopmass);
      Muon_mujet_WTopmass.push_back(muwtopmass);
      //Mu IP 
      GlobalVector mujetgv(mujetx,mujety,mujetz);
      double IP3D_val  = -9999;
      double IP3D_err  = -9999;
      double IP3D_sig  = -9999;
      double sIP3D_val = -9999;
      double sIP3D_err = -9999;
      double sIP3D_sig = -9999;
      double IP2D_val  = -9999;
      double IP2D_err  = -9999;
      double IP2D_sig  = -9999;
      double sIP2D_val = -9999;
      double sIP2D_err = -9999;
      double sIP2D_sig = -9999;
      double IP1D_val  = -9999;
      double IP1D_err  = -9999;
      double IP1D_sig  = -9999;
      double sIP1D_val = -9999;
      double sIP1D_err = -9999;
      double sIP1D_sig = -9999;
      if(mu->innerTrack().isNonnull()){
        TrackRef muit = mu->innerTrack();
        TransientTrack muttrk = ttrkbuilder->build(muit);
        IP3D2D(muttrk,firstGoodVertex,mujetgv,IP3D_val,IP3D_err,IP3D_sig,sIP3D_val,sIP3D_err,sIP3D_sig,IP2D_val,IP2D_err,IP2D_sig,sIP2D_val,sIP2D_err,sIP2D_sig);
        zIP1D(muttrk,firstGoodVertex,mujetgv,IP1D_val,IP1D_err,IP1D_sig,sIP1D_val,sIP1D_err,sIP1D_sig);
      }
      //Max Lep jet IP (the maximum IP for a tracks of the lepton jet)
      double lepjetMaxIP3D_val  = IP3D_val; 
      double lepjetMaxIP3D_sig  = IP3D_sig; 
      double lepjetMaxsIP3D_val = sIP3D_val; 
      double lepjetMaxsIP3D_sig = sIP3D_sig; 
      double lepjetMaxIP2D_val  = IP2D_val; 
      double lepjetMaxIP2D_sig  = IP2D_sig; 
      double lepjetMaxsIP2D_val = sIP2D_val; 
      double lepjetMaxsIP2D_sig = sIP2D_sig; 
      double lepjetMaxIP1D_val  = IP1D_val; 
      double lepjetMaxIP1D_sig  = IP1D_sig; 
      double lepjetMaxsIP1D_val = sIP1D_val; 
      double lepjetMaxsIP1D_sig = sIP1D_sig; 
      //Av Lep jet IP (the average IP of the lepton jet tracks)   
      double lepjetAvIP3D_val  = 0; double denlepjetAvIP3D_val  = 0;
      if(IP3D_val!=-9999){lepjetAvIP3D_val   = IP3D_val; denlepjetAvIP3D_val   = 1;} 
      double lepjetAvIP3D_sig  = 0; double denlepjetAvIP3D_sig  = 0;
      if(IP3D_sig!=-9999){lepjetAvIP3D_sig   = IP3D_sig; denlepjetAvIP3D_sig   = 1;}
      double lepjetAvsIP3D_val = 0; double denlepjetAvsIP3D_val = 0;
      if(sIP3D_val!=-9999){lepjetAvsIP3D_val = sIP3D_val; denlepjetAvsIP3D_val = 1;} 
      double lepjetAvsIP3D_sig = 0; double denlepjetAvsIP3D_sig = 0;
      if(sIP3D_sig!=-9999){lepjetAvsIP3D_sig = sIP3D_sig; denlepjetAvsIP3D_sig = 1;} 
      double lepjetAvIP2D_val  = 0; double denlepjetAvIP2D_val  = 0;
      if(IP2D_val!=-9999){lepjetAvIP2D_val   = IP2D_val; denlepjetAvIP2D_val   = 1;} 
      double lepjetAvIP2D_sig  = 0; double denlepjetAvIP2D_sig  = 0;
      if(IP2D_sig!=-9999){lepjetAvIP2D_sig   = IP2D_sig; denlepjetAvIP2D_sig   = 1;} 
      double lepjetAvsIP2D_val = 0; double denlepjetAvsIP2D_val = 0;
      if(sIP2D_val!=-9999){lepjetAvsIP2D_val = sIP2D_val; denlepjetAvsIP2D_val = 1;} 
      double lepjetAvsIP2D_sig = 0; double denlepjetAvsIP2D_sig = 0;
      if(sIP2D_sig!=-9999){lepjetAvsIP2D_sig = sIP2D_sig; denlepjetAvsIP2D_sig = 1;} 
      double lepjetAvIP1D_val  = 0; double denlepjetAvIP1D_val  = 0;
      if(IP1D_val!=-9999){lepjetAvIP1D_val   = IP1D_val; denlepjetAvIP1D_val   = 1;} 
      double lepjetAvIP1D_sig  = 0; double denlepjetAvIP1D_sig  = 0;
      if(IP1D_sig!=-9999){lepjetAvIP1D_sig   = IP1D_sig; denlepjetAvIP1D_sig   = 1;} 
      double lepjetAvsIP1D_val = 0; double denlepjetAvsIP1D_val = 0;
      if(sIP1D_val!=-9999){lepjetAvsIP1D_val = sIP1D_val; denlepjetAvsIP1D_val = 1;} 
      double lepjetAvsIP1D_sig = 0; double denlepjetAvsIP1D_sig = 0;
      if(sIP1D_sig!=-9999){lepjetAvsIP1D_sig = sIP1D_sig; denlepjetAvsIP1D_sig = 1;} 
      //Get values of Max and Av IP
      if(lepjetidx!=-1){
        const pat::Jet & lepjet = (*jets)[lepjetidx]; 
        lepjetIP(lepjet,firstGoodVertex,mujetgv,*ttrkbuilder,
                 lepjetMaxIP3D_val, lepjetMaxIP3D_sig, lepjetMaxsIP3D_val, lepjetMaxsIP3D_sig, lepjetMaxIP2D_val, lepjetMaxIP2D_sig, lepjetMaxsIP2D_val, lepjetMaxsIP2D_sig, lepjetMaxIP1D_val, lepjetMaxIP1D_sig, lepjetMaxsIP1D_val, lepjetMaxsIP1D_sig,
                 lepjetAvIP3D_val, lepjetAvIP3D_sig, lepjetAvsIP3D_val, lepjetAvsIP3D_sig, lepjetAvIP2D_val, lepjetAvIP2D_sig, lepjetAvsIP2D_val, lepjetAvsIP2D_sig, lepjetAvIP1D_val, lepjetAvIP1D_sig, lepjetAvsIP1D_val, lepjetAvsIP1D_sig,
                 denlepjetAvIP3D_val, denlepjetAvIP3D_sig, denlepjetAvsIP3D_val, denlepjetAvsIP3D_sig, denlepjetAvIP2D_val, denlepjetAvIP2D_sig, denlepjetAvsIP2D_val, denlepjetAvsIP2D_sig, denlepjetAvIP1D_val, denlepjetAvIP1D_sig, denlepjetAvsIP1D_val, denlepjetAvsIP1D_sig,
                 IP3D_val   
                );                                                                    
      }
      Muon_IP3D_val.push_back(IP3D_val);
      Muon_IP3D_err.push_back(IP3D_err);
      Muon_IP3D_sig.push_back(IP3D_sig);
      Muon_sIP3D_val.push_back(sIP3D_val);
      Muon_sIP3D_err.push_back(sIP3D_err);
      Muon_sIP3D_sig.push_back(sIP3D_sig);
      Muon_IP2D_val.push_back(IP2D_val);
      Muon_IP2D_err.push_back(IP2D_err);
      Muon_IP2D_sig.push_back(IP2D_sig);
      Muon_sIP2D_val.push_back(sIP2D_val);
      Muon_sIP2D_err.push_back(sIP2D_err);
      Muon_sIP2D_sig.push_back(sIP2D_sig);
      Muon_IP1D_val.push_back(IP1D_val);
      Muon_IP1D_err.push_back(IP1D_err);
      Muon_IP1D_sig.push_back(IP1D_sig);
      Muon_sIP1D_val.push_back(sIP1D_val);
      Muon_sIP1D_err.push_back(sIP1D_err);
      Muon_sIP1D_sig.push_back(sIP1D_sig);
      Muon_lepjetMaxIP3D_val.push_back(lepjetMaxIP3D_val);
      Muon_lepjetMaxIP3D_sig.push_back(lepjetMaxIP3D_sig);
      Muon_lepjetMaxsIP3D_val.push_back(lepjetMaxsIP3D_val);
      Muon_lepjetMaxsIP3D_sig.push_back(lepjetMaxsIP3D_sig);
      Muon_lepjetMaxIP2D_val.push_back(lepjetMaxIP2D_val);
      Muon_lepjetMaxIP2D_sig.push_back(lepjetMaxIP2D_sig);
      Muon_lepjetMaxsIP2D_val.push_back(lepjetMaxsIP2D_val);
      Muon_lepjetMaxsIP2D_sig.push_back(lepjetMaxsIP2D_sig);
      Muon_lepjetMaxIP1D_val.push_back(lepjetMaxIP1D_val);
      Muon_lepjetMaxIP1D_sig.push_back(lepjetMaxIP1D_sig);
      Muon_lepjetMaxsIP1D_val.push_back(lepjetMaxsIP1D_val);
      Muon_lepjetMaxsIP1D_sig.push_back(lepjetMaxsIP1D_sig);
      Muon_lepjetAvIP3D_val.push_back(denlepjetAvIP3D_val!=0   ? lepjetAvIP3D_val/denlepjetAvIP3D_val   : IP3D_val);
      Muon_lepjetAvIP3D_sig.push_back(denlepjetAvIP3D_sig!=0   ? lepjetAvIP3D_sig/denlepjetAvIP3D_sig   : IP3D_sig);
      Muon_lepjetAvsIP3D_val.push_back(denlepjetAvsIP3D_val!=0 ? lepjetAvsIP3D_val/denlepjetAvsIP3D_val : sIP3D_val);
      Muon_lepjetAvsIP3D_sig.push_back(denlepjetAvsIP3D_sig!=0 ? lepjetAvsIP3D_sig/denlepjetAvsIP3D_sig : sIP3D_sig);
      Muon_lepjetAvIP2D_val.push_back(denlepjetAvIP2D_val!=0   ? lepjetAvIP2D_val/denlepjetAvIP2D_val   : IP2D_val);
      Muon_lepjetAvIP2D_sig.push_back(denlepjetAvIP2D_sig!=0   ? lepjetAvIP2D_sig/denlepjetAvIP2D_sig   : IP2D_sig);
      Muon_lepjetAvsIP2D_val.push_back(denlepjetAvsIP2D_val!=0 ? lepjetAvsIP2D_val/denlepjetAvsIP2D_val : sIP2D_val);
      Muon_lepjetAvsIP2D_sig.push_back(denlepjetAvsIP2D_sig!=0 ? lepjetAvsIP2D_sig/denlepjetAvsIP2D_sig : sIP2D_sig);
      Muon_lepjetAvIP1D_val.push_back(denlepjetAvIP1D_val!=0   ? lepjetAvIP1D_val/denlepjetAvIP1D_val   : IP1D_val);
      Muon_lepjetAvIP1D_sig.push_back(denlepjetAvIP1D_sig!=0   ? lepjetAvIP1D_sig/denlepjetAvIP1D_sig   : IP1D_sig);
      Muon_lepjetAvsIP1D_val.push_back(denlepjetAvsIP1D_val!=0 ? lepjetAvsIP1D_val/denlepjetAvsIP1D_val : sIP1D_val);
      Muon_lepjetAvsIP1D_sig.push_back(denlepjetAvsIP1D_sig!=0 ? lepjetAvsIP1D_sig/denlepjetAvsIP1D_sig : sIP1D_sig);
      //Lep jet trks
      double lepjetchtrks      = 0;
      double lepjetpvchtrks    = 0;
      double lepjetnonpvchtrks = 0;
      double lepjetndaus       = 0;
      if(lepjetidx!=-1){
        const pat::Jet & lepjet = (*jets)[lepjetidx];
        lepjetTrks(lepjet, firstGoodVertex, lepjetchtrks, lepjetpvchtrks, lepjetnonpvchtrks, lepjetndaus);
      }
      Muon_lepjetchtrks.push_back(lepjetchtrks);
      Muon_lepjetpvchtrks.push_back(lepjetpvchtrks);
      Muon_lepjetnonpvchtrks.push_back(lepjetnonpvchtrks);
      Muon_lepjetndaus.push_back(lepjetndaus);
      //Lep jet vtx compatibility
      double lepjetpvchi2    = 0;
      double lepjetnumno2trk = 0;
      if(lepjetidx!=-1){
        const pat::Jet & lepjet = (*jets)[lepjetidx];
        lepjetVtxCompatibility(lepjet, firstGoodVertex, *ttrkbuilder, lepjetpvchi2, lepjetnumno2trk);
      }
      Muon_lepjetpvchi2.push_back(lepjetpvchi2);
      Muon_lepjetnumno2trk.push_back(lepjetnumno2trk);
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
  }
}
void MuonSelector::SetBranches(){
  if(debug_)    std::cout<<"setting branches: calling AddBranch of baseTree"<<std::endl;
  //Kinematics
  AddBranch(&Muon_pt                ,"Muon_pt");
  AddBranch(&Muon_eta               ,"Muon_eta");
  AddBranch(&Muon_phi               ,"Muon_phi");
  AddBranch(&Muon_energy            ,"Muon_energy");
  AddBranch(&Muon_px                ,"Muon_px");
  AddBranch(&Muon_py                ,"Muon_py");
  AddBranch(&Muon_pz                ,"Muon_pz");
  AddBranch(&Muon_p                 ,"Muon_p");
  AddBranch(&Muon_dB                ,"Muon_dB");
  AddBranch(&Muon_pt_it             ,"Muon_pt_it");
  AddBranch(&Muon_ptErr_it          ,"Muon_ptErr_it");
  AddBranch(&Muon_pTErrOVpT_it      ,"Muon_pTErrOVpT_it");
  AddBranch(&Muon_pt_bt             ,"Muon_pt_bt");
  AddBranch(&Muon_ptErr_bt          ,"Muon_ptErr_bt");
  AddBranch(&Muon_pTErrOVpT_bt      ,"Muon_pTErrOVpT_bt");
  AddBranch(&Muon_pt_tunePbt        ,"Muon_pt_tunePbt");
  //Charge
  AddBranch(&Muon_charge             ,"Muon_charge");
  //ID
  AddBranch(&Muon_soft               ,"Muon_soft");
  AddBranch(&Muon_loose              ,"Muon_loose");
  AddBranch(&Muon_medium             ,"Muon_medium");
  AddBranch(&Muon_tight              ,"Muon_tight");
  AddBranch(&Muon_isHighPt           ,"Muon_isHighPt");
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
  AddBranch(&Muon_relIsoDeltaBetaR04  ,"Muon_relIsoDeltaBetaR04");
  AddBranch(&Muon_isoR04CharParPt     ,"Muon_isoR04CharParPt");
  AddBranch(&Muon_isoR03Charged       ,"Muon_isoR03Charged");
  AddBranch(&Muon_isoR03NeutralHadron ,"Muon_isoR03NeutralHadron");
  AddBranch(&Muon_isoR03Photon        ,"Muon_isoR03Photon");
  AddBranch(&Muon_isoR03PU            ,"Muon_isoR03PU");
  AddBranch(&Muon_relIsoDeltaBetaR03  ,"Muon_relIsoDeltaBetaR03");
  AddBranch(&Muon_isoR03CharParPt     ,"Muon_isoR03CharParPt");
  AddBranch(&Muon_trackIso            ,"Muon_trackIso");
  AddBranch(&Muon_TrackerIso          ,"Muon_TrackerIso");
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
  if(_AJVar){
    AddBranch(&Muon_track_PCAx_bs          ,"Muon_track_PCAx_bs");
    AddBranch(&Muon_track_PCAy_bs          ,"Muon_track_PCAy_bs");
    AddBranch(&Muon_track_PCAz_bs          ,"Muon_track_PCAz_bs");
    AddBranch(&Muon_track_PCAx_pv          ,"Muon_track_PCAx_pv");
    AddBranch(&Muon_track_PCAy_pv          ,"Muon_track_PCAy_pv");
    AddBranch(&Muon_track_PCAz_pv          ,"Muon_track_PCAz_pv");
    AddBranch(&Muon_trackFitErrorMatrix_00 ,"Muon_trackFitErrorMatrix_00");
    AddBranch(&Muon_trackFitErrorMatrix_01 ,"Muon_trackFitErrorMatrix_01");
    AddBranch(&Muon_trackFitErrorMatrix_02 ,"Muon_trackFitErrorMatrix_02");
    AddBranch(&Muon_trackFitErrorMatrix_11 ,"Muon_trackFitErrorMatrix_11");
    AddBranch(&Muon_trackFitErrorMatrix_12 ,"Muon_trackFitErrorMatrix_12");
    AddBranch(&Muon_trackFitErrorMatrix_22 ,"Muon_trackFitErrorMatrix_22");
    AddBranch(&Muon_dz_bt                  ,"Muon_dz_bt");
    AddBranch(&Muon_dxy_bt                 ,"Muon_dxy_bt");
  }
  //TTH
  if(_tthlepVar){
    AddBranch(&Muon_miniIsoRel        ,"Muon_miniIsoRel");
    AddBranch(&Muon_miniIsoCh         ,"Muon_miniIsoCh");
    AddBranch(&Muon_miniIsoNeu        ,"Muon_miniIsoNeu");
    AddBranch(&Muon_miniIsoPUsub      ,"Muon_miniIsoPUsub");
    AddBranch(&Muon_jetdr             ,"Muon_jetdr");
    AddBranch(&Muon_jetpt             ,"Muon_jetpt");
    AddBranch(&Muon_jetptratio        ,"Muon_jetptratio");
    AddBranch(&Muon_jetcsv            ,"Muon_jetcsv");
    AddBranch(&Muon_ptrel             ,"Muon_ptrel");
    AddBranch(&Muon_parel             ,"Muon_parel");
    AddBranch(&Muon_parelsub          ,"Muon_parelsub");
    AddBranch(&Muon_etarel            ,"Muon_etarel");
    AddBranch(&Muon_etarelsub         ,"Muon_etarelsub");
    AddBranch(&Muon_etapar            ,"Muon_etapar");
    AddBranch(&Muon_etaparsub         ,"Muon_etaparsub");
    AddBranch(&Muon_IP3Dsig_it        ,"Muon_IP3Dsig_it");
    AddBranch(&Muon_pvass             ,"Muon_pvass");
    AddBranch(&Muon_ptOVen            ,"Muon_ptOVen");
    AddBranch(&Muon_mujet_pfJetProbabilityBJetTag ,"Muon_mujet_pfJetProbabilityBJetTag");
    AddBranch(&Muon_mujet_pfCombinedMVABJetTags   ,"Muon_mujet_pfCombinedMVABJetTags");
    AddBranch(&Muon_mujet_qgl         ,"Muon_mujet_qgl");
    AddBranch(&Muon_mumass            ,"Muon_mumass");
    AddBranch(&Muon_mujet_mass        ,"Muon_mujet_mass");
    AddBranch(&Muon_mujet_Wmass       ,"Muon_mujet_Wmass");
    AddBranch(&Muon_mujet_Topmass     ,"Muon_mujet_Topmass");
    AddBranch(&Muon_mujet_WTopmass    ,"Muon_mujet_WTopmass");
    //Lep jet IP ntrks
    AddBranch(&Muon_IP3D_val           ,"Muon_IP3D_val");
    AddBranch(&Muon_IP3D_err           ,"Muon_IP3D_err");
    AddBranch(&Muon_IP3D_sig           ,"Muon_IP3D_sig");
    AddBranch(&Muon_IP2D_val           ,"Muon_IP2D_val");
    AddBranch(&Muon_IP2D_err           ,"Muon_IP2D_err");
    AddBranch(&Muon_IP2D_sig           ,"Muon_IP2D_sig");
    AddBranch(&Muon_sIP3D_val          ,"Muon_sIP3D_val");
    AddBranch(&Muon_sIP3D_err          ,"Muon_sIP3D_err");
    AddBranch(&Muon_sIP3D_sig          ,"Muon_sIP3D_sig");
    AddBranch(&Muon_sIP2D_val          ,"Muon_sIP2D_val");
    AddBranch(&Muon_sIP2D_err          ,"Muon_sIP2D_err");
    AddBranch(&Muon_sIP2D_sig          ,"Muon_sIP2D_sig");
    AddBranch(&Muon_IP1D_val           ,"Muon_IP1D_val");
    AddBranch(&Muon_IP1D_err           ,"Muon_IP1D_err");
    AddBranch(&Muon_IP1D_sig           ,"Muon_IP1D_sig");
    AddBranch(&Muon_sIP1D_val          ,"Muon_sIP1D_val");
    AddBranch(&Muon_sIP1D_err          ,"Muon_sIP1D_err");
    AddBranch(&Muon_sIP1D_sig          ,"Muon_sIP1D_sig");
    AddBranch(&Muon_lepjetMaxIP3D_val  ,"Muon_lepjetMaxIP3D_val");
    AddBranch(&Muon_lepjetMaxIP3D_sig  ,"Muon_lepjetMaxIP3D_sig");
    AddBranch(&Muon_lepjetMaxsIP3D_val ,"Muon_lepjetMaxsIP3D_val");
    AddBranch(&Muon_lepjetMaxsIP3D_sig ,"Muon_lepjetMaxsIP3D_sig");
    AddBranch(&Muon_lepjetMaxIP2D_val  ,"Muon_lepjetMaxIP2D_val");
    AddBranch(&Muon_lepjetMaxIP2D_sig  ,"Muon_lepjetMaxIP2D_sig");
    AddBranch(&Muon_lepjetMaxsIP2D_val ,"Muon_lepjetMaxsIP2D_val");
    AddBranch(&Muon_lepjetMaxsIP2D_sig ,"Muon_lepjetMaxsIP2D_sig");
    AddBranch(&Muon_lepjetMaxIP1D_val  ,"Muon_lepjetMaxIP1D_val");
    AddBranch(&Muon_lepjetMaxIP1D_sig  ,"Muon_lepjetMaxIP1D_sig");
    AddBranch(&Muon_lepjetMaxsIP1D_val ,"Muon_lepjetMaxsIP1D_val");
    AddBranch(&Muon_lepjetMaxsIP1D_sig ,"Muon_lepjetMaxsIP1D_sig");
    AddBranch(&Muon_lepjetAvIP3D_val   ,"Muon_lepjetAvIP3D_val");
    AddBranch(&Muon_lepjetAvIP3D_sig   ,"Muon_lepjetAvIP3D_sig");
    AddBranch(&Muon_lepjetAvsIP3D_val  ,"Muon_lepjetAvsIP3D_val");
    AddBranch(&Muon_lepjetAvsIP3D_sig  ,"Muon_lepjetAvsIP3D_sig");
    AddBranch(&Muon_lepjetAvIP2D_val   ,"Muon_lepjetAvIP2D_val");
    AddBranch(&Muon_lepjetAvIP2D_sig   ,"Muon_lepjetAvIP2D_sig");
    AddBranch(&Muon_lepjetAvsIP2D_val  ,"Muon_lepjetAvsIP2D_val");
    AddBranch(&Muon_lepjetAvsIP2D_sig  ,"Muon_lepjetAvsIP2D_sig");
    AddBranch(&Muon_lepjetAvIP1D_val   ,"Muon_lepjetAvIP1D_val");
    AddBranch(&Muon_lepjetAvIP1D_sig   ,"Muon_lepjetAvIP1D_sig");
    AddBranch(&Muon_lepjetAvsIP1D_val  ,"Muon_lepjetAvsIP1D_val");
    AddBranch(&Muon_lepjetAvsIP1D_sig  ,"Muon_lepjetAvsIP1D_sig");
    AddBranch(&Muon_lepjetchtrks       ,"Muon_lepjetchtrks");
    AddBranch(&Muon_lepjetpvchtrks     ,"Muon_lepjetpvchtrks");
    AddBranch(&Muon_lepjetnonpvchtrks  ,"Muon_lepjetnonpvchtrks");
    AddBranch(&Muon_lepjetndaus        ,"Muon_lepjetndaus");
    AddBranch(&Muon_lepjetpvchi2       ,"Muon_lepjetpvchi2");
    AddBranch(&Muon_lepjetnumno2trk    ,"Muon_lepjetnumno2trk");
  }
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
  Muon_pt_it.clear();
  Muon_ptErr_it.clear();
  Muon_pTErrOVpT_it.clear();
  Muon_pt_bt.clear();
  Muon_ptErr_bt.clear();
  Muon_pTErrOVpT_bt.clear();
  Muon_pt_tunePbt.clear();
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
  if(_AJVar){
    Muon_track_PCAx_bs.clear();
    Muon_track_PCAy_bs.clear();
    Muon_track_PCAz_bs.clear();
    Muon_track_PCAx_pv.clear();
    Muon_track_PCAy_pv.clear();
    Muon_track_PCAz_pv.clear();
    Muon_trackFitErrorMatrix_00.clear();
    Muon_trackFitErrorMatrix_01.clear();
    Muon_trackFitErrorMatrix_02.clear();
    Muon_trackFitErrorMatrix_11.clear();
    Muon_trackFitErrorMatrix_12.clear();
    Muon_trackFitErrorMatrix_22.clear();
    Muon_dz_bt.clear();
    Muon_dxy_bt.clear();
  }
  //TTH
  if(_tthlepVar){
    Muon_miniIsoRel.clear();
    Muon_miniIsoCh.clear();
    Muon_miniIsoNeu.clear();
    Muon_miniIsoPUsub.clear();
    Muon_jetdr.clear();
    Muon_jetpt.clear();
    Muon_jetptratio.clear();
    Muon_jetcsv.clear();
    Muon_ptrel.clear();
    Muon_parel.clear();
    Muon_parelsub.clear();
    Muon_etarel.clear();
    Muon_etarelsub.clear();
    Muon_etapar.clear();
    Muon_etaparsub.clear();
    Muon_IP3Dsig_it.clear();
    Muon_pvass.clear();
    Muon_ptOVen.clear();
    Muon_mujet_pfJetProbabilityBJetTag.clear();
    Muon_mujet_pfCombinedMVABJetTags.clear();
    Muon_mujet_qgl.clear();
    Muon_mumass.clear();
    Muon_mujet_mass.clear();
    Muon_mujet_Wmass.clear();
    Muon_mujet_Topmass.clear();
    Muon_mujet_WTopmass.clear();
    //Lep jet IP ntrks
    Muon_IP3D_val.clear();
    Muon_IP3D_err.clear();
    Muon_IP3D_sig.clear();
    Muon_IP2D_val.clear();
    Muon_IP2D_err.clear();
    Muon_IP2D_sig.clear();
    Muon_sIP3D_val.clear();
    Muon_sIP3D_err.clear();
    Muon_sIP3D_sig.clear();
    Muon_sIP2D_val.clear();
    Muon_sIP2D_err.clear();
    Muon_sIP2D_sig.clear();
    Muon_IP1D_val.clear();
    Muon_IP1D_err.clear();
    Muon_IP1D_sig.clear();
    Muon_sIP1D_val.clear();
    Muon_sIP1D_err.clear();
    Muon_sIP1D_sig.clear();
    Muon_lepjetMaxIP3D_val.clear();
    Muon_lepjetMaxIP3D_sig.clear();
    Muon_lepjetMaxsIP3D_val.clear();
    Muon_lepjetMaxsIP3D_sig.clear();
    Muon_lepjetMaxIP2D_val.clear();
    Muon_lepjetMaxIP2D_sig.clear();
    Muon_lepjetMaxsIP2D_val.clear();
    Muon_lepjetMaxsIP2D_sig.clear();
    Muon_lepjetMaxIP1D_val.clear();
    Muon_lepjetMaxIP1D_sig.clear();
    Muon_lepjetMaxsIP1D_val.clear();
    Muon_lepjetMaxsIP1D_sig.clear();
    Muon_lepjetAvIP3D_val.clear();
    Muon_lepjetAvIP3D_sig.clear();
    Muon_lepjetAvsIP3D_val.clear();
    Muon_lepjetAvsIP3D_sig.clear();
    Muon_lepjetAvIP2D_val.clear();
    Muon_lepjetAvIP2D_sig.clear();
    Muon_lepjetAvsIP2D_val.clear();
    Muon_lepjetAvsIP2D_sig.clear();
    Muon_lepjetAvIP1D_val.clear();
    Muon_lepjetAvIP1D_sig.clear();
    Muon_lepjetAvsIP1D_val.clear();
    Muon_lepjetAvsIP1D_sig.clear();
    Muon_lepjetchtrks.clear();
    Muon_lepjetpvchtrks.clear();
    Muon_lepjetnonpvchtrks.clear();
    Muon_lepjetndaus.clear();
    Muon_lepjetpvchi2.clear();
    Muon_lepjetnumno2trk.clear();
  }
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
}
bool MuonSelector::isGoodVertex(const reco::Vertex& vtx){
  if(vtx.isFake())                                   return false;
  if(vtx.ndof()<_vtx_ndof_min)                       return false;
  if(vtx.position().Rho()>_vtx_rho_max)              return false;
  if(fabs(vtx.position().Z()) > _vtx_position_z_max) return false;
  return true;
}
void MuonSelector::get_muminiIso_info(const pat::PackedCandidateCollection& pcc, double rho, const pat::Muon& cand, double& miniIso, double& miniIsoCh, double& miniIsoNeu, double& miniIsoPUsub){
  double miniIsoConeSize = 10.0/min(max(cand.pt(), 50.),200.);
  vector<const pat::PackedCandidate *> pfc_all; pfc_all.clear();
  vector<const pat::PackedCandidate *> pfc_ch;  pfc_ch.clear();
  vector<const pat::PackedCandidate *> pfc_neu; pfc_neu.clear();
  vector<const pat::PackedCandidate *> pfc_pu;  pfc_pu.clear();
  get_chneupu_pcc(pcc,pfc_all,pfc_ch,pfc_neu,pfc_pu);
  miniIsoCh  = get_isosumraw(pfc_ch,  cand, miniIsoConeSize, 0.0001, 0, SelfVetoPolicyMu::selfVetoAll, 0);
  miniIsoNeu = get_isosumraw(pfc_neu, cand, miniIsoConeSize, 0.01, 0.5, SelfVetoPolicyMu::selfVetoAll, 0);
  double effarea    = get_effarea(cand.eta());
  double correction = rho*effarea*pow((miniIsoConeSize/0.3),2);
  miniIsoPUsub = std::max(0.0, miniIsoNeu-correction);
  miniIso = miniIsoCh+miniIsoPUsub;
}
void MuonSelector::get_chneupu_pcc(const pat::PackedCandidateCollection& pcc,vector<const pat::PackedCandidate *>& pfc_all,vector<const pat::PackedCandidate *>& pfc_ch,vector<const pat::PackedCandidate *>& pfc_neu,vector<const pat::PackedCandidate *>&pfc_pu){
  for(const pat::PackedCandidate &p : pcc){
    pfc_all.push_back(&p);
    if(p.charge()==0){
      pfc_neu.push_back(&p);
    }else{
      if((abs(p.pdgId())==211)){// || ((abs(p.pdgId()) == 11 ) || (abs(p.pdgId()) == 13 ))){
        if(p.fromPV()>1 && fabs(p.dz())<9999){
          pfc_ch.push_back(&p);
        }else{
          pfc_pu.push_back(&p);
        }
      }
    }
  }
  std::sort(pfc_ch.begin(),  pfc_ch.end(),  ByEta());
  std::sort(pfc_neu.begin(), pfc_neu.end(), ByEta());
  std::sort(pfc_pu.begin(),  pfc_pu.end(),  ByEta());
}
double MuonSelector::get_isosumraw(const std::vector<const pat::PackedCandidate *> & pcc, const pat::Muon& cand, double miniIsoConeSize, double innerR, double ptTh, SelfVetoPolicyMu::SelfVetoPolicyMu selfVeto, int pdgId){
  //Look for cand sources
  std::vector<const reco::Candidate *> vetos; vetos.clear();
  for(uint i=0, n=cand.numberOfSourceCandidatePtrs(); i<n; ++i){
    if(selfVeto == SelfVetoPolicyMu::selfVetoNone) break;
    const reco::CandidatePtr &cp = cand.sourceCandidatePtr(i);
    if(cp.isNonnull() && cp.isAvailable()){
      vetos.push_back(&*cp);
      if(selfVeto == SelfVetoPolicyMu::selfVetoFirst) break;
    }
  }
  //Get the isolation
  double isosum = 0;
  float miniIsoConeSize2 = miniIsoConeSize*miniIsoConeSize, innerR2 = innerR*innerR;
  typedef std::vector<const pat::PackedCandidate *>::const_iterator IT;
  IT pccbegin = std::lower_bound(pcc.begin(), pcc.end(), cand.eta() - miniIsoConeSize, ByEta());
  IT pccend   = std::upper_bound(pccbegin,    pcc.end(), cand.eta() + miniIsoConeSize, ByEta());
  for(IT pc = pccbegin; pc<pccend; ++pc){
    //pdgId veto
    if(pdgId>0 && abs((*pc)->pdgId())!=pdgId) continue;
    //pT requirement 
    if(ptTh>0 && (*pc)->pt()<ptTh) continue;
    //cone region
    double dr2 = reco::deltaR2(**pc, cand);
    if(miniIsoConeSize2<dr2 || dr2<innerR2) continue;
    //itself veto
    if(std::find(vetos.begin(), vetos.end(),*pc)!=vetos.end()) continue;
    //add to sum
    isosum += (*pc)->pt();
  }
  return isosum;
}
double MuonSelector::get_effarea(double eta){
  double effarea = -1;
  if(abs(eta) < 0.8)      effarea = 0.0735;
  else if(abs(eta) < 1.3) effarea = 0.0619;
  else if(abs(eta) < 2.0) effarea = 0.0465;
  else if(abs(eta) < 2.2) effarea = 0.0433;
  else                    effarea = 0.0577;
  return effarea;
}
void MuonSelector::get_mujet_info(const pat::Muon& mu, const edm::Event& iEvent, const edm::EventSetup& iSetup,
                                  double& mujet_mindr, double& mujet_pt, double& muptOVmujetpt,
                                  double& mujet_pfCombinedInclusiveSecondaryVertexV2BJetTags, double& mujet_pfJetProbabilityBJetTags, double& mujet_pfCombinedMVABJetTags, double& mujet_qgl,
                                  double& jx, double& jy, double& jz, double& muptrel, double& muparel, double& muparelsub, int& lepjetidx){
  //Look for jet associated to mu
  edm::Handle<pat::JetCollection> jets;
  iEvent.getByToken(jets_, jets);
  pat::Jet mujet;
  int currjetpos = 0;
  for(const pat::Jet &j : *jets){
    pat::Jet jet = j;
    double dr = deltaR(mu.p4(),jet.p4());
    if(dr<mujet_mindr){
      mujet_mindr = dr;
      mujet       = jet;
      lepjetidx = currjetpos;
    }
    currjetpos++;
  }
  //Get the other info
  if(mujet.jecSetsAvailable()){
    double L2L3_corr = mujet.p4().E()/mujet.correctedJet(1).p4().E(); 
    mujet.setP4(((mujet.correctedJet(1).p4()-mu.p4())*L2L3_corr)+mu.p4());
  }
  mujet_pt       = mujet.pt();
  muptOVmujetpt  = min(mu.pt()/mujet.pt(), 1.5);
  mujet_pfCombinedInclusiveSecondaryVertexV2BJetTags = max(double(mujet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags")), 0.0);
  if(mujet_pfCombinedInclusiveSecondaryVertexV2BJetTags!=mujet_pfCombinedInclusiveSecondaryVertexV2BJetTags) mujet_pfCombinedInclusiveSecondaryVertexV2BJetTags = -996;
  mujet_pfJetProbabilityBJetTags = mujet.bDiscriminator("pfJetProbabilityBJetTags");
  if(mujet_pfJetProbabilityBJetTags!=mujet_pfJetProbabilityBJetTags) mujet_pfJetProbabilityBJetTags = -996;
  mujet_pfCombinedMVABJetTags    = mujet.bDiscriminator("pfCombinedMVABJetTags"); 
  if(mujet_pfCombinedMVABJetTags!=mujet_pfCombinedMVABJetTags) mujet_pfCombinedMVABJetTags = -996;
  if(_qglVar){
    edm::Handle<edm::View<pat::Jet>> jets_QGL;
    iEvent.getByToken(jetsToken,jets_QGL);
    edm::Handle<edm::ValueMap<float>> qgHandle;
    iEvent.getByToken(qgToken_, qgHandle);
    for(auto jet = jets_QGL->begin();  jet != jets_QGL->end(); ++jet){
      if(distance(jets_QGL->begin(),jet)!=lepjetidx) continue;
      edm::RefToBase<pat::Jet> jetRef(edm::Ref<edm::View<pat::Jet> >(jets_QGL, jet - jets_QGL->begin()));
      mujet_qgl = (*qgHandle)[jetRef];
      break;
    }
  }
  jx = mujet.px();
  jy = mujet.py();
  jz = mujet.pz();
  TLorentzVector mu_lv    = TLorentzVector(mu.px(),mu.py(),mu.pz(),mu.p4().E());
  TLorentzVector mujet_lv = TLorentzVector(mujet.px(),mujet.py(),mujet.pz(),mujet.p4().E());
  muptrel = mu_lv.Perp((mujet_lv-mu_lv).Vect());
  muparel = (mu_lv.Vect()).Dot(mujet_lv.Vect().Unit());
  muparelsub = (mu_lv.Vect()).Dot((mujet_lv-mu_lv).Vect().Unit());
}
int MuonSelector::pvassociation(const pat::Muon& mu, const pat::PackedCandidateCollection& pcc){
  int pvass = -1;
  double mindr = 0.3;
  for(const pat::PackedCandidate &cpf : pcc){ 
    if(deltaR(mu.p4(),cpf.p4())<mindr             //dR is the standard geometrical way to associate 
       && (fabs(mu.pt()-cpf.pt())/mu.pt())<0.05   //Check in pT, because ele,tau are usually faked by jets (many trks) and dR may not be enough
       && cpf.charge()!=0 && cpf.numberOfHits()>0 //Leptons are charged and built from tracks, also to be consistent with PV tracks  
    ){
      mindr = deltaR(mu.p4(),cpf.p4());
      pvass = cpf.fromPV();
    }
  }
  return pvass;  
}
double MuonSelector::relEta(const math::XYZVector& dir, const math::XYZVector& track){
 double etarel = -9999;
 double momPar = track.Dot(dir.Unit());
 double energy = std::sqrt(track.Mag2() + pow(0.105658,2));
 if((energy - momPar)!=0 && (energy + momPar)!=0) etarel = 0.5 * log((energy + momPar) / (energy - momPar));
 return etarel;
//Old def
//const math::XYZVector& vector, const math::XYZVector& axis){
//  double etarel = 15; //Take this as a default value and in the end use min(etarel,15)
//  double mag = vector.r() * axis.r();
//  double dot = vector.Dot(axis);
//  if((mag-dot)!=0 && (mag+dot)!=0) etarel = -log((mag-dot)/(mag+dot)) / 2;
//  return etarel;  
}
double MuonSelector::parEta(const math::XYZVector& dir, const math::XYZVector& track){
 double etapar  = -9999;
 TVector3 trkdir(track.x(),track.y(),track.z());
 TVector3 dirdir(dir.x(),dir.y(),dir.z());
 double momPerp = trkdir.Perp(dirdir);
 double energy  = std::sqrt(track.Mag2() + pow(0.105658,2));
 if((energy - momPerp)!=0 && (energy + momPerp)!=0) etapar = 0.5 * log((energy + momPerp) / (energy - momPerp));
 return etapar;
//Old def
//const math::XYZVector& vector, const math::XYZVector& axis){
//  double etarel = 15; //Take this as a default value and in the end use min(etarel,15)
//  double mag = vector.r() * axis.r();
//  double dot = vector.Dot(axis);
//  if((mag-dot)!=0 && (mag+dot)!=0) etarel = -log((mag-dot)/(mag+dot)) / 2;
//  return etarel;  
}
double MuonSelector::get_lepWmass(const pat::Muon& mu, const edm::Event& iEvent, int& lepjetidx){
  double lepWmass = 0;
  double minchi2  = 999;
  //Take lepjet lv 
  edm::Handle<pat::JetCollection> jets;
  iEvent.getByToken(jets_, jets);
  TLorentzVector lepjet;
  if(lepjetidx!=-1){
    const pat::Jet & jet = (*jets)[lepjetidx];
    lepjet.SetPtEtaPhiE(jet.pt(),jet.eta(),jet.phi(),jet.energy()); 
  }else{
    lepjet.SetPtEtaPhiE(mu.pt(),mu.eta(),mu.phi(),mu.energy());
  }
  //Start combinatorial with other jets
  for(const pat::Jet &j : *jets){
    pat::Jet jet = j;
    if(!(j.pt()>25 && fabs(j.eta())<2.4 && deltaR(j.p4(),mu.p4())>0.4)) continue;
    if(!(j.neutralHadronEnergyFraction()<0.99 && j.neutralEmEnergyFraction()<0.99 && (j.chargedMultiplicity() + j.neutralMultiplicity())>1
       && j.chargedHadronEnergyFraction()>0 && j.chargedEmEnergyFraction()<0.99 && j.chargedMultiplicity()>0
       && j.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags")<0.89)//0.605)
      ) continue;
    TLorentzVector jet_lv = TLorentzVector(j.px(),j.py(),j.pz(),j.p4().E());
    double currmass = (lepjet+jet_lv).M();
    double currchi2 = pow((currmass-80.385)/25,2); 
    if(currchi2<minchi2){
      minchi2 = currchi2;
      lepWmass = currmass;
    }
  }
  return lepWmass;
};
double MuonSelector::get_lepTopmass(const pat::Muon& mu, const edm::Event& iEvent, int& lepjetidx){
  double lepTopmass = 0;
  double minchi2  = 999;
  //Take lepjet lv 
  edm::Handle<pat::JetCollection> jets;
  iEvent.getByToken(jets_, jets);
  TLorentzVector lepjet;
  if(lepjetidx!=-1){
    const pat::Jet & jet = (*jets)[lepjetidx];
    lepjet.SetPtEtaPhiE(jet.pt(),jet.eta(),jet.phi(),jet.energy()); 
  }else{
    lepjet.SetPtEtaPhiE(mu.pt(),mu.eta(),mu.phi(),mu.energy());
  }
  //Start combinatorial with other jets (you can assume any type of jets. Instead asking non b jet would favor the case where lepjet is a b)
  for(uint j1 = 0; j1 < jets->size(); j1++){
    const pat::Jet & jet1 = (*jets)[j1];
    if(!(jet1.pt()>25 && fabs(jet1.eta())<2.4 && deltaR(jet1.p4(),mu.p4())>0.4)) continue;
    if(!(jet1.neutralHadronEnergyFraction()<0.99 && jet1.neutralEmEnergyFraction()<0.99 && (jet1.chargedMultiplicity() + jet1.neutralMultiplicity())>1
       && jet1.chargedHadronEnergyFraction()>0 && jet1.chargedEmEnergyFraction()<0.99 && jet1.chargedMultiplicity()>0
       && jet1.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags")<0.89
        )
      ) continue;
    for(uint j2 = j1+1; j2 < jets->size(); j2++){
      const pat::Jet & jet2 = (*jets)[j2];
      if(!(jet2.pt()>25 && fabs(jet2.eta())<2.4 && deltaR(jet2.p4(),mu.p4())>0.4)) continue;
      if(!(jet2.neutralHadronEnergyFraction()<0.99 && jet2.neutralEmEnergyFraction()<0.99 && (jet2.chargedMultiplicity() + jet2.neutralMultiplicity())>1
         && jet2.chargedHadronEnergyFraction()>0 && jet2.chargedEmEnergyFraction()<0.99 && jet2.chargedMultiplicity()>0
         && jet2.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags")<0.89
          )
        ) continue;
      TLorentzVector jet1_lv = TLorentzVector(jet1.px(),jet1.py(),jet1.pz(),jet1.p4().E());      
      TLorentzVector jet2_lv = TLorentzVector(jet2.px(),jet2.py(),jet2.pz(),jet2.p4().E());      
      TLorentzVector W_lv   = jet1_lv+jet2_lv;
      TLorentzVector Top_lv = W_lv+lepjet;    
      double currmassw   = W_lv.M();  
      double currmasstop = Top_lv.M(); 
      double currchi2 = pow((currmassw-80.385)/25,2) + pow((currmasstop-173.34)/50,2);
      if(currchi2<minchi2){
        minchi2 = currchi2;
        lepTopmass = currmasstop;
      }    
    }
  }
  return lepTopmass;
};
double MuonSelector::get_lepWTopmass(const pat::Muon& mu, const edm::Event& iEvent, int& lepjetidx){
  double lepWTopmass = 0;
  double minchi2  = 999;
  //Take lepjet lv 
  edm::Handle<pat::JetCollection> jets;
  iEvent.getByToken(jets_, jets);
  TLorentzVector lepjet;
  if(lepjetidx!=-1){
    const pat::Jet & jet = (*jets)[lepjetidx];
    lepjet.SetPtEtaPhiE(jet.pt(),jet.eta(),jet.phi(),jet.energy()); 
  }else{
    lepjet.SetPtEtaPhiE(mu.pt(),mu.eta(),mu.phi(),mu.energy());
  }
  //Start combinatorial with other jets (you can assume any type of jets. Instead asking non b jet would favor the case where lepjet is a b)
  for(uint j1 = 0; j1 < jets->size(); j1++){
    const pat::Jet & jet1 = (*jets)[j1];
    if(!(jet1.pt()>25 && fabs(jet1.eta())<2.4 && deltaR(jet1.p4(),mu.p4())>0.4)) continue;
    if(!(jet1.neutralHadronEnergyFraction()<0.99 && jet1.neutralEmEnergyFraction()<0.99 && (jet1.chargedMultiplicity() + jet1.neutralMultiplicity())>1
       && jet1.chargedHadronEnergyFraction()>0 && jet1.chargedEmEnergyFraction()<0.99 && jet1.chargedMultiplicity()>0
       && jet1.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags")<0.89
        )
      ) continue;
    TLorentzVector jet1_lv = TLorentzVector(jet1.px(),jet1.py(),jet1.pz(),jet1.p4().E());         
    TLorentzVector W_lv    = jet1_lv+lepjet;
    for(uint j2 = j1+1; j2 < jets->size(); j2++){
      const pat::Jet & jet2 = (*jets)[j2];
      if(!(jet2.pt()>25 && fabs(jet2.eta())<2.4 && deltaR(jet2.p4(),mu.p4())>0.4)) continue;
      if(!(jet2.neutralHadronEnergyFraction()<0.99 && jet2.neutralEmEnergyFraction()<0.99 && (jet2.chargedMultiplicity() + jet2.neutralMultiplicity())>1
         && jet2.chargedHadronEnergyFraction()>0 && jet2.chargedEmEnergyFraction()<0.99 && jet2.chargedMultiplicity()>0
         && jet2.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags")>0.605
          )
        ) continue;
      TLorentzVector jet2_lv = TLorentzVector(jet2.px(),jet2.py(),jet2.pz(),jet2.p4().E());      
      TLorentzVector Top_lv = W_lv+jet2_lv;    
      double currmassw   = W_lv.M();  
      double currmasstop = Top_lv.M(); 
      double currchi2 = pow((currmassw-80.385)/25,2) + pow((currmasstop-173.34)/50,2);
      if(currchi2<minchi2){
        minchi2 = currchi2;
        lepWTopmass = currmassw;
      }    
    }
  }
  return lepWTopmass;
};
void MuonSelector::IP3D2D(TransientTrack ttrk, const reco::Vertex& vtx, GlobalVector gv, double& IP3D_val,double& IP3D_err,double& IP3D_sig, double& sIP3D_val,double& sIP3D_err,double& sIP3D_sig, double& IP2D_val,double& IP2D_err,double& IP2D_sig, double& sIP2D_val,double& sIP2D_err,double& sIP2D_sig){
 pair<bool, Measurement1D> currIP;
 //3D
 currIP = IPTools::absoluteImpactParameter3D(ttrk,vtx);
 if(currIP.first){
   if(currIP.second.value()==currIP.second.value())               IP3D_val = currIP.second.value();
   if(currIP.second.error()==currIP.second.error())               IP3D_err = currIP.second.error();
   if(currIP.second.significance()==currIP.second.significance()) IP3D_sig = currIP.second.significance();
 }
 //s3D
 currIP = IPTools::signedImpactParameter3D(ttrk,gv,vtx);
 if(currIP.first){
   if(currIP.second.value()==currIP.second.value())               sIP3D_val = currIP.second.value();
   if(currIP.second.error()==currIP.second.error())               sIP3D_err = currIP.second.error();
   if(currIP.second.significance()==currIP.second.significance()) sIP3D_sig = currIP.second.significance();
 }
 //2D 
 currIP = IPTools::absoluteTransverseImpactParameter(ttrk,vtx);
 if(currIP.first){
   if(currIP.second.value()==currIP.second.value())               IP2D_val = currIP.second.value();
   if(currIP.second.error()==currIP.second.error())               IP2D_err = currIP.second.error();
   if(currIP.second.significance()==currIP.second.significance()) IP2D_sig = currIP.second.significance();
 }
 //s2D 
 currIP = IPTools::signedTransverseImpactParameter(ttrk,gv,vtx);
 if(currIP.first){
   if(currIP.second.value()==currIP.second.value())               sIP2D_val = currIP.second.value();
   if(currIP.second.error()==currIP.second.error())               sIP2D_err = currIP.second.error();
   if(currIP.second.significance()==currIP.second.significance()) sIP2D_sig = currIP.second.significance();
 }
}
void MuonSelector::zIP1D(TransientTrack ttrk, const reco::Vertex& vtx, GlobalVector gv, double& IP1D_val,double& IP1D_err,double& IP1D_sig, double& sIP1D_val,double& sIP1D_err,double& sIP1D_sig){
 SignedTransverseImpactParameter stip;
 pair<bool, Measurement1D> currIP = stip.zImpactParameter(ttrk,gv,vtx);
 if(currIP.first){
   if(currIP.second.value()==currIP.second.value())               IP1D_val = fabs(currIP.second.value());
   if(currIP.second.error()==currIP.second.error())               IP1D_err = fabs(currIP.second.error());
   if(currIP.second.significance()==currIP.second.significance()) IP1D_sig = fabs(currIP.second.significance());
   if(currIP.second.value()==currIP.second.value())               sIP1D_val = currIP.second.value();
   if(currIP.second.error()==currIP.second.error())               sIP1D_err = currIP.second.error();
   if(currIP.second.significance()==currIP.second.significance()) sIP1D_sig = currIP.second.significance();
 }
}
void MuonSelector::lepjetIP(const pat::Jet& jet, const reco::Vertex& vtx, GlobalVector lepjetgv, const TransientTrackBuilder& ttrkbuilder,
                    double& lepjetMaxIP3D_val, double& lepjetMaxIP3D_sig, double& lepjetMaxsIP3D_val, double& lepjetMaxsIP3D_sig, double& lepjetMaxIP2D_val, double& lepjetMaxIP2D_sig, double& lepjetMaxsIP2D_val, double& lepjetMaxsIP2D_sig, double& lepjetMaxIP1D_val, double& lepjetMaxIP1D_sig, double& lepjetMaxsIP1D_val, double& lepjetMaxsIP1D_sig,
                    double& lepjetAvIP3D_val, double& lepjetAvIP3D_sig, double& lepjetAvsIP3D_val, double& lepjetAvsIP3D_sig, double& lepjetAvIP2D_val, double& lepjetAvIP2D_sig, double& lepjetAvsIP2D_val, double& lepjetAvsIP2D_sig, double& lepjetAvIP1D_val, double& lepjetAvIP1D_sig, double& lepjetAvsIP1D_val, double& lepjetAvsIP1D_sig,
                    double& denlepjetAvIP3D_val, double& denlepjetAvIP3D_sig, double& denlepjetAvsIP3D_val, double& denlepjetAvsIP3D_sig, double& denlepjetAvIP2D_val, double& denlepjetAvIP2D_sig, double& denlepjetAvsIP2D_val, double& denlepjetAvsIP2D_sig, double& denlepjetAvIP1D_val, double& denlepjetAvIP1D_sig, double& denlepjetAvsIP1D_val, double& denlepjetAvsIP1D_sig,
                    double& Lep_IP3D_val
){
 //Access jet daughters
 vector<CandidatePtr> jdaus(jet.daughterPtrVector());
 sort(jdaus.begin(), jdaus.end(), [](const reco::CandidatePtr &p1, const reco::CandidatePtr &p2) {return p1->pt() > p2->pt();});
 for(uint jd=0; jd<jdaus.size(); jd++){
  const pat::PackedCandidate &jcand = dynamic_cast<const pat::PackedCandidate &>(*jdaus[jd]);
  if(deltaR(jcand.p4(),jet.p4())>0.4) continue;
  Track trk = Track(jcand.pseudoTrack());
  bool isgoodtrk = is_goodtrk(trk,vtx);
  //Minimal conditions for a track 
  if(isgoodtrk && jcand.charge()!=0 && jcand.fromPV()>1){
    TransientTrack ttrk = ttrkbuilder.build(&trk);
    //Current IP values 
    double IP3D_val  = -9999;
    double IP3D_err  = -9999;
    double IP3D_sig  = -9999;
    double sIP3D_val = -9999;
    double sIP3D_err = -9999;
    double sIP3D_sig = -9999;
    double IP2D_val  = -9999;
    double IP2D_err  = -9999;
    double IP2D_sig  = -9999;
    double sIP2D_val = -9999;
    double sIP2D_err = -9999;
    double sIP2D_sig = -9999;
    double IP1D_val  = -9999;
    double IP1D_err  = -9999;
    double IP1D_sig  = -9999;
    double sIP1D_val = -9999;
    double sIP1D_err = -9999;
    double sIP1D_sig = -9999;
    IP3D2D(ttrk,vtx,lepjetgv,IP3D_val,IP3D_err,IP3D_sig,sIP3D_val,sIP3D_err,sIP3D_sig,IP2D_val,IP2D_err,IP2D_sig,sIP2D_val,sIP2D_err,sIP2D_sig);
    zIP1D(ttrk,vtx,lepjetgv,IP1D_val,IP1D_err,IP1D_sig,sIP1D_val,sIP1D_err,sIP1D_sig);
    //Max Lep jet IP
    if(IP3D_val>lepjetMaxIP3D_val)   lepjetMaxIP3D_val  = IP3D_val;
    if(IP3D_sig>lepjetMaxIP3D_sig)   lepjetMaxIP3D_sig  = IP3D_sig;
    if(sIP3D_val>lepjetMaxsIP3D_val) lepjetMaxsIP3D_val = sIP3D_val;
    if(sIP3D_sig>lepjetMaxsIP3D_sig) lepjetMaxsIP3D_sig = sIP3D_sig;
    if(IP2D_val>lepjetMaxIP2D_val)   lepjetMaxIP2D_val  = IP2D_val;
    if(IP2D_sig>lepjetMaxIP2D_sig)   lepjetMaxIP2D_sig  = IP2D_sig;
    if(sIP2D_val>lepjetMaxsIP2D_val) lepjetMaxsIP2D_val = sIP2D_val;
    if(sIP2D_sig>lepjetMaxsIP2D_sig) lepjetMaxsIP2D_sig = sIP2D_sig;
    if(IP1D_val>lepjetMaxIP1D_val)   lepjetMaxIP1D_val  = IP1D_val;
    if(IP1D_sig>lepjetMaxIP1D_sig)   lepjetMaxIP1D_sig  = IP1D_sig;
    if(sIP1D_val>lepjetMaxsIP1D_val) lepjetMaxsIP1D_val = sIP1D_val;
    if(sIP1D_sig>lepjetMaxsIP1D_sig) lepjetMaxsIP1D_sig = sIP1D_sig;
    //Max Lep jet IP
    if((fabs(IP3D_val-Lep_IP3D_val)/IP3D_val)>0.01){
      if(IP3D_val!=-9999 ){
        lepjetAvIP3D_val     += IP3D_val;
        denlepjetAvIP3D_val  += 1;
      }
      if(IP3D_sig!=-9999 ){
        lepjetAvIP3D_sig     += IP3D_sig;
        denlepjetAvIP3D_sig  += 1;
      }
      if(sIP3D_val!=-9999){
        lepjetAvsIP3D_val    += sIP3D_val;
        denlepjetAvsIP3D_val += 1;
      }
      if(sIP3D_sig!=-9999){
        lepjetAvsIP3D_sig    += sIP3D_sig;
        denlepjetAvsIP3D_sig += 1;
      }
      if(IP2D_val!=-9999 ){
        lepjetAvIP2D_val     += IP2D_val;
        denlepjetAvIP2D_val  += 1;
      }
      if(IP2D_sig!=-9999 ){
        lepjetAvIP2D_sig     += IP2D_sig;
        denlepjetAvIP2D_sig  += 1;
      }
      if(sIP2D_val!=-9999){
        lepjetAvsIP2D_val    += sIP2D_val;
        denlepjetAvsIP2D_val += 1;
      }
      if(sIP2D_sig!=-9999){
        lepjetAvsIP2D_sig    += sIP2D_sig;
        denlepjetAvsIP2D_sig += 1;
      }
      if(IP1D_val!=-9999 ){
        lepjetAvIP1D_val     += IP1D_val;
        denlepjetAvIP1D_val  += 1;
      }
      if(IP1D_sig!=-9999 ){
        lepjetAvIP1D_sig     += IP1D_sig;
        denlepjetAvIP1D_sig  += 1;
      }
      if(sIP1D_val!=-9999){
        lepjetAvsIP1D_val    += sIP1D_val;
        denlepjetAvsIP1D_val += 1;
      }
      if(sIP1D_sig!=-9999){
        lepjetAvsIP1D_sig    += sIP1D_sig;
        denlepjetAvsIP1D_sig += 1;
      }
    }
  }//Ch trks 
 }//Loop on jet daus
}
bool MuonSelector::is_goodtrk(Track trk,const reco::Vertex& vtx){
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
void MuonSelector::lepjetTrks(const pat::Jet& jet, const reco::Vertex& vtx, double& lepjetchtrks, double& lepjetpvchtrks, double& lepjetnonpvchtrks, double& lepjetndaus){
 //Access jet daughters
 vector<CandidatePtr> jdaus(jet.daughterPtrVector());
 sort(jdaus.begin(), jdaus.end(), [](const reco::CandidatePtr &p1, const reco::CandidatePtr &p2) {return p1->pt() > p2->pt();});
 for(uint jd=0; jd<jdaus.size(); jd++){
  const pat::PackedCandidate &jcand = dynamic_cast<const pat::PackedCandidate &>(*jdaus[jd]);
  if(deltaR(jcand.p4(),jet.p4())>0.4) continue;
  Track trk = Track(jcand.pseudoTrack());
  bool isgoodtrk = is_goodtrk(trk,vtx);
  //Minimal conditions for a track 
  if(isgoodtrk && jcand.charge()!=0 && jcand.fromPV()>1){
    //Get jet trk num
    lepjetchtrks++;
    if(jcand.fromPV()==pat::PackedCandidate::PVUsedInFit){
      lepjetpvchtrks++;
    }else{
      lepjetnonpvchtrks++;
    }
  }//Ch trks 
 }//Loop on jet daus
 lepjetndaus = jdaus.size();
}
void MuonSelector::lepjetVtxCompatibility(const pat::Jet& jet, const reco::Vertex& vtx, const TransientTrackBuilder& ttrkbuilder, double& lepjetpvchi2, double& lepjetnumno2tr){
 //Access jet daughters
 vector<TransientTrack> jetttrks;
 vector<CandidatePtr> jdaus(jet.daughterPtrVector());
 sort(jdaus.begin(), jdaus.end(), [](const reco::CandidatePtr &p1, const reco::CandidatePtr &p2) {return p1->pt() > p2->pt();});
 for(uint jd=0; jd<jdaus.size(); jd++){
  const pat::PackedCandidate &jcand = dynamic_cast<const pat::PackedCandidate &>(*jdaus[jd]);
  if(deltaR(jcand.p4(),jet.p4())>0.4) continue;
  Track trk = Track(jcand.pseudoTrack());
  bool isgoodtrk = is_goodtrk(trk,vtx);
  //Minimal conditions for a track 
  if(isgoodtrk && jcand.charge()!=0 && jcand.fromPV()>1){
    TransientTrack ttrk = ttrkbuilder.build(&trk);
    jetttrks.push_back(ttrk);
  }//Ch trks 
 }//Loop on jet daus
 //LepJet vertex chi2
 TransientVertex tv;
 if(jetttrks.size()>=2) tv = vertexfittermu.vertex(jetttrks);
 if(tv.isValid()) lepjetpvchi2 = tv.totalChiSquared()/tv.degreesOfFreedom();
 double num2v = 0; double numno2v = 0; 
 get_2trksinfo(jetttrks,  num2v,  numno2v);
 if((numno2v+num2v)!=0) lepjetnumno2tr = numno2v/(numno2v+num2v); 
 else                   lepjetnumno2tr = 0;
}
void MuonSelector::get_2trksinfo(vector<TransientTrack> ttrks, double& num2v, double& numno2v){
 for(uint t=0; t<ttrks.size(); t++){
  for(uint t2=t+1; t2<ttrks.size(); t2++){
   vector<TransientTrack> twotrks;
   twotrks.push_back(ttrks[t]);
   twotrks.push_back(ttrks[t2]);
   TransientVertex tv;
   if(ttrks.size()>=2) tv = vertexfittermu.vertex(ttrks);
   if(tv.isValid() && TMath::Prob(tv.totalChiSquared(),tv.degreesOfFreedom())>0.05){
    num2v++;
   }else{
    numno2v++;
   }
  }
 }
}
//TTHLep synch
  //const JetCorrector* corrector = JetCorrector::getJetCorrector( "ak4PFCHSL1L2L3Residual", iSetup );
  //const JetCorrector* corrector = JetCorrector::getJetCorrector( "ak4PFchsL1L2L3", iSetup );
  //  pat::Jet jet = j;//j.correctedJet(0);
    //double scale = corrector->correction(jet, iEvent, iSetup);
    //jet.scaleEnergy(scale);
  //Some info
  //cout<<"Jet info "<<mujet<<endl;
  //cout<<"Corrected (L0)"<<setw(20)<<mujet.correctedJet(0).p4().E()<<endl; 
  //cout<<"Corrected (L1)"<<setw(20)<<mujet.correctedJet(1).p4().E()<<endl;
  //cout<<"Corrected (L2)"<<setw(20)<<mujet.correctedJet(2).p4().E()<<endl;
  //cout<<"Corrected (L3)"<<setw(20)<<mujet.correctedJet(3).p4().E()<<endl;
  //cout<<"Corrected Fina"<<setw(20)<<mujet.p4().E()<<endl; 
//GenPart
        //cout<<setw(20)<<"pT"<<setw(20)<<"eta"<<setw(20)<<"phi"<<setw(20)<<"energy"<<endl;
        //cout<<setw(20)<<genpart->pt()<<setw(20)<<genpart->eta()<<setw(20)<<genpart->phi()<<setw(20)<<genpart->energy()<<endl;
        //cout<<setw(20)<<"isPrompt"<<setw(20)<<"isfromtau"<<setw(20)<<"pdgId"<<endl;
        //cout<<setw(20)<<genpart->isPromptFinalState()<<setw(20)<<genpart->isDirectPromptTauDecayProductFinalState()<<setw(20)<<genpart->pdgId()<<endl;
//Synch TTHLep
      //Print info
      //cout<<setiosflags(ios::fixed)<<setprecision(5);
      /*
      if(!amu &&
        mu->innerTrack().isNonnull()
        && mu->pt()>5 && fabs(mu->eta())<2.4
        && fabs(mu->innerTrack()->dxy(firstGoodVertex.position()))<=0.05 && fabs(mu->innerTrack()->dz(firstGoodVertex.position()))<0.1
        && miniIso/mu->pt()<0.4 && fabs(mu->dB(pat::Muon::PV3D))/mu->edB(pat::Muon::PV3D)<8 && mu->isLooseMuon()
      )
      //if(1)
      {
        //cout<<setw(10)<<"event"<<setw(10)<<"pT"<<setw(10)<<"Eta"<<setw(10)<<"Phi"<<setw(10)<<"E"<<setw(5)<<"pdgID"<<setw(5)<<"charge"<<setw(15)<<"miniIso"<<setw(15)<<"miniIsoCharged"<<setw(15)<<"miniIsoNeutral"<<setw(10)<<"jetPtRel"<<setw(10)<<"jetCSV"<<setw(10)<<"jetPtRatio"<<setw(10)<<"sip3D"<<setw(10)<<"dxy"<<setw(10)<<"dz"<<setw(21)<<"segmentCompatibility"<<endl;
	cout<<"Mu"<<setw(10)<<iEvent.id().event()<<setw(10)<<mu->pt()<<setw(10)<<mu->eta()<<setw(10)<<mu->phi()<<setw(10)<<mu->energy()<<setw(5)<<mu->pdgId()<<setw(5)<<mu->charge()<<setw(15)<<miniIso/mu->pt()<<setw(15)<<miniIsoCh<<setw(15)<<miniIsoPUsub<<setw(10)<<muptrel<<setw(10)<<mujet_pfCombinedInclusiveSecondaryVertexV2BJetTags<<setw(10)<<mu->pt()/mujet_pt;
	if(mu->innerTrack().isNonnull()){
	cout<<setw(10)<<fabs(mu->dB(pat::Muon::PV3D))/mu->edB(pat::Muon::PV3D)<<setw(10)<<fabs(mu->innerTrack()->dxy(firstGoodVertex.position()))<<setw(10)<<fabs(mu->innerTrack()->dz(firstGoodVertex.position()))<<setw(21)<<mu->segmentCompatibility()<<endl;
	}else{
	 cout<<setw(15)<<"no track for Mu IP"<<endl; 
	}
	amu = true;
      }*/
//BOOM->Scan("Muon_pt[0]:Muon_eta[0]:Muon_phi[0]:Muon_energy[0]:Muon_pdgId[0]:Muon_charge[0]:Muon_miniIsoRel[0]:Muon_miniIsoCh[0]:Muon_miniIsoPUsub[0]:Muon_ptrel[0]:Muon_jetcsv[0]:Muon_jetptratio[0]:Muon_IP3D_sig[0]:Muon_IP2D_val[0]:Muon_IP1D_val[0]:Muon_segmentCompatibility[0]")
