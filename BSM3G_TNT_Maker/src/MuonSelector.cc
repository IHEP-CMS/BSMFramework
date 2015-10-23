#include "BSMFramework/BSM3G_TNT_Maker/interface/MuonSelector.h"
MuonSelector::MuonSelector(std::string name, TTree* tree, bool debug, const pset& iConfig, edm::ConsumesCollector && ic):
  baseTree(name,tree,debug)
{
  _muonToken          = iConfig.getParameter<edm::InputTag>("muons");
  _vertexInputTag     = iConfig.getParameter<edm::InputTag>("vertices");
  _beamSpot           = iConfig.getParameter<edm::InputTag>("beamSpot");
  _Muon_pt_min        = iConfig.getParameter<double>("Muon_pt_min");
  _Muon_eta_max       = iConfig.getParameter<double>("Muon_eta_max");
  _vtx_ndof_min       = iConfig.getParameter<int>("vtx_ndof_min");
  _vtx_rho_max        = iConfig.getParameter<int>("vtx_rho_max");
  _vtx_position_z_max = iConfig.getParameter<double>("vtx_position_z_max");
  _super_TNT          = iConfig.getParameter<bool>("super_TNT");
  pfToken_            = ic.consumes<pat::PackedCandidateCollection>(edm::InputTag("packedPFCandidates"));
  jetToken_           = iConfig.getParameter<edm::InputTag>("jets");
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
  edm::Handle<edm::View<pat::Muon> > muon_h;
  iEvent.getByLabel(_muonToken, muon_h);
  edm::Handle<reco::VertexCollection> vtx_h;
  iEvent.getByLabel(_vertexInputTag, vtx_h);
  reco::BeamSpot beamSpot;
  edm::Handle<reco::BeamSpot> beamSpotHandle;
  iEvent.getByLabel(_beamSpot, beamSpotHandle);
  edm::ESHandle<TransientTrackBuilder> ttrkbuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",ttrkbuilder);
  edm::Handle<double> rhoHandle;
  iEvent.getByLabel("fixedGridRhoFastjetCentralNeutral",rhoHandle);
  double rho = *rhoHandle;
  edm::Handle<pat::PackedCandidateCollection> pcc;
  iEvent.getByToken(pfToken_, pcc);
  /////
  //   Require a good vertex 
  ///// 
  //reco::VertexCollection::const_iterator firstGoodVertex = vtx_h->end();
  //for(reco::VertexCollection::const_iterator it = vtx_h->begin(); it != firstGoodVertex; it++){
  //  isGoodVertex(*it);
  //  firstGoodVertex = it;
  //  break;
  //}
  //if(firstGoodVertex == vtx_h->end()) return;
  if(vtx_h->empty()) return; // skip the event if no PV found
  const reco::Vertex &firstGoodVertex = vtx_h->front();  
  bool isgoodvtx = isGoodVertex(firstGoodVertex);
  if(!isgoodvtx) return;
  ////
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
    Muon_ecalIso.push_back(mu->ecalIso());
    Muon_hcalIso.push_back(mu->hcalIso());
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
    /////
    //   TTH variables
    ///// 
    //cout<<setw(20)<<"Event num MUON"<<setw(20)<<iEvent.id().event()<<endl;
    //cout<<setw(20)<<"Mu kin"<<setw(20)<<mu->pt()<<setw(20)<<mu->eta()<<setw(20)<<mu->phi()<<setw(20)<<mu->energy()<<setw(20)<<mu->pdgId()<<setw(20)<<mu->charge()<<endl;  
    double miniIso      = 999;
    double miniIsoCh    = 999;
    double miniIsoNeu   = 999;
    double miniIsoPUsub = 999;
    get_muminiIso_info(*pcc,rho,*mu,miniIso,miniIsoCh,miniIsoNeu,miniIsoPUsub);
    //cout<<setw(20)<<"Mu iso"<<setw(20)<<miniIso/mu->pt()<<setw(20)<<miniIsoCh/mu->pt()<<setw(20)<<miniIsoNeu/mu->pt()<<endl;
    //cout<<setw(20)<<"Mu iso"<<setw(20)<<miniIso/mu->pt()<<setw(20)<<miniIsoCh<<setw(20)<<miniIsoNeu<<endl;
    double mujet_mindr    = 999;
    double mujet_pt       = -1;
    double muptOVmujetpt  = -1;
    double mujet_btagdisc = -1;
    double mujetx  = -999;
    double mujety  = -999;
    double mujetz  = -999;
    double muptrel = -999;
    get_mujet_info(*mu,iEvent,iSetup,mujet_mindr,mujet_pt,muptOVmujetpt,mujet_btagdisc,mujetx,mujety,mujetz,muptrel);
    //cout<<setw(20)<<"Jet related var"<<setw(20)<<muptrel<<setw(20)<<mujet_btagdisc<<setw(20)<<mu->pt()/mujet_pt<<endl;
    if(mu->innerTrack().isNonnull()){
     //cout<<setw(20)<<"Mu IP"<<setw(20)<<fabs(mu->dB(pat::Muon::PV3D))/mu->edB(pat::Muon::PV3D)<<setw(20)<<fabs(mu->innerTrack()->dxy(firstGoodVertex.position()))<<setw(20)<<fabs(mu->innerTrack()->dz(firstGoodVertex.position()))<<setw(20)<<mu->segmentCompatibility()<<endl;
    }else{
     //cout<<setw(20)<<"no track for Mu IP"<<endl; 
    }
    Muon_miniIsoRel.push_back(miniIso/mu->pt());
    Muon_miniIsoCh.push_back(miniIsoCh);
    Muon_miniIsoNeu.push_back(miniIsoNeu);
    Muon_miniIsoPUsub.push_back(miniIsoPUsub);
    Muon_jetdr.push_back(mujet_mindr);
    Muon_jetpt.push_back(mujet_pt);
    Muon_jetptratio.push_back(muptOVmujetpt);
    Muon_jetcsv.push_back(mujet_btagdisc);
    Muon_ptrel.push_back(muptrel);
    Muon_IP3Dsig_it.push_back(fabs(mu->dB(pat::Muon::PV3D))/mu->edB(pat::Muon::PV3D));
    //break;
  }
}
void MuonSelector::SetBranches(){
  if(debug_)    std::cout<<"setting branches: calling AddBranch of baseTree"<<std::endl;
  //Kinematics
  AddBranch(&Muon_pt                ,"Muon_pt");
  AddBranch(&Muon_eta               ,"Muon_eta");
  AddBranch(&Muon_phi               ,"Muon_phi");
  AddBranch(&Muon_energy            ,"Muon_energy");
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
  AddBranch(&Muon_ecalIso             ,"Muon_ecalIso");
  AddBranch(&Muon_hcalIso             ,"Muon_hcalIso");
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
  //TTH
  AddBranch(&Muon_miniIsoRel        ,"Muon_miniIsoRel");
  AddBranch(&Muon_miniIsoCh         ,"Muon_miniIsoCh");
  AddBranch(&Muon_miniIsoNeu        ,"Muon_miniIsoNeu");
  AddBranch(&Muon_miniIsoPUsub      ,"Muon_miniIsoPUsub");
  AddBranch(&Muon_jetdr             ,"Muon_jetdr");
  AddBranch(&Muon_jetpt             ,"Muon_jetpt");
  AddBranch(&Muon_jetptratio        ,"Muon_jetptratio");
  AddBranch(&Muon_jetcsv            ,"Muon_jetcsv");
  AddBranch(&Muon_ptrel             ,"Muon_ptrel");
  AddBranch(&Muon_IP3Dsig_it        ,"Muon_IP3Dsig_it");

  if(debug_) std::cout<<"set branches"<<std::endl;
}
void MuonSelector::Clear(){
  //Kinematics  
  Muon_pt.clear();
  Muon_eta.clear();
  Muon_phi.clear();
  Muon_energy.clear();
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
  Muon_ecalIso.clear();
  Muon_hcalIso.clear(); 
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
  //TTH
  Muon_miniIsoRel.clear();
  Muon_miniIsoCh.clear();
  Muon_miniIsoNeu.clear();
  Muon_miniIsoPUsub.clear();
  Muon_jetdr.clear();
  Muon_jetpt.clear();
  Muon_jetptratio.clear();
  Muon_jetcsv.clear();
  Muon_ptrel.clear();
  Muon_IP3Dsig_it.clear();
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
  //miniIsoCh  = get_isosumraw(pfc_ch,  cand, miniIsoConeSize, 0.0001, 0.0, 0);
  miniIsoCh  = get_isosumraw(pfc_ch,  cand, miniIsoConeSize, 0.01, 0.5, 0);
  miniIsoNeu = get_isosumraw(pfc_neu, cand, miniIsoConeSize, 0.01, 0.5, 0);
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
}
double MuonSelector::get_isosumraw(const std::vector<const pat::PackedCandidate *> & pcc, const pat::Muon& cand, double miniIsoConeSize, double innerR, double ptTh, int pdgId){
  //Look for cand sources
  std::vector<const reco::Candidate *> vetos; vetos.clear();
  for(uint i=0, n=cand.numberOfSourceCandidatePtrs(); i<n; ++i){
    const reco::CandidatePtr &cp = cand.sourceCandidatePtr(i);
    if(cp.isNonnull() && cp.isAvailable()){
      vetos.push_back(&*cp);
    }
  }     
  //Get the isolation
  double isosum = 0;
  for(std::vector<const pat::PackedCandidate *>::const_iterator pc = pcc.begin(); pc<pcc.end(); ++pc){
    //pdgId veto
    if(pdgId>0 && abs((*pc)->pdgId())!=pdgId) continue;
    //pT requirement 
    if(ptTh>0 && (*pc)->pt()<ptTh) continue;
    //cone region
    double dr = reco::deltaR(**pc, cand);
    if(dr<innerR || dr>=miniIsoConeSize) continue;
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
//double MuonSelector::get_iso_rho(const pat::Muon& mu, double& rho){
// double miniIsoCh  = mu.pfIsolationR03().sumChargedHadronPt;
// double miniIsoNeu = mu.pfIsolationR03().sumNeutralHadronEt + mu.pfIsolationR03().sumPhotonEt;
// double effarea  = get_effarea(mu.eta());
// double correction = rho*effarea;
	// double pfIsoPUsub = std::max( 0.0, miniIsoNeu - correction);
// double iso = (miniIsoCh + pfIsoPUsub)/mu.pt();
// return iso;
//}
void MuonSelector::get_mujet_info(const pat::Muon& mu, const edm::Event& iEvent, const edm::EventSetup& iSetup, double& mujet_mindr, double& mujet_pt, double& muptOVmujetpt, double& mujet_btagdisc, double& jx, double& jy, double& jz, double& muptrel){
 //Look for jet associated to mu
 edm::Handle<pat::JetCollection> jets;
 iEvent.getByLabel(jetToken_, jets);
 pat::Jet mujet;
 //const JetCorrector* corrector = JetCorrector::getJetCorrector( "ak4PFCHSL1L2L3Residual", iSetup );
 //const JetCorrector* corrector = JetCorrector::getJetCorrector( "ak4PFchsL1L2L3", iSetup );
 for(const pat::Jet &j : *jets){
  pat::Jet jet = j;//j.correctedJet(0);
  //double scale = corrector->correction(jet, iEvent, iSetup);
  //jet.scaleEnergy(scale);
  double dr = deltaR(mu.p4(),jet.p4());
  if(dr<mujet_mindr){
   mujet_mindr = dr;
   mujet       = jet;
  }
 }
 //Some info
 //cout<<"Jet info "<<mujet<<endl;
 //cout<<"Corrected (L0)"<<setw(20)<<mujet.correctedJet(0).p4().E()<<endl; 
 //cout<<"Corrected (L1)"<<setw(20)<<mujet.correctedJet(1).p4().E()<<endl;
 //cout<<"Corrected (L2)"<<setw(20)<<mujet.correctedJet(2).p4().E()<<endl;
 //cout<<"Corrected (L3)"<<setw(20)<<mujet.correctedJet(3).p4().E()<<endl;
 //cout<<"Corrected Fina"<<setw(20)<<mujet.p4().E()<<endl; 
 //Get info
 double L2L3_corr = mujet.p4().E()/mujet.correctedJet(1).p4().E(); 
 //cout<<"L2L3_corr"<<setw(20)<<L2L3_corr<<endl;
 mujet.setP4(((mujet.correctedJet(1).p4()-mu.p4())*L2L3_corr)+mu.p4());
 mujet_pt       = mujet.pt();
 muptOVmujetpt  = min(mu.pt()/mujet.pt(), 1.5);
 mujet_btagdisc = max(double(mujet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags")), 0.0);
 jx = mujet.px();
 jy = mujet.py();
 jz = mujet.pz();
 TLorentzVector mu_lv    = TLorentzVector(mu.px(),mu.py(),mu.pz(),mu.p4().E());
 TLorentzVector mujet_lv = TLorentzVector(mujet.px(),mujet.py(),mujet.pz(),mujet.p4().E());
 muptrel = mu_lv.Perp((mujet_lv-mu_lv).Vect());
}
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

