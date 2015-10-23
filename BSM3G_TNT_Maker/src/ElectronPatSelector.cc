#include "BSMFramework/BSM3G_TNT_Maker/interface/ElectronPatSelector.h"
ElectronPatSelector::ElectronPatSelector(std::string name, TTree* tree, bool debug, const pset& iConfig, edm::ConsumesCollector && ic): 
  baseTree(name,tree,debug),
  electronVetoIdMapToken_(ic.consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electronVetoIdMap"))),
  electronLooseIdMapToken_(ic.consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electronLooseIdMap"))),
  electronMediumIdMapToken_(ic.consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electronMediumIdMap"))),
  electronTightIdMapToken_(ic.consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electronTightIdMap"))),
  eleHEEPIdMapToken_(ic.consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleHEEPIdMap")))
{
  _vertexInputTag      = iConfig.getParameter<edm::InputTag>("vertices");
  _beamSpot            = iConfig.getParameter<edm::InputTag>("beamSpot");
  _patElectronToken    = iConfig.getParameter<edm::InputTag>("patElectrons");
  _patElectron_pt_min  = iConfig.getParameter<double>("patElectron_pt_min");
  _patElectron_eta_max = iConfig.getParameter<double>("patElectron_eta_max");
  _vtx_ndof_min        = iConfig.getParameter<int>("vtx_ndof_min");
  _vtx_rho_max         = iConfig.getParameter<int>("vtx_rho_max");
  _vtx_position_z_max  = iConfig.getParameter<double>("vtx_position_z_max");
  pfToken_             = ic.consumes<pat::PackedCandidateCollection>(edm::InputTag("packedPFCandidates"));
  jetToken_            = iConfig.getParameter<edm::InputTag>("jets");
  SetBranches();
}
ElectronPatSelector::~ElectronPatSelector(){
  delete tree_;
}
void ElectronPatSelector::Fill(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  Clear();
  /////
  //   Recall collections
  /////  
  edm::Handle<reco::VertexCollection> vtx_h;
  iEvent.getByLabel(_vertexInputTag, vtx_h);
  edm::Handle<edm::View<pat::Electron> > electron_pat;
  iEvent.getByLabel(_patElectronToken, electron_pat);
  reco::BeamSpot beamSpot;
  edm::Handle<reco::BeamSpot> beamSpotHandle;
  iEvent.getByLabel(_beamSpot, beamSpotHandle);
  edm::ESHandle<TransientTrackBuilder> ttrkbuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",ttrkbuilder);
  edm::Handle<edm::ValueMap<bool> > veto_id_decisions;
  edm::Handle<edm::ValueMap<bool> > loose_id_decisions;
  edm::Handle<edm::ValueMap<bool> > medium_id_decisions;
  edm::Handle<edm::ValueMap<bool> > tight_id_decisions;
  edm::Handle<edm::ValueMap<bool> > heep_id_decisions;
  iEvent.getByToken(electronVetoIdMapToken_,veto_id_decisions);
  iEvent.getByToken(electronLooseIdMapToken_,loose_id_decisions);
  iEvent.getByToken(electronMediumIdMapToken_,medium_id_decisions);
  iEvent.getByToken(electronTightIdMapToken_,tight_id_decisions);  
  iEvent.getByToken(eleHEEPIdMapToken_, heep_id_decisions);
  edm::Handle<double> rhopogHandle;
  iEvent.getByLabel("fixedGridRhoFastjetAll",rhopogHandle);
  double rhopog = *rhopogHandle;
  edm::Handle<double> rhotthHandle;
  iEvent.getByLabel("fixedGridRhoFastjetCentralNeutral",rhotthHandle);
  double rhotth = *rhotthHandle;
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
  /////
  //   Get electron information 
  /////
  for(edm::View<pat::Electron>::const_iterator el = electron_pat->begin(); el != electron_pat->end(); el++){
    //Acceptance
    if(el->pt() < _patElectron_pt_min)         continue;
    if(fabs(el->eta()) > _patElectron_eta_max) continue;  
    //Kinematics    
    patElectron_pt.push_back(el->pt());
    patElectron_eta.push_back(el->eta());
    patElectron_phi.push_back(el->phi());
    patElectron_energy.push_back(el->energy());
    patElectron_Et.push_back(el->caloEnergy()*sin(el->p4().theta()));
    double absEleSCeta = el->superCluster()->position().eta();
    patElectron_SCeta.push_back(absEleSCeta);
    bool inCrack  = (absEleSCeta>1.4442 && absEleSCeta<1.5660);
    patElectron_inCrack.push_back(inCrack);
    //Charge
    patElectron_charge.push_back(el->charge());
    //ID
    const Ptr<pat::Electron> elPtr(electron_pat, el - electron_pat->begin() );
    bool isPassVeto   = (*veto_id_decisions)  [ elPtr ];
    bool isPassLoose  = (*loose_id_decisions) [ elPtr ];
    bool isPassMedium = (*medium_id_decisions)[ elPtr ];
    bool isPassTight  = (*tight_id_decisions) [ elPtr ];
    bool isHEEPId     = (*heep_id_decisions)  [ elPtr ];
    passVetoId_.push_back  ( isPassVeto   );
    passLooseId_.push_back ( isPassLoose  );
    passMediumId_.push_back( isPassMedium );
    passTightId_.push_back ( isPassTight  );
    passHEEPId_.push_back  ( isHEEPId     );   
    patElectron_pdgId.push_back(el->pdgId());
    patElectron_isEcalDriven.push_back(el->ecalDriven());
    //Isolation
    //reco::GsfElectron::PflowIsolationVariables pfIso = el->pfIsolationVariables();
    double SumChHadPt       = el->pfIsolationVariables().sumChargedHadronPt;
    double SumNeuHadEt      = el->pfIsolationVariables().sumNeutralHadronEt;
    double SumPhotonEt      = el->pfIsolationVariables().sumPhotonEt; 
    double SumPU            = el->pfIsolationVariables().sumPUPt;
    patElectron_isoChargedHadrons.push_back( SumChHadPt );
    patElectron_isoNeutralHadrons.push_back( SumNeuHadEt );
    patElectron_isoPhotons.push_back( SumPhotonEt );
    patElectron_isoPU.push_back( SumPU );
    double SumNeutralCorrEt = std::max( 0.0, SumNeuHadEt+SumPhotonEt - 0.5*SumPU );
    double relIsoDeltaBeta = (SumChHadPt + SumNeutralCorrEt)/el->pt();
    patElectron_relIsoDeltaBeta.push_back(relIsoDeltaBeta);
    double eleEta = fabs(el->eta());
    double EffArea = 0;
    if (eleEta >= 0. && eleEta < 0.8)        EffArea = 0.1013;
    else if (eleEta >= 0.8 && eleEta < 1.3)  EffArea = 0.0988;
    else if (eleEta >= 1.3 && eleEta < 2.0)  EffArea = 0.0572;
    else if (eleEta >= 2.0 && eleEta < 2.2)  EffArea = 0.0842;
    else if (eleEta >= 2.2 && eleEta <= 2.5) EffArea = 0.1530;
    SumNeutralCorrEt = std::max( 0.0, SumNeuHadEt+SumPhotonEt - rhopog*EffArea );
    double relIsoRhoEA = (SumChHadPt + SumNeutralCorrEt)/el->pt();
    patElectron_relIsoRhoEA.push_back(relIsoRhoEA);
    patElectron_dr03EcalRecHitSumEt.push_back(el->dr03EcalRecHitSumEt());
    patElectron_dr03HcalDepth1TowerSumEt.push_back(el->dr03HcalDepth1TowerSumEt());
    patElectron_isolPtTracks.push_back(el->dr03TkSumPt());
    //Shape, Track related variables, other prop
    double dEtaIn = el->deltaEtaSuperClusterTrackAtVtx();
    double dPhiIn = el->deltaPhiSuperClusterTrackAtVtx();
    double full5x5_sigmaIetaIeta = el->full5x5_sigmaIetaIeta();
    double hOverE = el->hcalOverEcal();
    double ooEmooP = -999;
    if(el->ecalEnergy()==0)                   ooEmooP = 1e30;
    else if(!std::isfinite(el->ecalEnergy())) ooEmooP = 1e30;
    else                                      ooEmooP = fabs(1.0/el->ecalEnergy() - el->eSuperClusterOverP()/el->ecalEnergy() );
    patElectron_dEtaIn.push_back(dEtaIn);
    patElectron_dPhiIn.push_back(dPhiIn);
    patElectron_full5x5_sigmaIetaIeta.push_back(full5x5_sigmaIetaIeta);
    patElectron_full5x5_e2x5Max.push_back(el->full5x5_e2x5Max());
    patElectron_full5x5_e5x5.push_back(el->full5x5_e5x5());
    patElectron_full5x5_e1x5.push_back(el->full5x5_e1x5());
    patElectron_hOverE.push_back(hOverE);
    patElectron_ooEmooP.push_back(ooEmooP);
    passConversionVeto_.push_back(el->passConversionVeto());
    if(el->gsfTrack().isNonnull()){
      expectedMissingInnerHits.push_back(el->gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS));
      patElectron_gsfTrack_normChi2.push_back(el->gsfTrack()->normalizedChi2());
      patElectron_gsfTrack_ndof.push_back(el->gsfTrack()->ndof());
    }else{
      expectedMissingInnerHits.push_back(-999);
      patElectron_gsfTrack_normChi2.push_back(-999);
      patElectron_gsfTrack_ndof.push_back(-999);
    }
    //IP
    if(el->gsfTrack().isNonnull()){
      patElectron_gsfTrack_dz_pv.push_back(fabs(el->gsfTrack()->dz(firstGoodVertex.position())));
      patElectron_gsfTrack_dxy_pv.push_back(fabs(el->gsfTrack()->dxy(firstGoodVertex.position())));
      patElectron_d0.push_back((-1) * el->gsfTrack()->dxy(firstGoodVertex.position()));
      patElectron_dzError.push_back(el->gsfTrack()->dzError());
      patElectron_dxyError.push_back(el->gsfTrack()->d0Error());
      patElectron_gsfTrack_vtx.push_back(el->gsfTrack()->vx());
      patElectron_gsfTrack_vty.push_back(el->gsfTrack()->vy());
      patElectron_gsfTrack_vtz.push_back(el->gsfTrack()->vz());
      if(beamSpotHandle.isValid() && el->closestCtfTrackRef().isNonnull()){//AJ vars (both pv and bs are in this if condition, tought for pv is not mandatory)
        beamSpot = *beamSpotHandle;
        math::XYZPoint point(beamSpot.x0(),beamSpot.y0(), beamSpot.z0());
        patElectron_gsfTrack_dz_bs.push_back(el->gsfTrack()->dz(point));
        patElectron_gsfTrack_dxy_bs.push_back(el->gsfTrack()->dxy(point));
        GlobalPoint thebs(beamSpot.x0(),beamSpot.y0(),beamSpot.z0()); 
        GlobalPoint thepv(firstGoodVertex.position().x(),firstGoodVertex.position().y(),firstGoodVertex.position().z());
        TrackRef eletr = el->closestCtfTrackRef(); 
        TransientTrack elecTransTkPtr = ttrkbuilder->build(eletr);
        GlobalPoint patElectron_pca_pv = elecTransTkPtr.trajectoryStateClosestToPoint(thepv).position();
        GlobalPoint patElectron_pca_bs = elecTransTkPtr.trajectoryStateClosestToPoint(thebs).position();
        patElectron_gsfTrack_PCAx_pv.push_back(patElectron_pca_pv.x());
        patElectron_gsfTrack_PCAy_pv.push_back(patElectron_pca_pv.y());
        patElectron_gsfTrack_PCAz_pv.push_back(patElectron_pca_pv.z());
        patElectron_gsfTrack_PCAx_bs.push_back(patElectron_pca_bs.x());
        patElectron_gsfTrack_PCAy_bs.push_back(patElectron_pca_bs.y());
        patElectron_gsfTrack_PCAz_bs.push_back(patElectron_pca_bs.z());
        const float elecMass = 0.000510998928;
        float elecSigma      = elecMass*1e-6;
        float chi2 = 0.0;
        float ndf  = 0.0;
        KinematicParticleFactoryFromTransientTrack pFactory;
        RefCountedKinematicParticle elecParticle = pFactory.particle(elecTransTkPtr, elecMass, chi2, ndf, elecSigma);
        patElectron_gsfTrackFitErrorMatrix_00.push_back(elecParticle->stateAtPoint(patElectron_pca_bs).kinematicParametersError().matrix()(0,0));
        patElectron_gsfTrackFitErrorMatrix_01.push_back(elecParticle->stateAtPoint(patElectron_pca_bs).kinematicParametersError().matrix()(0,1));
        patElectron_gsfTrackFitErrorMatrix_02.push_back(elecParticle->stateAtPoint(patElectron_pca_bs).kinematicParametersError().matrix()(0,2));
        patElectron_gsfTrackFitErrorMatrix_11.push_back(elecParticle->stateAtPoint(patElectron_pca_bs).kinematicParametersError().matrix()(1,1));
        patElectron_gsfTrackFitErrorMatrix_12.push_back(elecParticle->stateAtPoint(patElectron_pca_bs).kinematicParametersError().matrix()(1,2));
        patElectron_gsfTrackFitErrorMatrix_22.push_back(elecParticle->stateAtPoint(patElectron_pca_bs).kinematicParametersError().matrix()(2,2));
      }else{
        patElectron_gsfTrack_dz_bs.push_back(-998);
        patElectron_gsfTrack_dxy_bs.push_back(-998);
        patElectron_gsfTrack_PCAx_pv.push_back(-998);
        patElectron_gsfTrack_PCAy_pv.push_back(-998);
        patElectron_gsfTrack_PCAz_pv.push_back(-998);
        patElectron_gsfTrack_PCAx_bs.push_back(-998);
        patElectron_gsfTrack_PCAy_bs.push_back(-998);
        patElectron_gsfTrack_PCAz_bs.push_back(-998);
        patElectron_gsfTrackFitErrorMatrix_00.push_back(-998);
        patElectron_gsfTrackFitErrorMatrix_01.push_back(-998);
        patElectron_gsfTrackFitErrorMatrix_02.push_back(-998);
        patElectron_gsfTrackFitErrorMatrix_11.push_back(-998);
        patElectron_gsfTrackFitErrorMatrix_12.push_back(-998);
        patElectron_gsfTrackFitErrorMatrix_22.push_back(-998);
      }
    }else{
      patElectron_gsfTrack_dz_pv.push_back(-999);
      patElectron_gsfTrack_dxy_pv.push_back(-999);
      patElectron_d0.push_back(-999);
      patElectron_dzError.push_back(-999);
      patElectron_dxyError.push_back(-999);
      patElectron_gsfTrack_vtx.push_back(-999);
      patElectron_gsfTrack_vty.push_back(-999);
      patElectron_gsfTrack_vtz.push_back(-999);
      patElectron_gsfTrack_dz_bs.push_back(-999);
      patElectron_gsfTrack_dxy_bs.push_back(-999);
      patElectron_gsfTrack_PCAx_pv.push_back(-999);
      patElectron_gsfTrack_PCAy_pv.push_back(-999);
      patElectron_gsfTrack_PCAz_pv.push_back(-999);
      patElectron_gsfTrack_PCAx_bs.push_back(-999);
      patElectron_gsfTrack_PCAy_bs.push_back(-999);
      patElectron_gsfTrack_PCAz_bs.push_back(-999);
      patElectron_gsfTrackFitErrorMatrix_00.push_back(-999);
      patElectron_gsfTrackFitErrorMatrix_01.push_back(-999);
      patElectron_gsfTrackFitErrorMatrix_02.push_back(-999);
      patElectron_gsfTrackFitErrorMatrix_11.push_back(-999);
      patElectron_gsfTrackFitErrorMatrix_12.push_back(-999);
      patElectron_gsfTrackFitErrorMatrix_22.push_back(-999);
    }
    /////
    //   TTH variables
    ///// 
    //cout<<setw(20)<<"Event num ELECTRON"<<setw(20)<<iEvent.id().event()<<endl;
    //cout<<setw(20)<<"Ele kin"<<setw(20)<<el->pt()<<setw(20)<<el->eta()<<setw(20)<<el->phi()<<setw(20)<<el->pdgId()<<setw(20)<<el->charge()<<endl;     
    double miniIso      = 999;
    double miniIsoCh    = 999;
    double miniIsoNeu   = 999;
    double miniIsoPUsub = 999;
    get_eleminiIso_info(*pcc,rhotth,*el,miniIso,miniIsoCh,miniIsoNeu,miniIsoPUsub);
    //cout<<setw(20)<<"Ele iso"<<setw(20)<<miniIso/el->pt()<<setw(20)<<miniIsoCh/el->pt()<<setw(20)<<miniIsoNeu/el->pt()<<endl;
    //cout<<setw(20)<<"Ele iso"<<setw(20)<<miniIso/el->pt()<<setw(20)<<miniIsoCh<<setw(20)<<miniIsoNeu<<endl;
    double elejet_mindr    = 9999;
    double elejet_pt       = -1;
    double eleptOVelejetpt = -1;
    double elejet_btagdisc = -1;
    double elejetx  = -9999;
    double elejety  = -9999;
    double elejetz  = -9999;
    double eleptrel = -9999;
    get_elejet_info(el,iEvent,iSetup,elejet_mindr,elejet_pt,eleptOVelejetpt,elejet_btagdisc,elejetx,elejety,elejetz,eleptrel);
    //cout<<setw(20)<<"Jet related var"<<setw(20)<<eleptrel<<setw(20)<<elejet_btagdisc<<setw(20)<<el->pt()/elejet_pt<<endl;
    //cout<<setw(20)<<"Ele IP"<<setw(20)<<fabs(el->dB(pat::Electron::PV3D))/el->edB(pat::Electron::PV3D)<<setw(20)<<fabs(el->gsfTrack()->dxy(firstGoodVertex.position()))<<setw(20)<<fabs(el->gsfTrack()->dz(firstGoodVertex.position()))<<endl;//setw(20)<<endl;//el->segmentCompatibility()<<endl;
    //if(el->hasUserFloat("ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Values")) cout<<setw(20)<<"Ele MVA"<<setw(20)<<el->userFloat("ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Values")<<endl;
    //else cout<<setw(20)<<"No Ele MVA"<<endl;
    patElectron_miniIsoRel.push_back(miniIso/el->pt());
    patElectron_miniIsoCh.push_back(miniIsoCh);
    patElectron_miniIsoNeu.push_back(miniIsoNeu);
    patElectron_miniIsoPUsub.push_back(miniIsoPUsub);
    patElectron_jetdr.push_back(elejet_mindr);
    patElectron_jetpt.push_back(elejet_pt);
    patElectron_jetptratio.push_back(eleptOVelejetpt);
    patElectron_jetcsv.push_back(elejet_btagdisc);
    patElectron_ptrel.push_back(eleptrel);
    patElectron_IP3Dsig.push_back(fabs(el->dB(pat::Electron::PV3D))/el->edB(pat::Electron::PV3D));
    patElectron_dxy.push_back(fabs(el->gsfTrack()->dxy(firstGoodVertex.position())));
    if(el->hasUserFloat("ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Values")) patElectron_eleMVASpring15NonTrig25ns.push_back(el->userFloat("ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Values"));   
    else                                                                        patElectron_eleMVASpring15NonTrig25ns.push_back(-999); 
  }
  //cout<<endl;
}
void ElectronPatSelector::SetBranches(){
  if(debug_) std::cout<<"setting branches: calling AddBranch of baseTree"<<std::endl;
  //Kinematics     
  AddBranch(&patElectron_pt           ,"patElectron_pt");
  AddBranch(&patElectron_eta          ,"patElectron_eta");
  AddBranch(&patElectron_phi          ,"patElectron_phi");
  AddBranch(&patElectron_energy       ,"patElectron_energy");
  AddBranch(&patElectron_Et           ,"patElectron_Et");
  AddBranch(&patElectron_SCeta        ,"patElectron_SCeta");
  AddBranch(&patElectron_inCrack      ,"patElectron_inCrack");
  //Charge
  AddBranch(&patElectron_charge       ,"patElectron_charge");
  //ID
  AddBranch(&passVetoId_              ,"patElectron_isPassVeto");          
  AddBranch(&passLooseId_             ,"patElectron_isPassLoose");
  AddBranch(&passMediumId_            ,"patElectron_isPassMedium");
  AddBranch(&passTightId_             ,"patElectron_isPassTight");
  AddBranch(&passHEEPId_              ,"patElectron_isPassHEEPId");
  AddBranch(&patElectron_pdgId        ,"patElectron_pdgId");
  AddBranch(&patElectron_isEcalDriven ,"patElectron_isEcalDriven");
  //Isolation
  AddBranch(&patElectron_isoChargedHadrons        ,"patElectron_isoChargedHadrons");
  AddBranch(&patElectron_isoNeutralHadrons        ,"patElectron_isoNeutralHadrons");
  AddBranch(&patElectron_isoPhotons               ,"patElectron_isoPhotons");
  AddBranch(&patElectron_isoPU                    ,"patElectron_isoPU");
  AddBranch(&patElectron_relIsoDeltaBeta          ,"patElectron_relIsoDeltaBeta");
  AddBranch(&patElectron_relIsoRhoEA              ,"patElectron_relIsoRhoEA");
  AddBranch(&patElectron_dr03EcalRecHitSumEt      ,"patElectron_dr03EcalRecHitSumEt");
  AddBranch(&patElectron_dr03HcalDepth1TowerSumEt ,"patElectron_dr03HcalDepth1TowerSumEt");
  AddBranch(&patElectron_isolPtTracks             ,"patElectron_isolPtTracks");
  //Shape, Track related variables, other prop
  AddBranch(&patElectron_dEtaIn                ,"patElectron_dEtaIn");
  AddBranch(&patElectron_dPhiIn                ,"patElectron_dPhiIn");
  AddBranch(&patElectron_full5x5_sigmaIetaIeta ,"patElectron_full5x5_sigmaIetaIeta");
  AddBranch(&patElectron_full5x5_e2x5Max       ,"patElectron_full5x5_e2x5Max");
  AddBranch(&patElectron_full5x5_e5x5          ,"patElectron_full5x5_e5x5");
  AddBranch(&patElectron_full5x5_e1x5          ,"patElectron_full5x5_e1x5");
  AddBranch(&patElectron_hOverE                ,"patElectron_hOverE");
  AddBranch(&patElectron_ooEmooP               ,"patElectron_ooEmooP");
  AddBranch(&passConversionVeto_               ,"patElectron_passConversionVeto"); 
  AddBranch(&expectedMissingInnerHits          ,"patElectron_expectedMissingInnerHits");
  AddBranch(&patElectron_gsfTrack_ndof         ,"patElectron_gsfTrack_ndof");
  AddBranch(&patElectron_gsfTrack_normChi2     ,"patElectron_gsfTrack_normChi2");
  //Vertex compatibility
  AddBranch(&patElectron_d0             ,"patElectron_d0");
  AddBranch(&patElectron_gsfTrack_dz_pv ,"patElectron_gsfTrack_dz_pv");
  //TTH
  AddBranch(&patElectron_miniIsoRel                ,"patElectron_miniIsoRel");
  AddBranch(&patElectron_miniIsoCh                 ,"patElectron_miniIsoCh");
  AddBranch(&patElectron_miniIsoNeu                ,"patElectron_miniIsoNeu");
  AddBranch(&patElectron_miniIsoPUsub              ,"patElectron_miniIsoPUsub");
  AddBranch(&patElectron_jetdr                     ,"patElectron_jetdr");
  AddBranch(&patElectron_jetpt                     ,"patElectron_jetpt");
  AddBranch(&patElectron_jetptratio                ,"patElectron_jetptratio");
  AddBranch(&patElectron_jetcsv                    ,"patElectron_jetcsv");
  AddBranch(&patElectron_ptrel                     ,"patElectron_ptrel");
  AddBranch(&patElectron_IP3Dsig                   ,"patElectron_IP3Dsig");
  AddBranch(&patElectron_dxy                       ,"patElectron_dxy");
  AddBranch(&patElectron_eleMVASpring15NonTrig25ns ,"patElectron_eleMVASpring15NonTrig25ns");
  if(debug_) std::cout<<"set branches"<<std::endl;
}
void ElectronPatSelector::Clear(){
  //Kinematics     
  patElectron_pt.clear();
  patElectron_eta.clear();
  patElectron_phi.clear();
  patElectron_energy.clear();
  patElectron_Et.clear();
  patElectron_SCeta.clear();
  patElectron_inCrack.clear();
  //Charge
  patElectron_charge.clear(); 
  //ID
  passVetoId_.clear();
  passLooseId_.clear();
  passMediumId_.clear();
  passTightId_.clear();  
  passHEEPId_.clear();
  patElectron_pdgId.clear();
  patElectron_isEcalDriven.clear();
  //Isolation
  patElectron_isoChargedHadrons.clear();
  patElectron_isoNeutralHadrons.clear();
  patElectron_isoPhotons.clear();
  patElectron_isoPU.clear();
  patElectron_relIsoDeltaBeta.clear();
  patElectron_relIsoRhoEA.clear();
  patElectron_dr03EcalRecHitSumEt.clear();
  patElectron_dr03HcalDepth1TowerSumEt.clear();
  patElectron_isolPtTracks.clear();
  //Shape
  patElectron_dEtaIn.clear();
  patElectron_dPhiIn.clear();
  patElectron_full5x5_sigmaIetaIeta.clear();
  patElectron_full5x5_e2x5Max.clear();
  patElectron_full5x5_e5x5.clear();
  patElectron_full5x5_e1x5.clear();
  patElectron_hOverE.clear();
  patElectron_ooEmooP.clear();
  passConversionVeto_.clear();
  expectedMissingInnerHits.clear();
  patElectron_gsfTrack_ndof.clear();
  patElectron_gsfTrack_normChi2.clear();
  //Vertex compatibility
  patElectron_d0.clear();
  patElectron_gsfTrack_dz_pv.clear();
  //TTH
  patElectron_miniIsoRel.clear();
  patElectron_miniIsoCh.clear();
  patElectron_miniIsoNeu.clear();
  patElectron_miniIsoPUsub.clear();
  patElectron_jetdr.clear();
  patElectron_jetpt.clear();
  patElectron_jetptratio.clear();
  patElectron_jetcsv.clear();
  patElectron_ptrel.clear();
  patElectron_IP3Dsig.clear();
  patElectron_dxy.clear();
  patElectron_eleMVASpring15NonTrig25ns.clear();
}
bool ElectronPatSelector::isGoodVertex(const reco::Vertex& vtx){
  if(vtx.isFake())                                   return false;
  if(vtx.ndof()<_vtx_ndof_min)                       return false;
  if(vtx.position().Rho()>_vtx_rho_max)              return false;
  if(fabs(vtx.position().Z()) > _vtx_position_z_max) return false;
  return true;
}
void ElectronPatSelector::get_eleminiIso_info(const pat::PackedCandidateCollection& pcc,double rhotth, const pat::Electron& cand, double& miniIso, double& miniIsoCh, double& miniIsoNeu, double& miniIsoPUsub){
  double miniIsoConeSize = 10.0/min(max(cand.pt(), 50.),200.);
  vector<const pat::PackedCandidate *> pfc_all; pfc_all.clear();
  vector<const pat::PackedCandidate *> pfc_ch;  pfc_ch.clear();
  vector<const pat::PackedCandidate *> pfc_neu; pfc_neu.clear();
  vector<const pat::PackedCandidate *> pfc_pu;  pfc_pu.clear();
  get_chneupu_pcc(pcc,pfc_all,pfc_ch,pfc_neu,pfc_pu);
  double innerR_ch;
  double innerR_neu;
  if(cand.isEB()){
    innerR_ch  = 0.0;
    innerR_neu = 0.0;
  }else{
    innerR_ch  = 0.015;
    innerR_neu = 0.08;
  }
  miniIsoCh  = get_isosumraw(pfc_ch,  cand, miniIsoConeSize, innerR_ch,  0.5, 0);
  miniIsoNeu = get_isosumraw(pfc_neu, cand, miniIsoConeSize, innerR_neu, 0.5, 0);
  double effarea    = get_effarea(cand.eta());
  double correction = rhotth*effarea*pow((miniIsoConeSize/0.3),2);
  miniIsoPUsub = std::max(0.0, miniIsoNeu-correction);
  miniIso = miniIsoCh+miniIsoPUsub;
}
void ElectronPatSelector::get_chneupu_pcc(const pat::PackedCandidateCollection& pcc,vector<const pat::PackedCandidate *>& pfc_all,vector<const pat::PackedCandidate *>& pfc_ch,vector<const pat::PackedCandidate *>& pfc_neu,vector<const pat::PackedCandidate *>&pfc_pu){
  for(const pat::PackedCandidate &p : pcc){
    pfc_all.push_back(&p);
    if(p.charge()==0){
     pfc_neu.push_back(&p);
    }else{
      if(1){//(abs(p.pdgId())==211)){// || ((abs(p.pdgId()) == 11 ) || (abs(p.pdgId()) == 13 )) )
        if(p.fromPV()>1 && fabs(p.dz())<9999){
          pfc_ch.push_back(&p);
        }else{
          pfc_pu.push_back(&p);
        }
      }
    }
  }
}
double ElectronPatSelector::get_isosumraw(const std::vector<const pat::PackedCandidate *> & pcc, const pat::Electron& cand, double miniIsoConeSize, double innerR, double ptTh, int pdgId){
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
    //if(ptTh==0.5) cout<<"pt is "<<setw(20)<<(*pc)->pt()<<setw(20)<<deltaR(**pc, cand)<<endl;
    //itself veto
    if(std::find(vetos.begin(), vetos.end(),*pc)!=vetos.end()) continue;
    //add to sum
    isosum += (*pc)->pt();
  }
  return isosum;
}
double ElectronPatSelector::get_effarea(double eta){
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
void ElectronPatSelector::get_elejet_info(edm::View<pat::Electron>::const_iterator& ele, const edm::Event& iEvent, const edm::EventSetup& iSetup, double& elejet_mindr, double& elejet_pt, double& eleptOVelejetpt, double& elejet_btagdisc, double& jx, double& jy, double& jz, double& eleptrel){
 //Look for jet associated to ele
 edm::Handle<pat::JetCollection> jets;
 iEvent.getByLabel(jetToken_, jets);
 pat::Jet elejet;
 //const JetCorrector* corrector = JetCorrector::getJetCorrector( "ak4PFchsL1L2L3", iSetup );
 for(const pat::Jet &j : *jets){
  pat::Jet jet = j;//j.correctedJet(0);
  //double scale = corrector->correction(jet, iEvent, iSetup);
  //jet.scaleEnergy(scale);
  double dr = deltaR(ele->p4(),jet.p4());
  if(dr<elejet_mindr){
   elejet_mindr = dr;
   elejet       = jet;
  }
 }
 //Get info
 double L2L3_SF = elejet.p4().E()/elejet.correctedJet(1).p4().E();
 elejet.setP4(((elejet.correctedJet(1).p4()-ele->p4())*L2L3_SF)+ele->p4());
 elejet_pt       = elejet.pt();
 eleptOVelejetpt = min(ele->pt()/elejet.pt(), 1.5);
 elejet_btagdisc = max(double(elejet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags")), 0.0);
 jx = elejet.px();
 jy = elejet.py();
 jz = elejet.pz();
 TLorentzVector ele_lv = TLorentzVector(ele->px(),ele->py(),ele->pz(),ele->p4().E());
 TLorentzVector jet_lv = TLorentzVector(elejet.px(),elejet.py(),elejet.pz(),elejet.p4().E());
 eleptrel = ele_lv.Perp((jet_lv-ele_lv).Vect());
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

