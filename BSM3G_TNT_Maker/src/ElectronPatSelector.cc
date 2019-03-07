#include "BSMFramework/BSM3G_TNT_Maker/interface/ElectronPatSelector.h"
ElectronPatSelector::ElectronPatSelector(std::string name, TTree* tree, bool debug, const pset& iConfig, edm::ConsumesCollector && ic): 
  baseTree(name,tree,debug),
  triggerBits_(ic.consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"))),
  triggerObjects_(ic.consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("objects"))),
  ebRecHitsToken_(ic.consumes<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>>>(iConfig.getParameter<edm::InputTag>("ebRecHits")))
{
  vtx_h_               = ic.consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"));
  beamSpot_            = ic.consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpot"));
  electron_pat_        = ic.consumes<edm::View<pat::Electron> >(iConfig.getParameter<edm::InputTag>("patElectrons"));
  pfToken_             = ic.consumes<pat::PackedCandidateCollection>(edm::InputTag("packedPFCandidates"));
  rhopogHandle_        = ic.consumes<double>(edm::InputTag("fixedGridRhoFastjetAll"));
  _patElectron_pt_min  = iConfig.getParameter<double>("patElectron_pt_min");
  _patElectron_eta_max = iConfig.getParameter<double>("patElectron_eta_max");
  _vtx_ndof_min        = iConfig.getParameter<int>("vtx_ndof_min");
  _vtx_rho_max         = iConfig.getParameter<int>("vtx_rho_max");
  _vtx_position_z_max  = iConfig.getParameter<double>("vtx_position_z_max");
  _is_data             = iConfig.getParameter<bool>("is_data");
  _reduced             = iConfig.getParameter<bool>("reduced");
  SetBranches();
}
ElectronPatSelector::~ElectronPatSelector(){
  delete tree_;
}
KalmanVertexFitter vertexfitterele(true);
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
void ElectronPatSelector::Fill(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  Clear();
  /////
  //   Recall collections
  /////  
  edm::Handle<reco::VertexCollection> vtx_h;
  iEvent.getByToken(vtx_h_, vtx_h);
  reco::BeamSpot beamSpot;
  edm::Handle<reco::BeamSpot> beamSpotHandle;
  iEvent.getByToken(beamSpot_, beamSpotHandle);
  edm::Handle<edm::View<pat::Electron> > electron_pat;
  iEvent.getByToken(electron_pat_, electron_pat);
  edm::Handle<pat::PackedCandidateCollection> pcc;
  iEvent.getByToken(pfToken_, pcc);
  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
  iEvent.getByToken(triggerObjects_, triggerObjects);
  edm::Handle<edm::TriggerResults> triggerBits;
  iEvent.getByToken(triggerBits_, triggerBits);
  edm::Handle<double> rhopogHandle;
  iEvent.getByToken(rhopogHandle_,rhopogHandle);
  double rhopog = *rhopogHandle;
  edm::ESHandle<TransientTrackBuilder> ttrkbuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",ttrkbuilder);
  iEvent.getByToken(ebRecHitsToken_, _ebRecHits);
  /////
  //   Require a good vertex 
  /////
  //if(vtx_h->empty()) return; // skip the event if no PV found
  const reco::Vertex &firstGoodVertex = vtx_h->front();  
  //bool isgoodvtx = isGoodVertex(firstGoodVertex);
  //if(!isgoodvtx) return;
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
    patElectron_px.push_back(el->px());
    patElectron_py.push_back(el->py());
    patElectron_pz.push_back(el->pz());
    patElectron_p.push_back(el->p());
    patElectron_Et.push_back(el->caloEnergy()*sin(el->p4().theta()));
    double EleSCeta = el->superCluster()->position().eta();
    patElectron_SCeta.push_back(EleSCeta);
    bool inCrack  = 1.4442<fabs(EleSCeta) && fabs(EleSCeta)<1.5660;
    patElectron_inCrack.push_back(inCrack);
    //Energy smearer
    if(el->energy()!=0) patElectron_energySF.push_back(el->userFloat("ecalTrkEnergyPostCorr")/el->energy());
    else patElectron_energySF.push_back(1.0);
    //Charge
    patElectron_charge.push_back(el->charge());
    //ID
    const Ptr<pat::Electron> elPtr(electron_pat, el - electron_pat->begin() );
    bool isPassVeto    = el->electronID("cutBasedElectronID-Fall17-94X-V1-veto");
    bool isPassLoose   = el->electronID("cutBasedElectronID-Fall17-94X-V1-loose");
    bool isPassMedium  = el->electronID("cutBasedElectronID-Fall17-94X-V1-medium");
    bool isPassTight   = el->electronID("cutBasedElectronID-Fall17-94X-V1-tight");
    bool isPassMvatrig = 0; //FIX THIS IF YOU NEED// (*mvatrig_id_decisions) [ elPtr ];
    bool isHEEPId      = 0; //FIX THIS IF YOU NEED// (*heep_id_decisions)  [ elPtr ];
    passVetoId_.push_back  ( isPassVeto   );
    passLooseId_.push_back ( isPassLoose  );
    passMediumId_.push_back( isPassMedium );
    passTightId_.push_back ( isPassTight  );
    passMvatrigId_.push_back( isPassMvatrig );
    passHEEPId_.push_back  ( isHEEPId     );
    //Isolation
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
    double EffArea = get_effarea(el->superCluster()->position().eta());
    SumNeutralCorrEt = std::max( 0.0, SumNeuHadEt+SumPhotonEt - rhopog*EffArea );
    double relIsoRhoEA = (SumChHadPt + SumNeutralCorrEt)/el->pt();
    patElectron_relIsoRhoEA.push_back(relIsoRhoEA);
    patElectron_dr03EcalRecHitSumEt.push_back(el->dr03EcalRecHitSumEt());
    patElectron_dr03HcalDepth1TowerSumEt.push_back(el->dr03HcalDepth1TowerSumEt());
    patElectron_isolPtTracks.push_back(el->dr03TkSumPt());
    patElectron_ecalPFClusterIso.push_back(el->ecalPFClusterIso());
    patElectron_hcalPFClusterIso.push_back(el->hcalPFClusterIso());
    //Shape, Track related variables, other prop
    double dEtaIn = el->deltaEtaSuperClusterTrackAtVtx();
    double dPhiIn = el->deltaPhiSuperClusterTrackAtVtx();
    double full5x5_sigmaIetaIeta = el->full5x5_sigmaIetaIeta();
    double hOverE = el->hcalOverEcal();
    double ooEmooP = -999;
    if(el->ecalEnergy()==0)                   ooEmooP = 1e30;
    else if(!std::isfinite(el->ecalEnergy())) ooEmooP = 1e30;
    else                                      ooEmooP = 1.0/el->ecalEnergy() - el->eSuperClusterOverP()/el->ecalEnergy();
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
      expectedMissingInnerHits.push_back(el->gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS));
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
    }else{
      patElectron_gsfTrack_dz_pv.push_back(-999);
      patElectron_gsfTrack_dxy_pv.push_back(-999);
      patElectron_d0.push_back(-999);
      patElectron_dzError.push_back(-999);
      patElectron_dxyError.push_back(-999);
      patElectron_gsfTrack_vtx.push_back(-999);
      patElectron_gsfTrack_vty.push_back(-999);
      patElectron_gsfTrack_vtz.push_back(-999);
    }
    /////
    //   MC info
    /////
    if(!_is_data){
      const reco::GenParticle * genpart = el->genParticle();
      if(genpart){
        patElectron_gen_pt.push_back(genpart->pt());
        patElectron_gen_eta.push_back(genpart->eta());
        patElectron_gen_phi.push_back(genpart->phi());
        patElectron_gen_en.push_back(genpart->energy());
        patElectron_gen_pdgId.push_back(genpart->pdgId());
        patElectron_gen_isPromptFinalState.push_back(genpart->isPromptFinalState());
        patElectron_gen_isDirectPromptTauDecayProductFinalState.push_back(genpart->isDirectPromptTauDecayProductFinalState());
      }else{
        patElectron_gen_pt.push_back(-999);
        patElectron_gen_eta.push_back(-999);
        patElectron_gen_phi.push_back(-999);
        patElectron_gen_en.push_back(-999);
        patElectron_gen_pdgId.push_back(-999);
        patElectron_gen_isPromptFinalState.push_back(-999);
        patElectron_gen_isDirectPromptTauDecayProductFinalState.push_back(-999);
      }
    }
    //TRIGGER MATCHING - for now only HLT_Ele115_CaloIdVT_GsfTrkIdT_v* is implemented
    int isMatchedToTrigger = 0;
    //if(!_is_data) isMatchedToTrigger = MatchingToTrigger(iEvent, triggerObjects, triggerBits, el->eta(), el->phi());
    patElectron_isMatchedToTrigger.push_back(isMatchedToTrigger);
  }
}
void ElectronPatSelector::SetBranches(){
  if(debug_) std::cout<<"setting branches: calling AddBranch of baseTree"<<std::endl;
  //Kinematics     
  AddBranch(&patElectron_pt           ,"patElectron_pt");
  AddBranch(&patElectron_eta          ,"patElectron_eta");
  AddBranch(&patElectron_phi          ,"patElectron_phi");
  AddBranch(&patElectron_energy       ,"patElectron_energy");
  AddBranch(&patElectron_energySF     ,"patElectron_energySF");
  if(!_reduced){
    AddBranch(&patElectron_energyCorr   ,"patElectron_energyCorr");
    AddBranch(&patElectron_px           ,"patElectron_px");
    AddBranch(&patElectron_py           ,"patElectron_py");
    AddBranch(&patElectron_pz           ,"patElectron_pz");
    AddBranch(&patElectron_p            ,"patElectron_p");
    AddBranch(&patElectron_Et           ,"patElectron_Et");
  }
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
  AddBranch(&passMvatrigId_           ,"patElectron_isPassMvatrig");
  //Isolation
  if(!_reduced){
    AddBranch(&patElectron_isoChargedHadrons        ,"patElectron_isoChargedHadrons");
    AddBranch(&patElectron_isoNeutralHadrons        ,"patElectron_isoNeutralHadrons");
    AddBranch(&patElectron_isoPhotons               ,"patElectron_isoPhotons");
    AddBranch(&patElectron_isoPU                    ,"patElectron_isoPU");
    AddBranch(&patElectron_relIsoDeltaBeta          ,"patElectron_relIsoDeltaBeta");
    AddBranch(&patElectron_relIsoRhoEA              ,"patElectron_relIsoRhoEA");
    AddBranch(&patElectron_dr03EcalRecHitSumEt      ,"patElectron_dr03EcalRecHitSumEt");
    AddBranch(&patElectron_dr03HcalDepth1TowerSumEt ,"patElectron_dr03HcalDepth1TowerSumEt");
    AddBranch(&patElectron_isolPtTracks             ,"patElectron_isolPtTracks");
    AddBranch(&patElectron_ecalPFClusterIso         ,"patElectron_ecalPFClusterIso");
    AddBranch(&patElectron_hcalPFClusterIso         ,"patElectron_hcalPFClusterIso");
    //Shape, Track related variables, other prop
    AddBranch(&patElectron_dEtaIn                ,"patElectron_dEtaIn");
    AddBranch(&patElectron_dPhiIn                ,"patElectron_dPhiIn");
    AddBranch(&patElectron_full5x5_sigmaIetaIeta ,"patElectron_full5x5_sigmaIetaIeta");
    AddBranch(&patElectron_full5x5_e2x5Max       ,"patElectron_full5x5_e2x5Max");
    AddBranch(&patElectron_full5x5_e5x5          ,"patElectron_full5x5_e5x5");
    AddBranch(&patElectron_full5x5_e1x5          ,"patElectron_full5x5_e1x5");
    AddBranch(&patElectron_hOverE                ,"patElectron_hOverE");
    AddBranch(&patElectron_ooEmooP               ,"patElectron_ooEmooP");
  }
  AddBranch(&passConversionVeto_               ,"patElectron_passConversionVeto");
  if(!_reduced){ 
    AddBranch(&expectedMissingInnerHits          ,"patElectron_expectedMissingInnerHits");
    AddBranch(&patElectron_gsfTrack_ndof         ,"patElectron_gsfTrack_ndof");
    AddBranch(&patElectron_gsfTrack_normChi2     ,"patElectron_gsfTrack_normChi2");
    //IP
    AddBranch(&patElectron_gsfTrack_dz_pv  ,"patElectron_gsfTrack_dz_pv");
    AddBranch(&patElectron_gsfTrack_dxy_pv ,"patElectron_gsfTrack_dxy_pv");
    AddBranch(&patElectron_d0              ,"patElectron_d0");
    AddBranch(&patElectron_dzError         ,"patElectron_dzError");
    AddBranch(&patElectron_dxyError        ,"patElectron_dxyError");
    AddBranch(&patElectron_gsfTrack_vtx    ,"patElectron_gsfTrack_vtx");
    AddBranch(&patElectron_gsfTrack_vty    ,"patElectron_gsfTrack_vty");
    AddBranch(&patElectron_gsfTrack_vtz    ,"patElectron_gsfTrack_vtz");
    //MC info
    if(!_is_data){
      AddBranch(&patElectron_gen_pt                                      ,"patElectron_gen_pt");
      AddBranch(&patElectron_gen_eta                                     ,"patElectron_gen_eta");
      AddBranch(&patElectron_gen_phi                                     ,"patElectron_gen_phi");
      AddBranch(&patElectron_gen_en                                      ,"patElectron_gen_en");
      AddBranch(&patElectron_gen_pdgId                                   ,"patElectron_gen_pdgId");
      AddBranch(&patElectron_gen_isPromptFinalState                      ,"patElectron_gen_isPromptFinalState");
      AddBranch(&patElectron_gen_isDirectPromptTauDecayProductFinalState ,"patElectron_gen_isDirectPromptTauDecayProductFinalState");
    }
  }
  AddBranch(&patElectron_isMatchedToTrigger,"patElectron_isMatchedToTrigger");
  if(debug_) std::cout<<"set branches"<<std::endl;
}
void ElectronPatSelector::Clear(){
  //Kinematics     
  patElectron_pt.clear();
  patElectron_eta.clear();
  patElectron_phi.clear();
  patElectron_energy.clear();
  patElectron_energyCorr.clear();
  patElectron_energySF.clear();
  patElectron_px.clear();
  patElectron_py.clear();
  patElectron_pz.clear();
  patElectron_p.clear();
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
  passMvatrigId_.clear();
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
  patElectron_ecalPFClusterIso.clear();
  patElectron_hcalPFClusterIso.clear();
  //Shape, Track related variables, other prop
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
  //IP
  patElectron_gsfTrack_dz_pv.clear();
  patElectron_gsfTrack_dxy_pv.clear();
  patElectron_d0.clear();
  patElectron_dzError.clear();
  patElectron_dxyError.clear();
  patElectron_gsfTrack_vtx.clear();
  patElectron_gsfTrack_vty.clear();
  patElectron_gsfTrack_vtz.clear();
  //MC info
  if(!_is_data){
    patElectron_gen_pt.clear();
    patElectron_gen_eta.clear();
    patElectron_gen_phi.clear();
    patElectron_gen_en.clear();
    patElectron_gen_pdgId.clear();
    patElectron_gen_isPromptFinalState.clear();
    patElectron_gen_isDirectPromptTauDecayProductFinalState.clear();
  }
  patElectron_isMatchedToTrigger.clear();
}
bool ElectronPatSelector::isGoodVertex(const reco::Vertex& vtx){
  if(vtx.isFake())                                   return false;
  if(vtx.ndof()<_vtx_ndof_min)                       return false;
  if(vtx.position().Rho()>_vtx_rho_max)              return false;
  if(fabs(vtx.position().Z()) > _vtx_position_z_max) return false;
  return true;
}
double ElectronPatSelector::get_effarea(double eta){
  double effarea = -1;
  if(abs(eta) < 1.0)        effarea = 0.1566;
  else if(abs(eta) < 1.479) effarea = 0.1626;
  else if(abs(eta) < 2.0)   effarea = 0.1073;
  else if(abs(eta) < 2.2)   effarea = 0.0854;
  else if(abs(eta) < 2.3)   effarea = 0.1051;
  else if(abs(eta) < 2.4)   effarea = 0.1204;
  else                      effarea = 0.1524;
  return effarea;
}
int ElectronPatSelector::MatchingToTrigger(const edm::Event& iEvent, edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects, edm::Handle<edm::TriggerResults> triggerBits, float eta, float phi){
  const edm::TriggerNames & triggerNames = iEvent.triggerNames(*triggerBits);
  float deltaRMin = 99.;
  for (pat::TriggerObjectStandAlone obj : *triggerObjects) { 
    obj.unpackPathNames(triggerNames);
    if (obj.hasFilterLabel("hltEle115CaloIdVTGsfTrkIdTGsfDphiFilter")) {
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
