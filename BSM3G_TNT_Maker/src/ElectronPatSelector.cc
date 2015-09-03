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
  _patElectronToken    = iConfig.getParameter<edm::InputTag>("patElectrons");
  _patElectron_pt_min  = iConfig.getParameter<double>("patElectron_pt_min");
  _patElectron_eta_max = iConfig.getParameter<double>("patElectron_eta_max");
  _vtx_ndof_min        = iConfig.getParameter<int>("vtx_ndof_min");
  _vtx_rho_max         = iConfig.getParameter<int>("vtx_rho_max");
  _vtx_position_z_max  = iConfig.getParameter<double>("vtx_position_z_max");
  SetBranches();
}
ElectronPatSelector::~ElectronPatSelector(){
  delete tree_;
}
void ElectronPatSelector::Fill(const edm::Event& iEvent){
  Clear();
  /////
  //   Recall collections
  /////  
  edm::Handle<reco::VertexCollection> vtx_h;
  iEvent.getByLabel(_vertexInputTag, vtx_h);
  if(vtx_h->empty()) return;
  edm::Handle<edm::View<pat::Electron> > electron_pat;
  iEvent.getByLabel(_patElectronToken, electron_pat);
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
  /////
  //   Require a good vertex 
  /////
  reco::VertexCollection::const_iterator firstGoodVertex = vtx_h->end();
  for(reco::VertexCollection::const_iterator it = vtx_h->begin(); it != firstGoodVertex; it++){
    isGoodVertex(*it);
    firstGoodVertex = it;
    break;
  }
  if(firstGoodVertex == vtx_h->end()) return;
  /////
  //   Get electron information 
  /////
  int elcoun = 0; //Only electrons in Acceptance
  for(edm::View<pat::Electron>::const_iterator el = electron_pat->begin(); el != electron_pat->end(); el++){
    //Acceptance
    if(el->pt()< _patElectron_pt_min)        continue;
    if(fabs(el->eta())>_patElectron_eta_max) continue;  
    //Kinematics    
    patElectron_pt.push_back(el->pt());
    patElectron_eta.push_back(el->eta());
    patElectron_SCeta.push_back(el->superCluster()->position().eta());
    patElectron_phi.push_back(el->phi());
    patElectron_energy.push_back(el->energy());
    //Charge
    patElectron_charge.push_back(el->charge());
    //ID
    //Look up the ID decision for this electron in
    //the ValueMap object and store it. We need a Ptr object as the key.
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
    //Isolation
    reco::GsfElectron::PflowIsolationVariables pfIso = el->pfIsolationVariables();
    isoChargedHadrons_.push_back( pfIso.sumChargedHadronPt );
    isoNeutralHadrons_.push_back( pfIso.sumNeutralHadronEt );
    isoPhotons_.push_back( pfIso.sumPhotonEt );
    isoPU_.push_back( pfIso.sumPUPt );
    //Relative isolation 
    double SumChargedHadronPt = el->pfIsolationVariables().sumChargedHadronPt;
    double SumNeutralEt       = el->pfIsolationVariables().sumNeutralHadronEt + el->pfIsolationVariables().sumPhotonEt;
    double SumPU              = 0.5*el->pfIsolationVariables().sumPUPt;
    double SumNeutralCorrEt   = std::max( 0.0, SumNeutralEt - SumPU );
    double relIso = (SumChargedHadronPt + SumNeutralCorrEt)/el->pt();
    patElectron_relIso.push_back(relIso);
    //Shape
    double absEleSCeta = fabs(el->superCluster()->position().eta());
    bool inCrack  = (absEleSCeta>1.4442 && absEleSCeta<1.5660);
    double dEtaIn = fabs(el->deltaEtaSuperClusterTrackAtVtx());
    double dPhiIn = fabs(el->deltaPhiSuperClusterTrackAtVtx());
    double full5x5_sigmaIetaIeta = el->full5x5_sigmaIetaIeta();
    double hOverE = el->hcalOverEcal();
    double ooEmooP = -999;
    if(el->ecalEnergy()==0)                   ooEmooP = 1e30;
    else if(!std::isfinite(el->ecalEnergy())) ooEmooP = 1e30;
    else                                      ooEmooP = fabs(1.0/el->ecalEnergy() - el->eSuperClusterOverP()/el->ecalEnergy() );
    patElectron_inCrack.push_back(inCrack);
    patElectron_dEtaIn.push_back(dEtaIn);
    patElectron_dPhiIn.push_back(dPhiIn);
    patElectron_full5x5_sigmaIetaIeta.push_back(full5x5_sigmaIetaIeta);
    patElectron_hOverE.push_back(hOverE);
    patElectron_ooEmooP.push_back(ooEmooP);
    //Vertex compatibility
    patElectron_d0.push_back((-1) * el->gsfTrack()->dxy(firstGoodVertex->position()));
    patElectron_dz.push_back(el->gsfTrack()->dz( firstGoodVertex->position()));
    expectedMissingInnerHits.push_back(el->gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS));
    passConversionVeto_.push_back(el->passConversionVeto());
    //Counter
    elcoun++;
  }
}
void ElectronPatSelector::SetBranches(){
  if(debug_) std::cout<<"setting branches: calling AddBranch of baseTree"<<std::endl;
  //Kinematics     
  AddBranch(&patElectron_pt           ,"patElectron_pt");
  AddBranch(&patElectron_eta          ,"patElectron_eta");
  AddBranch(&patElectron_SCeta        ,"patElectron_SCeta");
  AddBranch(&patElectron_phi          ,"patElectron_phi");
  AddBranch(&patElectron_energy       ,"patElectron_energy");
  //Charge
  AddBranch(&patElectron_charge       ,"patElectron_charge");
  //ID
  AddBranch(&passVetoId_              ,"patElectron_isPassVeto");          
  AddBranch(&passLooseId_             ,"patElectron_isPassLoose");
  AddBranch(&passMediumId_            ,"patElectron_isPassMedium");
  AddBranch(&passTightId_             ,"patElectron_isPassTight");
  AddBranch(&passHEEPId_              ,"patElectron_isPassHEEPId");
  AddBranch(&patElectron_pdgId                   ,"patElectron_pdgId");
  //Isolation
  AddBranch(&isoChargedHadrons_       ,"patElectron_isoChargedHadrons");
  AddBranch(&isoNeutralHadrons_       ,"patElectron_isoNeutralHadrons");
  AddBranch(&isoPhotons_              ,"patElectron_isoPhotons");
  AddBranch(&isoPU_                   ,"patElectron_isoPU");
  //Relative isolation 
  AddBranch(&patElectron_relIso                  ,"patElectron_relIso");
  //Shape
  AddBranch(&patElectron_inCrack                 ,"patElectron_inCrack");
  AddBranch(&patElectron_dEtaIn                  ,"patElectron_dEtaIn");
  AddBranch(&patElectron_dPhiIn                  ,"patElectron_dPhiIn");
  AddBranch(&patElectron_full5x5_sigmaIetaIeta   ,"patElectron_full5x5_sigmaIetaIeta");
  AddBranch(&patElectron_hOverE                  ,"patElectron_hOverE");
  AddBranch(&patElectron_ooEmooP                 ,"patElectron_ooEmooP");
  //Vertex compatibility
  AddBranch(&patElectron_d0           ,"patElectron_d0");
  AddBranch(&patElectron_dz           ,"patElectron_dz");
  AddBranch(&expectedMissingInnerHits ,"patElectron_expectedMissingInnerHits");
  AddBranch(&passConversionVeto_      ,"patElectron_passConversionVeto"); 
  if(debug_) std::cout<<"set branches"<<std::endl;
}
void ElectronPatSelector::Clear(){
  //Kinematics     
  patElectron_pt.clear();
  patElectron_eta.clear();
  patElectron_SCeta.clear();
  patElectron_phi.clear();
  patElectron_energy.clear();
  //Charge
  patElectron_charge.clear(); 
  //ID
  passVetoId_.clear();
  passLooseId_.clear();
  passMediumId_.clear();
  passTightId_.clear();  
  passHEEPId_.clear();
  patElectron_pdgId.clear();
  //Isolation
  isoChargedHadrons_.clear();
  isoNeutralHadrons_.clear();
  isoPhotons_.clear();
  isoPU_.clear();
  //Relative isolation 
  patElectron_relIso.clear();
  //Shape
  patElectron_inCrack.clear();
  patElectron_dEtaIn.clear();
  patElectron_dPhiIn.clear();
  patElectron_full5x5_sigmaIetaIeta.clear();
  patElectron_hOverE.clear();
  patElectron_ooEmooP.clear();
  //Vertex compatibility
  patElectron_d0.clear();
  patElectron_dz.clear();
  expectedMissingInnerHits.clear();
  passConversionVeto_.clear();
}
bool ElectronPatSelector::isGoodVertex(const reco::Vertex& vtx){
  if(vtx.isFake())                                   return false;
  if(vtx.ndof()<_vtx_ndof_min)                       return false;
  if(vtx.position().Rho()>_vtx_rho_max)              return false;
  if(fabs(vtx.position().Z()) > _vtx_position_z_max) return false;
  return true;
}
