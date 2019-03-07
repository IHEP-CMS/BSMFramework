#include "BSMFramework/BSM3G_TNT_Maker/interface/TauboostedSelector.h"
TauboostedSelector::TauboostedSelector(std::string name, TTree* tree, bool debug, const pset& iConfig, edm::ConsumesCollector && ic):
  baseTree(name,tree,debug)
{
  vtx_h_                  = ic.consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"));
  beamSpot_               = ic.consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpot"));
  tausboosted_            = ic.consumes<edm::View<pat::Tau> >(iConfig.getParameter<edm::InputTag>("tausboosted"));
  pfToken_                = ic.consumes<pat::PackedCandidateCollection>(edm::InputTag("packedPFCandidates"));
  _Tauboosted_pt_min     	  = iConfig.getParameter<double>("Tau_pt_min");
  _Tauboosted_eta_max    	  = iConfig.getParameter<double>("Tau_eta_max");
  _Tauboosted_vtx_ndof_min       = iConfig.getParameter<int>("vtx_ndof_min");
  _Tauboosted_vtx_rho_max        = iConfig.getParameter<int>("vtx_rho_max");
  _Tauboosted_vtx_position_z_max = iConfig.getParameter<double>("vtx_position_z_max");
  _reduced                = iConfig.getParameter<bool>("reduced");
  SetBranches();
}
TauboostedSelector::~TauboostedSelector(){
  delete tree_;
}
void TauboostedSelector::Fill(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  Clear();
  /////
  //   Recall collections
  /////
  edm::Handle<reco::VertexCollection> vtx_h;
  iEvent.getByToken(vtx_h_, vtx_h);
  reco::BeamSpot beamSpot;
  edm::Handle<reco::BeamSpot> beamSpotHandle;
  iEvent.getByToken(beamSpot_, beamSpotHandle);
  edm::Handle<edm::View<pat::Tau> > tausboosted;
  iEvent.getByToken(tausboosted_, tausboosted);
  edm::Handle<pat::PackedCandidateCollection> pfs;
  iEvent.getByToken(pfToken_, pfs);
  edm::ESHandle<TransientTrackBuilder> theB;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);
  /////
  //   Require a good vertex 
  /////
  //reco::VertexCollection::const_iterator firstGoodVertex = vtx->end();
  //for(reco::VertexCollection::const_iterator it = vtx->begin(); it != vtx->end(); it++){
  //  if(isGoodVertex(*it)){
  //    firstGoodVertex = it;
  //    break;
  //  }
  //}
  //if(firstGoodVertex == vtx->end()) return; // skip event if there are no good PVs
  //if(vtx_h->empty()) return; // skip the event if no PV found
  const reco::Vertex &firstGoodVertex = vtx_h->front();
  //bool isgoodvtx = isGoodVertex(firstGoodVertex);
  //if(!isgoodvtx) return;
  GlobalPoint thepv(firstGoodVertex.position().x(),firstGoodVertex.position().y(),firstGoodVertex.position().z());
  /////
  //   Get tauboosted information 
  /////
  for(edm::View<pat::Tau>::const_iterator tauboosted = tausboosted->begin(); tauboosted != tausboosted->end(); tauboosted++){
    //Acceptance 
    if(tauboosted->pt() < _Tauboosted_pt_min) continue;
    if(fabs(tauboosted->eta()) > _Tauboosted_eta_max) continue;
    //if(!(tauboosted->leadChargedHadrCand().isNonnull())) continue;
    //Kinematic
    Tauboosted_pt.push_back(tauboosted->pt());
    Tauboosted_eta.push_back(tauboosted->eta());
    Tauboosted_phi.push_back(tauboosted->phi());
    Tauboosted_energy.push_back(tauboosted->energy());
    Tauboosted_px.push_back(tauboosted->px());
    Tauboosted_py.push_back(tauboosted->py());
    Tauboosted_pz.push_back(tauboosted->pz());
    Tauboosted_p.push_back(tauboosted->p());
    //Charge
    Tauboosted_charge.push_back(tauboosted->charge());
    Tauboosted_leadChargedCandCharge.push_back(tauboosted->leadChargedHadrCand().isNonnull() ? tauboosted->leadChargedHadrCand()->charge() : -999);   
    //Decay mode finding
    //Tauboosted_decayModeFindingOldDMs.push_back(tauboosted->tauID("decayModeFindingOldDMs"));
    Tauboosted_decayModeFinding.push_back(tauboosted->tauID("decayModeFinding"));
    Tauboosted_decayModeFindingNewDMs.push_back(tauboosted->tauID("decayModeFindingNewDMs"));
    //Against Muon
    Tauboosted_againstMuonLoose3.push_back(tauboosted->tauID("againstMuonLoose3"));
    Tauboosted_againstMuonTight3.push_back(tauboosted->tauID("againstMuonTight3"));
    //Against Electron
    Tauboosted_againstElectronVLooseMVA6.push_back(tauboosted->tauID("againstElectronVLooseMVA6"));
    Tauboosted_againstElectronLooseMVA6.push_back(tauboosted->tauID("againstElectronLooseMVA6"));
    Tauboosted_againstElectronMediumMVA6.push_back(tauboosted->tauID("againstElectronMediumMVA6"));
    Tauboosted_againstElectronTightMVA6.push_back(tauboosted->tauID("againstElectronTightMVA6"));
    //Isolation
    Tauboosted_chargedIsoPtSum.push_back(tauboosted->tauID("chargedIsoPtSum"));
    Tauboosted_neutralIsoPtSum.push_back(tauboosted->tauID("neutralIsoPtSum"));
    Tauboosted_puCorrPtSum.push_back(tauboosted->tauID("puCorrPtSum")); 
    Tauboosted_byLooseCombinedIsolationDeltaBetaCorr3Hits.push_back(tauboosted->tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits"));
    Tauboosted_byMediumCombinedIsolationDeltaBetaCorr3Hits.push_back(tauboosted->tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits"));
    Tauboosted_byTightCombinedIsolationDeltaBetaCorr3Hits.push_back(tauboosted->tauID("byTightCombinedIsolationDeltaBetaCorr3Hits"));
    Tauboosted_byVLooseIsolationMVArun2v1DBoldDMwLT.push_back(tauboosted->tauID("byVLooseIsolationMVArun2v1DBoldDMwLT"));
    Tauboosted_byLooseIsolationMVArun2v1DBoldDMwLT.push_back(tauboosted->tauID("byLooseIsolationMVArun2v1DBoldDMwLT"));
    Tauboosted_byMediumIsolationMVArun2v1DBoldDMwLT.push_back(tauboosted->tauID("byMediumIsolationMVArun2v1DBoldDMwLT"));
    Tauboosted_byTightIsolationMVArun2v1DBoldDMwLT.push_back(tauboosted->tauID("byTightIsolationMVArun2v1DBoldDMwLT"));
    Tauboosted_byVTightIsolationMVArun2v1DBoldDMwLT.push_back(tauboosted->tauID("byVTightIsolationMVArun2v1DBoldDMwLT"));
    Tauboosted_byVLooseIsolationMVArun2v1DBnewDMwLT.push_back(tauboosted->tauID("byVLooseIsolationMVArun2v1DBnewDMwLT"));
    Tauboosted_byLooseIsolationMVArun2v1DBnewDMwLT.push_back(tauboosted->tauID("byLooseIsolationMVArun2v1DBnewDMwLT"));
    Tauboosted_byMediumIsolationMVArun2v1DBnewDMwLT.push_back(tauboosted->tauID("byMediumIsolationMVArun2v1DBnewDMwLT"));
    Tauboosted_byTightIsolationMVArun2v1DBnewDMwLT.push_back(tauboosted->tauID("byTightIsolationMVArun2v1DBnewDMwLT"));
    Tauboosted_byVTightIsolationMVArun2v1DBnewDMwLT.push_back(tauboosted->tauID("byVTightIsolationMVArun2v1DBnewDMwLT"));
    //Tauboosted_byLooseCombinedIsolationDeltaBetaCorr3HitsdR03.push_back(tauboosted->tauID("byLooseCombinedIsolationDeltaBetaCorr3HitsdR03"));
    //Tauboosted_byMediumCombinedIsolationDeltaBetaCorr3HitsdR03.push_back(tauboosted->tauID("byMediumCombinedIsolationDeltaBetaCorr3HitsdR03"));
    //Tauboosted_byTightCombinedIsolationDeltaBetaCorr3HitsdR03.push_back(tauboosted->tauID("byTightCombinedIsolationDeltaBetaCorr3HitsdR03"));
    Tauboosted_byVLooseIsolationMVArun2v1DBdR03oldDMwLT.push_back(tauboosted->tauID("byVLooseIsolationMVArun2v1DBdR03oldDMwLT"));
    Tauboosted_byLooseIsolationMVArun2v1DBdR03oldDMwLT.push_back(tauboosted->tauID("byLooseIsolationMVArun2v1DBdR03oldDMwLT"));
    Tauboosted_byMediumIsolationMVArun2v1DBdR03oldDMwLT.push_back(tauboosted->tauID("byMediumIsolationMVArun2v1DBdR03oldDMwLT"));
    Tauboosted_byTightIsolationMVArun2v1DBdR03oldDMwLT.push_back(tauboosted->tauID("byTightIsolationMVArun2v1DBdR03oldDMwLT"));
    Tauboosted_byVTightIsolationMVArun2v1DBdR03oldDMwLT.push_back(tauboosted->tauID("byVTightIsolationMVArun2v1DBdR03oldDMwLT"));
    const reco::Track *leadTrack = 0;
    bool isBestTrackNonNull = false;
    bool leadPackedCandidateExists = false;
    if(tauboosted->leadChargedHadrCand().isNonnull()){
      const reco::CandidatePtr hadTauLeadChargedCand = tauboosted->leadChargedHadrCand();                                                                   
      Tauboosted_leadChargedCandPt.push_back(hadTauLeadChargedCand.isNonnull()  ? hadTauLeadChargedCand->pt()     : -999);      
      Tauboosted_leadChargedCandEta.push_back(hadTauLeadChargedCand.isNonnull() ? hadTauLeadChargedCand->eta()    : -999);
      Tauboosted_leadChargedCandPhi.push_back(hadTauLeadChargedCand.isNonnull() ? hadTauLeadChargedCand->phi()    : -999);
      Tauboosted_leadChargedCandE.push_back(hadTauLeadChargedCand.isNonnull()   ? hadTauLeadChargedCand->energy() : -999);
      //loop over packed PF candidates and find the one which matches the embedded packed candidate within the pat::Tau
      for(unsigned int iPF = 0, nPF = pfs->size(); iPF < nPF; ++iPF){
        const pat::PackedCandidate &pfCandidate = (*pfs)[iPF];
        if((hadTauLeadChargedCand->pt()  == pfCandidate.pt())  &&
           (hadTauLeadChargedCand->eta() == pfCandidate.eta()) && 
           (hadTauLeadChargedCand->phi() == pfCandidate.phi())){ // the packed PF candidate and embedded lead candidate within the pat::Tau should be the same
          leadPackedCandidateExists = true; // if there is a match between the packed PF candidate and embedded lead candidate within the pat::Tau
          if(pfCandidate.bestTrack() != NULL){isBestTrackNonNull = true; leadTrack = pfCandidate.bestTrack();} // grab the associated CTF track (if it exists)
        }
        if(leadPackedCandidateExists && isBestTrackNonNull) break;
      }
      //if(!(leadPackedCandidateExists)) continue; // throw away the tau if there was no matching packed PF candidate to the embedded lead candidate within the pat::Tau
      //if(!(isBestTrackNonNull)) continue; // throw away the tau if it's lead charged hadron has no associated CTF track
      if(isBestTrackNonNull && leadPackedCandidateExists){  
        Tauboosted_leadChargedCandTrack_pt.push_back(leadTrack->pt());
        Tauboosted_leadChargedCandTrack_ptError.push_back(leadTrack->ptError());
      }else{
        Tauboosted_leadChargedCandTrack_pt.push_back(-998);
        Tauboosted_leadChargedCandTrack_ptError.push_back(-998);
      }
    }else{
      Tauboosted_leadChargedCandPt.push_back(-999);  
      Tauboosted_leadChargedCandEta.push_back(-999);  
      Tauboosted_leadChargedCandPhi.push_back(-999);  
      Tauboosted_leadChargedCandE.push_back(-999);  
      Tauboosted_leadChargedCandTrack_pt.push_back(-999);  
      Tauboosted_leadChargedCandTrack_ptError.push_back(-999);  
    }
    //Other prop and Track related variables
    Tauboosted_nProngs.push_back(tauboosted->signalChargedHadrCands().size());
    if(isBestTrackNonNull && leadPackedCandidateExists){
      Tauboosted_leadChargedCandNdof.push_back(leadTrack->ndof());
      Tauboosted_leadChargedCandChi2.push_back(leadTrack->chi2());
      Tauboosted_leadChargedCandValidHits.push_back(leadTrack->numberOfValidHits());
    }else{
      Tauboosted_leadChargedCandNdof.push_back(-998);
      Tauboosted_leadChargedCandChi2.push_back(-998);
      Tauboosted_leadChargedCandValidHits.push_back(-998);
    }
  }
}

void TauboostedSelector::SetBranches(){
  if(debug_) std::cout<<"setting branches: calling AddBranch of baseTree"<<std::endl;
  //Kinematic
  AddBranch(&Tauboosted_pt                           ,"Tauboosted_pt");
  AddBranch(&Tauboosted_eta                          ,"Tauboosted_eta");
  AddBranch(&Tauboosted_phi                          ,"Tauboosted_phi");
  AddBranch(&Tauboosted_energy                       ,"Tauboosted_energy");
  if(!_reduced){
    AddBranch(&Tauboosted_px                           ,"Tauboosted_px");
    AddBranch(&Tauboosted_py                           ,"Tauboosted_py");
    AddBranch(&Tauboosted_pz                           ,"Tauboosted_pz");
    AddBranch(&Tauboosted_p                            ,"Tauboosted_p");
    AddBranch(&Tauboosted_leadChargedCandPt            ,"Tauboosted_leadChargedCandPt");
    AddBranch(&Tauboosted_leadChargedCandEta           ,"Tauboosted_leadChargedCandEta");
    AddBranch(&Tauboosted_leadChargedCandPhi           ,"Tauboosted_leadChargedCandPhi");
    AddBranch(&Tauboosted_leadChargedCandE             ,"Tauboosted_leadChargedCandE");
    AddBranch(&Tauboosted_leadChargedCandTrack_pt      ,"Tauboosted_leadChargedCandTrack_pt");
    AddBranch(&Tauboosted_leadChargedCandTrack_ptError ,"Tauboosted_leadChargedCandTrack_ptError");
  }
  //Charge
  AddBranch(&Tauboosted_charge                ,"Tauboosted_charge");
  if(!_reduced) AddBranch(&Tauboosted_leadChargedCandCharge ,"Tauboosted_leadChargedCandCharge");
  //Decay mode finding
  //AddBranch(&Tauboosted_decayModeFindingOldDMs ,"Tauboosted_decayModeFindingOldDMs");
  AddBranch(&Tauboosted_decayModeFinding       ,"Tauboosted_decayModeFinding");
  AddBranch(&Tauboosted_decayModeFindingNewDMs ,"Tauboosted_decayModeFindingNewDMs");
  //Against Muon
  AddBranch(&Tauboosted_againstMuonLoose3 ,"Tauboosted_againstMuonLoose3");
  AddBranch(&Tauboosted_againstMuonTight3 ,"Tauboosted_againstMuonTight3");
  //Against Electron
  AddBranch(&Tauboosted_againstElectronVLooseMVA6   ,"Tauboosted_againstElectronVLooseMVA6");
  AddBranch(&Tauboosted_againstElectronLooseMVA6    ,"Tauboosted_againstElectronLooseMVA6");
  AddBranch(&Tauboosted_againstElectronMediumMVA6   ,"Tauboosted_againstElectronMediumMVA6");
  AddBranch(&Tauboosted_againstElectronTightMVA6    ,"Tauboosted_againstElectronTightMVA6");
  //Isolation
  AddBranch(&Tauboosted_byLooseCombinedIsolationDeltaBetaCorr3Hits        ,"Tauboosted_byLooseCombinedIsolationDeltaBetaCorr3Hits");
  AddBranch(&Tauboosted_byMediumCombinedIsolationDeltaBetaCorr3Hits       ,"Tauboosted_byMediumCombinedIsolationDeltaBetaCorr3Hits");
  AddBranch(&Tauboosted_byTightCombinedIsolationDeltaBetaCorr3Hits        ,"Tauboosted_byTightCombinedIsolationDeltaBetaCorr3Hits");
  AddBranch(&Tauboosted_byVLooseIsolationMVArun2v1DBoldDMwLT              ,"Tauboosted_byVLooseIsolationMVArun2v1DBoldDMwLT");
  AddBranch(&Tauboosted_byLooseIsolationMVArun2v1DBoldDMwLT               ,"Tauboosted_byLooseIsolationMVArun2v1DBoldDMwLT");
  AddBranch(&Tauboosted_byMediumIsolationMVArun2v1DBoldDMwLT              ,"Tauboosted_byMediumIsolationMVArun2v1DBoldDMwLT");
  AddBranch(&Tauboosted_byTightIsolationMVArun2v1DBoldDMwLT               ,"Tauboosted_byTightIsolationMVArun2v1DBoldDMwLT");
  AddBranch(&Tauboosted_byVTightIsolationMVArun2v1DBoldDMwLT              ,"Tauboosted_byVTightIsolationMVArun2v1DBoldDMwLT");
  AddBranch(&Tauboosted_byVLooseIsolationMVArun2v1DBnewDMwLT              ,"Tauboosted_byVLooseIsolationMVArun2v1DBnewDMwLT");
  AddBranch(&Tauboosted_byLooseIsolationMVArun2v1DBnewDMwLT               ,"Tauboosted_byLooseIsolationMVArun2v1DBnewDMwLT");
  AddBranch(&Tauboosted_byMediumIsolationMVArun2v1DBnewDMwLT              ,"Tauboosted_byMediumIsolationMVArun2v1DBnewDMwLT");
  AddBranch(&Tauboosted_byTightIsolationMVArun2v1DBnewDMwLT               ,"Tauboosted_byTightIsolationMVArun2v1DBnewDMwLT");
  AddBranch(&Tauboosted_byVTightIsolationMVArun2v1DBnewDMwLT              ,"Tauboosted_byVTightIsolationMVArun2v1DBnewDMwLT");
  //AddBranch(&Tauboosted_byLooseCombinedIsolationDeltaBetaCorr3HitsdR03    ,"Tauboosted_byLooseCombinedIsolationDeltaBetaCorr3HitsdR03");
  //AddBranch(&Tauboosted_byMediumCombinedIsolationDeltaBetaCorr3HitsdR03   ,"Tauboosted_byMediumCombinedIsolationDeltaBetaCorr3HitsdR03");
  //AddBranch(&Tauboosted_byTightCombinedIsolationDeltaBetaCorr3HitsdR03    ,"Tauboosted_byTightCombinedIsolationDeltaBetaCorr3HitsdR03");
  AddBranch(&Tauboosted_byVLooseIsolationMVArun2v1DBdR03oldDMwLT          ,"Tauboosted_byVLooseIsolationMVArun2v1DBdR03oldDMwLT");
  AddBranch(&Tauboosted_byLooseIsolationMVArun2v1DBdR03oldDMwLT           ,"Tauboosted_byLooseIsolationMVArun2v1DBdR03oldDMwLT");
  AddBranch(&Tauboosted_byMediumIsolationMVArun2v1DBdR03oldDMwLT          ,"Tauboosted_byMediumIsolationMVArun2v1DBdR03oldDMwLT");
  AddBranch(&Tauboosted_byTightIsolationMVArun2v1DBdR03oldDMwLT           ,"Tauboosted_byTightIsolationMVArun2v1DBdR03oldDMwLT");
  AddBranch(&Tauboosted_byVTightIsolationMVArun2v1DBdR03oldDMwLT          ,"Tauboosted_byVTightIsolationMVArun2v1DBdR03oldDMwLT");
  if(!_reduced){
    AddBranch(&Tauboosted_chargedIsoPtSum                             ,"Tauboosted_chargedIsoPtSum");
    AddBranch(&Tauboosted_neutralIsoPtSum                             ,"Tauboosted_neutralIsoPtSum");
    AddBranch(&Tauboosted_puCorrPtSum                                 ,"Tauboosted_puCorrPtSum");
    //Other prop and Track related variables
    AddBranch(&Tauboosted_nProngs                  ,"Tauboosted_nProngs");
    AddBranch(&Tauboosted_leadChargedCandNdof      ,"Tauboosted_leadChargedCandNdof");
    AddBranch(&Tauboosted_leadChargedCandChi2      ,"Tauboosted_leadChargedCandChi2");
    AddBranch(&Tauboosted_leadChargedCandValidHits ,"Tauboosted_leadChargedCandValidHits");
  }
}

void TauboostedSelector::Clear(){
  //Kinematic
  Tauboosted_pt.clear();
  Tauboosted_eta.clear();
  Tauboosted_phi.clear();
  Tauboosted_energy.clear();
  Tauboosted_px.clear();
  Tauboosted_py.clear();
  Tauboosted_pz.clear();
  Tauboosted_p.clear();
  Tauboosted_leadChargedCandPt.clear();
  Tauboosted_leadChargedCandEta.clear();
  Tauboosted_leadChargedCandPhi.clear();
  Tauboosted_leadChargedCandE.clear();
  Tauboosted_leadChargedCandTrack_pt.clear();
  Tauboosted_leadChargedCandTrack_ptError.clear();
  //Charge
  Tauboosted_charge.clear();
  Tauboosted_leadChargedCandCharge.clear();
  //Decay mode finding
  //Tauboosted_decayModeFindingOldDMs.clear();
  Tauboosted_decayModeFinding.clear();
  Tauboosted_decayModeFindingNewDMs.clear();
  //Against Muon
  Tauboosted_againstMuonLoose3.clear();
  Tauboosted_againstMuonTight3.clear();
  //Against Electron
  Tauboosted_againstElectronVLooseMVA6.clear();
  Tauboosted_againstElectronLooseMVA6.clear();
  Tauboosted_againstElectronMediumMVA6.clear();
  Tauboosted_againstElectronTightMVA6.clear();
  //Isolation
  Tauboosted_byLooseCombinedIsolationDeltaBetaCorr3Hits.clear();
  Tauboosted_byMediumCombinedIsolationDeltaBetaCorr3Hits.clear();
  Tauboosted_byTightCombinedIsolationDeltaBetaCorr3Hits.clear();
  Tauboosted_byVLooseIsolationMVArun2v1DBoldDMwLT.clear();
  Tauboosted_byLooseIsolationMVArun2v1DBoldDMwLT.clear();
  Tauboosted_byMediumIsolationMVArun2v1DBoldDMwLT.clear();
  Tauboosted_byTightIsolationMVArun2v1DBoldDMwLT.clear();
  Tauboosted_byVTightIsolationMVArun2v1DBoldDMwLT.clear();
  Tauboosted_byVLooseIsolationMVArun2v1DBnewDMwLT.clear();
  Tauboosted_byLooseIsolationMVArun2v1DBnewDMwLT.clear();
  Tauboosted_byMediumIsolationMVArun2v1DBnewDMwLT.clear();
  Tauboosted_byTightIsolationMVArun2v1DBnewDMwLT.clear();
  Tauboosted_byVTightIsolationMVArun2v1DBnewDMwLT.clear();
  //Tauboosted_byLooseCombinedIsolationDeltaBetaCorr3HitsdR03.clear();
  //Tauboosted_byMediumCombinedIsolationDeltaBetaCorr3HitsdR03.clear();
  //Tauboosted_byTightCombinedIsolationDeltaBetaCorr3HitsdR03.clear();
  Tauboosted_byVLooseIsolationMVArun2v1DBdR03oldDMwLT.clear();
  Tauboosted_byLooseIsolationMVArun2v1DBdR03oldDMwLT.clear();
  Tauboosted_byMediumIsolationMVArun2v1DBdR03oldDMwLT.clear();
  Tauboosted_byTightIsolationMVArun2v1DBdR03oldDMwLT.clear();
  Tauboosted_byVTightIsolationMVArun2v1DBdR03oldDMwLT.clear();
  Tauboosted_chargedIsoPtSum.clear();
  Tauboosted_neutralIsoPtSum.clear();
  Tauboosted_puCorrPtSum.clear();
  //Other prop and Track related variables
  Tauboosted_nProngs.clear();
  Tauboosted_leadChargedCandNdof.clear();
  Tauboosted_leadChargedCandChi2.clear();
  Tauboosted_leadChargedCandValidHits.clear();
}
bool TauboostedSelector::isGoodVertex(const reco::Vertex& vtxxx) {
  if (vtxxx.isFake()) return false;
  if (vtxxx.ndof() < _Tauboosted_vtx_ndof_min) return false;
  if (vtxxx.position().Rho() > _Tauboosted_vtx_rho_max) return false;
  if (fabs(vtxxx.position().Z()) > _Tauboosted_vtx_position_z_max) return false;
  return true;
}
