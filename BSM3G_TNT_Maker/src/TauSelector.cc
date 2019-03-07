#include "BSMFramework/BSM3G_TNT_Maker/interface/TauSelector.h"
TauSelector::TauSelector(std::string name, TTree* tree, bool debug, const pset& iConfig, edm::ConsumesCollector && ic):
  baseTree(name,tree,debug)
{
  vtx_h_                  = ic.consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"));
  beamSpot_               = ic.consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpot"));
  taus_                   = ic.consumes<edm::View<pat::Tau> >(iConfig.getParameter<edm::InputTag>("taus"));
  pfToken_                = ic.consumes<pat::PackedCandidateCollection>(edm::InputTag("packedPFCandidates"));
  _Tau_pt_min     	  = iConfig.getParameter<double>("Tau_pt_min");
  _Tau_eta_max    	  = iConfig.getParameter<double>("Tau_eta_max");
  _Tau_vtx_ndof_min       = iConfig.getParameter<int>("vtx_ndof_min");
  _Tau_vtx_rho_max        = iConfig.getParameter<int>("vtx_rho_max");
  _Tau_vtx_position_z_max = iConfig.getParameter<double>("vtx_position_z_max");
  _reduced                = iConfig.getParameter<bool>("reduced");
  SetBranches();
}
TauSelector::~TauSelector(){
  delete tree_;
}
void TauSelector::Fill(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  Clear();
  /////
  //   Recall collections
  /////
  edm::Handle<reco::VertexCollection> vtx_h;
  iEvent.getByToken(vtx_h_, vtx_h);
  reco::BeamSpot beamSpot;
  edm::Handle<reco::BeamSpot> beamSpotHandle;
  iEvent.getByToken(beamSpot_, beamSpotHandle);
  edm::Handle<edm::View<pat::Tau> > taus;
  iEvent.getByToken(taus_, taus);
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
  //   Get tau information 
  /////
  for(edm::View<pat::Tau>::const_iterator tau = taus->begin(); tau != taus->end(); tau++){
    //Acceptance 
    if(tau->pt() < _Tau_pt_min) continue;
    if(fabs(tau->eta()) > _Tau_eta_max) continue;
    //if(!(tau->leadChargedHadrCand().isNonnull())) continue;
    //Kinematic
    Tau_pt.push_back(tau->pt());
    Tau_eta.push_back(tau->eta());
    Tau_phi.push_back(tau->phi());
    Tau_energy.push_back(tau->energy());
    Tau_px.push_back(tau->px());
    Tau_py.push_back(tau->py());
    Tau_pz.push_back(tau->pz());
    Tau_p.push_back(tau->p());
    //Charge
    Tau_charge.push_back(tau->charge());
    Tau_leadChargedCandCharge.push_back(tau->leadChargedHadrCand().isNonnull() ? tau->leadChargedHadrCand()->charge() : -999);   
    //Decay mode finding
    //Tau_decayModeFindingOldDMs.push_back(tau->tauID("decayModeFindingOldDMs"));
    Tau_decayModeFinding.push_back(tau->tauID("decayModeFinding"));
    Tau_decayModeFindingNewDMs.push_back(tau->tauID("decayModeFindingNewDMs"));
    //Against Muon
    Tau_againstMuonLoose3.push_back(tau->tauID("againstMuonLoose3"));
    Tau_againstMuonTight3.push_back(tau->tauID("againstMuonTight3"));
    //Against Electron
    Tau_againstElectronVLooseMVA6.push_back(tau->tauID("againstElectronVLooseMVA6"));
    Tau_againstElectronLooseMVA6.push_back(tau->tauID("againstElectronLooseMVA6"));
    Tau_againstElectronMediumMVA6.push_back(tau->tauID("againstElectronMediumMVA6"));
    Tau_againstElectronTightMVA6.push_back(tau->tauID("againstElectronTightMVA6"));
    //Isolation
    Tau_chargedIsoPtSum.push_back(tau->tauID("chargedIsoPtSum"));
    Tau_neutralIsoPtSum.push_back(tau->tauID("neutralIsoPtSum"));
    Tau_puCorrPtSum.push_back(tau->tauID("puCorrPtSum")); 
    //Isolation - 2016
    Tau_byLooseCombinedIsolationDeltaBetaCorr3Hits.push_back(tau->tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits"));
    Tau_byMediumCombinedIsolationDeltaBetaCorr3Hits.push_back(tau->tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits"));
    Tau_byTightCombinedIsolationDeltaBetaCorr3Hits.push_back(tau->tauID("byTightCombinedIsolationDeltaBetaCorr3Hits"));
    Tau_byVLooseIsolationMVArun2v1DBoldDMwLT.push_back(tau->tauID("byVLooseIsolationMVArun2v1DBoldDMwLT"));
    Tau_byLooseIsolationMVArun2v1DBoldDMwLT.push_back(tau->tauID("byLooseIsolationMVArun2v1DBoldDMwLT"));
    Tau_byMediumIsolationMVArun2v1DBoldDMwLT.push_back(tau->tauID("byMediumIsolationMVArun2v1DBoldDMwLT"));
    Tau_byTightIsolationMVArun2v1DBoldDMwLT.push_back(tau->tauID("byTightIsolationMVArun2v1DBoldDMwLT"));
    Tau_byVTightIsolationMVArun2v1DBoldDMwLT.push_back(tau->tauID("byVTightIsolationMVArun2v1DBoldDMwLT"));
    Tau_byVLooseIsolationMVArun2v1DBnewDMwLT.push_back(tau->tauID("byVLooseIsolationMVArun2v1DBnewDMwLT"));
    Tau_byLooseIsolationMVArun2v1DBnewDMwLT.push_back(tau->tauID("byLooseIsolationMVArun2v1DBnewDMwLT"));
    Tau_byMediumIsolationMVArun2v1DBnewDMwLT.push_back(tau->tauID("byMediumIsolationMVArun2v1DBnewDMwLT"));
    Tau_byTightIsolationMVArun2v1DBnewDMwLT.push_back(tau->tauID("byTightIsolationMVArun2v1DBnewDMwLT"));
    Tau_byVTightIsolationMVArun2v1DBnewDMwLT.push_back(tau->tauID("byVTightIsolationMVArun2v1DBnewDMwLT"));
    Tau_byVLooseIsolationMVArun2v1DBdR03oldDMwLT.push_back(tau->tauID("byVLooseIsolationMVArun2v1DBdR03oldDMwLT"));
    Tau_byLooseIsolationMVArun2v1DBdR03oldDMwLT.push_back(tau->tauID("byLooseIsolationMVArun2v1DBdR03oldDMwLT"));
    Tau_byMediumIsolationMVArun2v1DBdR03oldDMwLT.push_back(tau->tauID("byMediumIsolationMVArun2v1DBdR03oldDMwLT"));
    Tau_byTightIsolationMVArun2v1DBdR03oldDMwLT.push_back(tau->tauID("byTightIsolationMVArun2v1DBdR03oldDMwLT"));
    Tau_byVTightIsolationMVArun2v1DBdR03oldDMwLT.push_back(tau->tauID("byVTightIsolationMVArun2v1DBdR03oldDMwLT"));
    //Isolation - 2017
    Tau_byVVLooseIsolationMVArun2017v2DBoldDMwLT2017.push_back(tau->tauID("byVVLooseIsolationMVArun2017v2DBoldDMwLT2017"));
    Tau_byVLooseIsolationMVArun2017v2DBoldDMwLT2017.push_back(tau->tauID("byVLooseIsolationMVArun2017v2DBoldDMwLT2017"));
    Tau_byLooseIsolationMVArun2017v2DBoldDMwLT2017.push_back(tau->tauID("byLooseIsolationMVArun2017v2DBoldDMwLT2017"));
    Tau_byMediumIsolationMVArun2017v2DBoldDMwLT2017.push_back(tau->tauID("byMediumIsolationMVArun2017v2DBoldDMwLT2017"));
    Tau_byTightIsolationMVArun2017v2DBoldDMwLT2017.push_back(tau->tauID("byTightIsolationMVArun2017v2DBoldDMwLT2017"));
    Tau_byVTightIsolationMVArun2017v2DBoldDMwLT2017.push_back(tau->tauID("byVTightIsolationMVArun2017v2DBoldDMwLT2017"));
    Tau_byVVTightIsolationMVArun2017v2DBoldDMwLT2017.push_back(tau->tauID("byVVTightIsolationMVArun2017v2DBoldDMwLT2017"));    
    Tau_byVVLooseIsolationMVArun2017v2DBnewDMwLT2017.push_back(tau->tauID("byVVLooseIsolationMVArun2017v2DBnewDMwLT2017"));
    Tau_byVLooseIsolationMVArun2017v2DBnewDMwLT2017.push_back(tau->tauID("byVLooseIsolationMVArun2017v2DBnewDMwLT2017"));
    Tau_byLooseIsolationMVArun2017v2DBnewDMwLT2017.push_back(tau->tauID("byLooseIsolationMVArun2017v2DBnewDMwLT2017"));
    Tau_byMediumIsolationMVArun2017v2DBnewDMwLT2017.push_back(tau->tauID("byMediumIsolationMVArun2017v2DBnewDMwLT2017"));
    Tau_byTightIsolationMVArun2017v2DBnewDMwLT2017.push_back(tau->tauID("byTightIsolationMVArun2017v2DBnewDMwLT2017"));
    Tau_byVTightIsolationMVArun2017v2DBnewDMwLT2017.push_back(tau->tauID("byVTightIsolationMVArun2017v2DBnewDMwLT2017"));
    Tau_byVVTightIsolationMVArun2017v2DBnewDMwLT2017.push_back(tau->tauID("byVVTightIsolationMVArun2017v2DBnewDMwLT2017"));    
    Tau_byVVLooseIsolationMVArun2017v2DBoldDMdR0p3wLT2017.push_back(tau->tauID("byVVLooseIsolationMVArun2017v2DBoldDMdR0p3wLT2017"));
    Tau_byVLooseIsolationMVArun2017v2DBoldDMdR0p3wLT2017.push_back(tau->tauID("byVLooseIsolationMVArun2017v2DBoldDMdR0p3wLT2017"));
    Tau_byLooseIsolationMVArun2017v2DBoldDMdR0p3wLT2017.push_back(tau->tauID("byLooseIsolationMVArun2017v2DBoldDMdR0p3wLT2017"));
    Tau_byMediumIsolationMVArun2017v2DBoldDMdR0p3wLT2017.push_back(tau->tauID("byMediumIsolationMVArun2017v2DBoldDMdR0p3wLT2017"));
    Tau_byTightIsolationMVArun2017v2DBoldDMdR0p3wLT2017.push_back(tau->tauID("byTightIsolationMVArun2017v2DBoldDMdR0p3wLT2017"));
    Tau_byVTightIsolationMVArun2017v2DBoldDMdR0p3wLT2017.push_back(tau->tauID("byVTightIsolationMVArun2017v2DBoldDMdR0p3wLT2017"));
    Tau_byVVTightIsolationMVArun2017v2DBoldDMdR0p3wLT2017.push_back(tau->tauID("byVVTightIsolationMVArun2017v2DBoldDMdR0p3wLT2017"));
    const reco::Track *leadTrack = 0;
    bool isBestTrackNonNull = false;
    bool leadPackedCandidateExists = false;
    if(tau->leadChargedHadrCand().isNonnull()){
      const reco::CandidatePtr hadTauLeadChargedCand = tau->leadChargedHadrCand();
      Tau_leadChargedCandPt.push_back(hadTauLeadChargedCand.isNonnull()  ? hadTauLeadChargedCand->pt()     : -999);      
      Tau_leadChargedCandEta.push_back(hadTauLeadChargedCand.isNonnull() ? hadTauLeadChargedCand->eta()    : -999);
      Tau_leadChargedCandPhi.push_back(hadTauLeadChargedCand.isNonnull() ? hadTauLeadChargedCand->phi()    : -999);
      Tau_leadChargedCandE.push_back(hadTauLeadChargedCand.isNonnull()   ? hadTauLeadChargedCand->energy() : -999);
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
        Tau_leadChargedCandTrack_pt.push_back(leadTrack->pt());
        Tau_leadChargedCandTrack_ptError.push_back(leadTrack->ptError());
      }else{
        Tau_leadChargedCandTrack_pt.push_back(-998);
        Tau_leadChargedCandTrack_ptError.push_back(-998);
      }
    }else{
      Tau_leadChargedCandPt.push_back(-999);  
      Tau_leadChargedCandEta.push_back(-999);  
      Tau_leadChargedCandPhi.push_back(-999);  
      Tau_leadChargedCandE.push_back(-999);  
      Tau_leadChargedCandTrack_pt.push_back(-999);  
      Tau_leadChargedCandTrack_ptError.push_back(-999);  
    }
    //Other prop and Track related variables
    Tau_nProngs.push_back(tau->signalChargedHadrCands().size());
    if(isBestTrackNonNull && leadPackedCandidateExists){
      Tau_leadChargedCandNdof.push_back(leadTrack->ndof());
      Tau_leadChargedCandChi2.push_back(leadTrack->chi2());
      Tau_leadChargedCandValidHits.push_back(leadTrack->numberOfValidHits());
    }else{
      Tau_leadChargedCandNdof.push_back(-998);
      Tau_leadChargedCandChi2.push_back(-998);
      Tau_leadChargedCandValidHits.push_back(-998);
    }
  }
}

void TauSelector::SetBranches(){
  if(debug_) std::cout<<"setting branches: calling AddBranch of baseTree"<<std::endl;
  //Kinematic
  AddBranch(&Tau_pt                           ,"Tau_pt");
  AddBranch(&Tau_eta                          ,"Tau_eta");
  AddBranch(&Tau_phi                          ,"Tau_phi");
  AddBranch(&Tau_energy                       ,"Tau_energy");
  if(!_reduced){
    AddBranch(&Tau_px                           ,"Tau_px");
    AddBranch(&Tau_py                           ,"Tau_py");
    AddBranch(&Tau_pz                           ,"Tau_pz");
    AddBranch(&Tau_p                            ,"Tau_p");
    AddBranch(&Tau_leadChargedCandPt            ,"Tau_leadChargedCandPt");
    AddBranch(&Tau_leadChargedCandEta           ,"Tau_leadChargedCandEta");
    AddBranch(&Tau_leadChargedCandPhi           ,"Tau_leadChargedCandPhi");
    AddBranch(&Tau_leadChargedCandE             ,"Tau_leadChargedCandE");
    AddBranch(&Tau_leadChargedCandTrack_pt      ,"Tau_leadChargedCandTrack_pt");
    AddBranch(&Tau_leadChargedCandTrack_ptError ,"Tau_leadChargedCandTrack_ptError");
  }
  //Charge
  AddBranch(&Tau_charge                ,"Tau_charge");
  if(!_reduced) AddBranch(&Tau_leadChargedCandCharge ,"Tau_leadChargedCandCharge");
  //Decay mode finding
  //AddBranch(&Tau_decayModeFindingOldDMs ,"Tau_decayModeFindingOldDMs");
  AddBranch(&Tau_decayModeFinding       ,"Tau_decayModeFinding");
  AddBranch(&Tau_decayModeFindingNewDMs ,"Tau_decayModeFindingNewDMs");
  //Against Muon
  AddBranch(&Tau_againstMuonLoose3 ,"Tau_againstMuonLoose3");
  AddBranch(&Tau_againstMuonTight3 ,"Tau_againstMuonTight3");
  //Against Electron
  AddBranch(&Tau_againstElectronVLooseMVA6   ,"Tau_againstElectronVLooseMVA6");
  AddBranch(&Tau_againstElectronLooseMVA6    ,"Tau_againstElectronLooseMVA6");
  AddBranch(&Tau_againstElectronMediumMVA6   ,"Tau_againstElectronMediumMVA6");
  AddBranch(&Tau_againstElectronTightMVA6    ,"Tau_againstElectronTightMVA6");
  //Isolation
  AddBranch(&Tau_byLooseCombinedIsolationDeltaBetaCorr3Hits        ,"Tau_byLooseCombinedIsolationDeltaBetaCorr3Hits");
  AddBranch(&Tau_byMediumCombinedIsolationDeltaBetaCorr3Hits       ,"Tau_byMediumCombinedIsolationDeltaBetaCorr3Hits");
  AddBranch(&Tau_byTightCombinedIsolationDeltaBetaCorr3Hits        ,"Tau_byTightCombinedIsolationDeltaBetaCorr3Hits");
  AddBranch(&Tau_byVLooseIsolationMVArun2v1DBoldDMwLT              ,"Tau_byVLooseIsolationMVArun2v1DBoldDMwLT");
  AddBranch(&Tau_byLooseIsolationMVArun2v1DBoldDMwLT               ,"Tau_byLooseIsolationMVArun2v1DBoldDMwLT");
  AddBranch(&Tau_byMediumIsolationMVArun2v1DBoldDMwLT              ,"Tau_byMediumIsolationMVArun2v1DBoldDMwLT");
  AddBranch(&Tau_byTightIsolationMVArun2v1DBoldDMwLT               ,"Tau_byTightIsolationMVArun2v1DBoldDMwLT");
  AddBranch(&Tau_byVTightIsolationMVArun2v1DBoldDMwLT              ,"Tau_byVTightIsolationMVArun2v1DBoldDMwLT");
  AddBranch(&Tau_byVLooseIsolationMVArun2v1DBnewDMwLT              ,"Tau_byVLooseIsolationMVArun2v1DBnewDMwLT");
  AddBranch(&Tau_byLooseIsolationMVArun2v1DBnewDMwLT               ,"Tau_byLooseIsolationMVArun2v1DBnewDMwLT");
  AddBranch(&Tau_byMediumIsolationMVArun2v1DBnewDMwLT              ,"Tau_byMediumIsolationMVArun2v1DBnewDMwLT");
  AddBranch(&Tau_byTightIsolationMVArun2v1DBnewDMwLT               ,"Tau_byTightIsolationMVArun2v1DBnewDMwLT");
  AddBranch(&Tau_byVTightIsolationMVArun2v1DBnewDMwLT              ,"Tau_byVTightIsolationMVArun2v1DBnewDMwLT");
  AddBranch(&Tau_byVLooseIsolationMVArun2v1DBdR03oldDMwLT          ,"Tau_byVLooseIsolationMVArun2v1DBdR03oldDMwLT");
  AddBranch(&Tau_byLooseIsolationMVArun2v1DBdR03oldDMwLT           ,"Tau_byLooseIsolationMVArun2v1DBdR03oldDMwLT");
  AddBranch(&Tau_byMediumIsolationMVArun2v1DBdR03oldDMwLT          ,"Tau_byMediumIsolationMVArun2v1DBdR03oldDMwLT");
  AddBranch(&Tau_byTightIsolationMVArun2v1DBdR03oldDMwLT           ,"Tau_byTightIsolationMVArun2v1DBdR03oldDMwLT");
  AddBranch(&Tau_byVTightIsolationMVArun2v1DBdR03oldDMwLT          ,"Tau_byVTightIsolationMVArun2v1DBdR03oldDMwLT");
  AddBranch(&Tau_byVVLooseIsolationMVArun2017v2DBoldDMwLT2017      ,"byVVLooseIsolationMVArun2017v2DBoldDMwLT2017");
  AddBranch(&Tau_byVLooseIsolationMVArun2017v2DBoldDMwLT2017       ,"byVLooseIsolationMVArun2017v2DBoldDMwLT2017");
  AddBranch(&Tau_byLooseIsolationMVArun2017v2DBoldDMwLT2017        ,"byLooseIsolationMVArun2017v2DBoldDMwLT2017");
  AddBranch(&Tau_byMediumIsolationMVArun2017v2DBoldDMwLT2017       ,"byMediumIsolationMVArun2017v2DBoldDMwLT2017");
  AddBranch(&Tau_byTightIsolationMVArun2017v2DBoldDMwLT2017        ,"byTightIsolationMVArun2017v2DBoldDMwLT2017");
  AddBranch(&Tau_byVTightIsolationMVArun2017v2DBoldDMwLT2017       ,"byVTightIsolationMVArun2017v2DBoldDMwLT2017");
  AddBranch(&Tau_byVVTightIsolationMVArun2017v2DBoldDMwLT2017      ,"byVVTightIsolationMVArun2017v2DBoldDMwLT2017");
  AddBranch(&Tau_byVVLooseIsolationMVArun2017v2DBnewDMwLT2017      ,"byVVLooseIsolationMVArun2017v2DBnewDMwLT2017");
  AddBranch(&Tau_byVLooseIsolationMVArun2017v2DBnewDMwLT2017       ,"byVLooseIsolationMVArun2017v2DBnewDMwLT2017");
  AddBranch(&Tau_byLooseIsolationMVArun2017v2DBnewDMwLT2017        ,"byLooseIsolationMVArun2017v2DBnewDMwLT2017");
  AddBranch(&Tau_byMediumIsolationMVArun2017v2DBnewDMwLT2017       ,"byMediumIsolationMVArun2017v2DBnewDMwLT2017");
  AddBranch(&Tau_byTightIsolationMVArun2017v2DBnewDMwLT2017        ,"byTightIsolationMVArun2017v2DBnewDMwLT2017");
  AddBranch(&Tau_byVTightIsolationMVArun2017v2DBnewDMwLT2017       ,"byVTightIsolationMVArun2017v2DBnewDMwLT2017");
  AddBranch(&Tau_byVVTightIsolationMVArun2017v2DBnewDMwLT2017      ,"byVVTightIsolationMVArun2017v2DBnewDMwLT2017");
  AddBranch(&Tau_byVVLooseIsolationMVArun2017v2DBoldDMdR0p3wLT2017 ,"byVVLooseIsolationMVArun2017v2DBoldDMdR0p3wLT2017");
  AddBranch(&Tau_byVLooseIsolationMVArun2017v2DBoldDMdR0p3wLT2017  ,"byVLooseIsolationMVArun2017v2DBoldDMdR0p3wLT2017");
  AddBranch(&Tau_byLooseIsolationMVArun2017v2DBoldDMdR0p3wLT2017   ,"byLooseIsolationMVArun2017v2DBoldDMdR0p3wLT2017");
  AddBranch(&Tau_byMediumIsolationMVArun2017v2DBoldDMdR0p3wLT2017  ,"byMediumIsolationMVArun2017v2DBoldDMdR0p3wLT2017");
  AddBranch(&Tau_byTightIsolationMVArun2017v2DBoldDMdR0p3wLT2017   ,"byTightIsolationMVArun2017v2DBoldDMdR0p3wLT2017");
  AddBranch(&Tau_byVTightIsolationMVArun2017v2DBoldDMdR0p3wLT2017  ,"byVTightIsolationMVArun2017v2DBoldDMdR0p3wLT2017");
  AddBranch(&Tau_byVVTightIsolationMVArun2017v2DBoldDMdR0p3wLT2017 ,"byVVTightIsolationMVArun2017v2DBoldDMdR0p3wLT2017");
  if(!_reduced){
    AddBranch(&Tau_chargedIsoPtSum                             ,"Tau_chargedIsoPtSum");
    AddBranch(&Tau_neutralIsoPtSum                             ,"Tau_neutralIsoPtSum");
    AddBranch(&Tau_puCorrPtSum                                 ,"Tau_puCorrPtSum");
    //Other prop and Track related variables
    AddBranch(&Tau_nProngs                  ,"Tau_nProngs");
    AddBranch(&Tau_leadChargedCandNdof      ,"Tau_leadChargedCandNdof");
    AddBranch(&Tau_leadChargedCandChi2      ,"Tau_leadChargedCandChi2");
    AddBranch(&Tau_leadChargedCandValidHits ,"Tau_leadChargedCandValidHits");
  }
}

void TauSelector::Clear(){
  //Kinematic
  Tau_pt.clear();
  Tau_eta.clear();
  Tau_phi.clear();
  Tau_energy.clear();
  Tau_px.clear();
  Tau_py.clear();
  Tau_pz.clear();
  Tau_p.clear();
  Tau_leadChargedCandPt.clear();
  Tau_leadChargedCandEta.clear();
  Tau_leadChargedCandPhi.clear();
  Tau_leadChargedCandE.clear();
  Tau_leadChargedCandTrack_pt.clear();
  Tau_leadChargedCandTrack_ptError.clear();
  //Charge
  Tau_charge.clear();
  Tau_leadChargedCandCharge.clear();
  //Decay mode finding
  //Tau_decayModeFindingOldDMs.clear();
  Tau_decayModeFinding.clear();
  Tau_decayModeFindingNewDMs.clear();
  //Against Muon
  Tau_againstMuonLoose3.clear();
  Tau_againstMuonTight3.clear();
  //Against Electron
  Tau_againstElectronVLooseMVA6.clear();
  Tau_againstElectronLooseMVA6.clear();
  Tau_againstElectronMediumMVA6.clear();
  Tau_againstElectronTightMVA6.clear();
  //Isolation
  Tau_byLooseCombinedIsolationDeltaBetaCorr3Hits.clear();
  Tau_byMediumCombinedIsolationDeltaBetaCorr3Hits.clear();
  Tau_byTightCombinedIsolationDeltaBetaCorr3Hits.clear();
  Tau_byVLooseIsolationMVArun2v1DBoldDMwLT.clear();
  Tau_byLooseIsolationMVArun2v1DBoldDMwLT.clear();
  Tau_byMediumIsolationMVArun2v1DBoldDMwLT.clear();
  Tau_byTightIsolationMVArun2v1DBoldDMwLT.clear();
  Tau_byVTightIsolationMVArun2v1DBoldDMwLT.clear();
  Tau_byVLooseIsolationMVArun2v1DBnewDMwLT.clear();
  Tau_byLooseIsolationMVArun2v1DBnewDMwLT.clear();
  Tau_byMediumIsolationMVArun2v1DBnewDMwLT.clear();
  Tau_byTightIsolationMVArun2v1DBnewDMwLT.clear();
  Tau_byVTightIsolationMVArun2v1DBnewDMwLT.clear();
  Tau_byVLooseIsolationMVArun2v1DBdR03oldDMwLT.clear();
  Tau_byLooseIsolationMVArun2v1DBdR03oldDMwLT.clear();
  Tau_byMediumIsolationMVArun2v1DBdR03oldDMwLT.clear();
  Tau_byTightIsolationMVArun2v1DBdR03oldDMwLT.clear();
  Tau_byVTightIsolationMVArun2v1DBdR03oldDMwLT.clear();
  Tau_byVVLooseIsolationMVArun2017v2DBoldDMwLT2017.clear();
  Tau_byVLooseIsolationMVArun2017v2DBoldDMwLT2017.clear();
  Tau_byLooseIsolationMVArun2017v2DBoldDMwLT2017.clear();
  Tau_byMediumIsolationMVArun2017v2DBoldDMwLT2017.clear();
  Tau_byTightIsolationMVArun2017v2DBoldDMwLT2017.clear();
  Tau_byVTightIsolationMVArun2017v2DBoldDMwLT2017.clear();
  Tau_byVVTightIsolationMVArun2017v2DBoldDMwLT2017.clear();
  Tau_byVVLooseIsolationMVArun2017v2DBnewDMwLT2017.clear();
  Tau_byVLooseIsolationMVArun2017v2DBnewDMwLT2017.clear();
  Tau_byLooseIsolationMVArun2017v2DBnewDMwLT2017.clear();
  Tau_byMediumIsolationMVArun2017v2DBnewDMwLT2017.clear();
  Tau_byTightIsolationMVArun2017v2DBnewDMwLT2017.clear();
  Tau_byVTightIsolationMVArun2017v2DBnewDMwLT2017.clear();
  Tau_byVVTightIsolationMVArun2017v2DBnewDMwLT2017.clear();
  Tau_byVVLooseIsolationMVArun2017v2DBoldDMdR0p3wLT2017.clear();
  Tau_byVLooseIsolationMVArun2017v2DBoldDMdR0p3wLT2017.clear();
  Tau_byLooseIsolationMVArun2017v2DBoldDMdR0p3wLT2017.clear();
  Tau_byMediumIsolationMVArun2017v2DBoldDMdR0p3wLT2017.clear();
  Tau_byTightIsolationMVArun2017v2DBoldDMdR0p3wLT2017.clear();
  Tau_byVTightIsolationMVArun2017v2DBoldDMdR0p3wLT2017.clear();
  Tau_byVVTightIsolationMVArun2017v2DBoldDMdR0p3wLT2017.clear();

  Tau_chargedIsoPtSum.clear();
  Tau_neutralIsoPtSum.clear();
  Tau_puCorrPtSum.clear();
  //Other prop and Track related variables
  Tau_nProngs.clear();
  Tau_leadChargedCandNdof.clear();
  Tau_leadChargedCandChi2.clear();
  Tau_leadChargedCandValidHits.clear();
}
bool TauSelector::isGoodVertex(const reco::Vertex& vtxxx) {
  if (vtxxx.isFake()) return false;
  if (vtxxx.ndof() < _Tau_vtx_ndof_min) return false;
  if (vtxxx.position().Rho() > _Tau_vtx_rho_max) return false;
  if (fabs(vtxxx.position().Z()) > _Tau_vtx_position_z_max) return false;
  return true;
}
