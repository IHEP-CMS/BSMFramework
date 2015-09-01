#include "NtupleMaker/BSM3G_TNT_Maker/interface/TTHJetSelector.h"
TTHJetSelector::TTHJetSelector(std::string name, TTree* tree, bool debug, const pset& iConfig, edm::ConsumesCollector && ic):
  baseTree(name,tree,debug),
  electronMediumIdMapToken_(ic.consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electronMediumIdMap")))
{  
  jetToken_         = iConfig.getParameter<edm::InputTag>("jets");
  _vertexInputTag   = iConfig.getParameter<edm::InputTag>("vertices");
  _muonToken        = iConfig.getParameter<edm::InputTag>("muons");
  _patElectronToken = iConfig.getParameter<edm::InputTag>("patElectrons");
  _super_TNT        = iConfig.getParameter<bool>("super_TNT");
  SetBranches();
}
TTHJetSelector::~TTHJetSelector(){
  delete tree_;
}
void TTHJetSelector::Fill(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  Clear();
  /////
  //   Recall collections
  /////	
  edm::Handle<reco::VertexCollection> vtx_h;
  iEvent.getByLabel(_vertexInputTag, vtx_h);
  edm::Handle<pat::MuonCollection> muon_h;                                       
  iEvent.getByLabel(_muonToken, muon_h);  
  edm::Handle< vector< pat::Electron > > electron_pat;
  iEvent.getByLabel(_patElectronToken, electron_pat);
  edm::Handle<pat::JetCollection> jets;                                       
  iEvent.getByLabel(jetToken_, jets);  
  JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty("/afs/cern.ch/work/f/fromeo/CMSSW_7_4_1/src/NtupleMaker/BSM3G_TNT_Maker/textfiles/Summer13_V5_MC_Uncertainty_AK5PFchs.txt");
  /////
  //   Prepare jet collections
  /////		
  std::vector<pat::Jet> jetsNoCorr = GetUncorrectedJets(*jets);
  std::vector<pat::Jet> jetsNoMu   = RemoveMuons(*muon_h, jetsNoCorr, vtx_h);
  std::vector<pat::Jet> jetsNoEl   = RemoveElectrons(electron_pat, jetsNoMu, iEvent);
  std::vector<pat::Jet> jetCorr    = GetCorrectedJets(jetsNoEl, iEvent, iSetup);
  /////
  //   Get TTHJet information
  /////
  for(std::vector<pat::Jet>::const_iterator j = jetCorr.begin(); j != jetCorr.end(); j++){
    //Acceptance
    if(j->pt() < 30)       continue;
    if(fabs(j->eta())>2.4) continue;
    //ID requirements
    if(j->neutralHadronEnergyFraction() >= 0.99) continue;
    if(j->chargedEmEnergyFraction()     >= 0.99) continue;
    if(j->neutralEmEnergyFraction()     >= 0.99) continue;
    if(j->numberOfDaughters()           <= 1)    continue;
    if(j->chargedHadronEnergyFraction() <= 0.0)  continue;
    if(j->chargedMultiplicity()         <= 0.0)  continue;
    //Kinematics	  
    float JesUncertainties=0;
    GetJetUncertainties(*j, jecUnc, JesUncertainties);
    TTHJet_pt.push_back(j->pt());
    TTHJet_ptJesUp.push_back((1+JesUncertainties));         
    TTHJet_ptJesDown.push_back((1-JesUncertainties));         
    TTHJet_eta.push_back(j->eta());       
    TTHJet_phi.push_back(j->phi());       
    TTHJet_energy.push_back(j->energy());
    //ID
    TTHJet_bDiscriminator.push_back(j->bDiscriminator("combinedSecondaryVertexBJetTags"));
    TTHJet_bDiscriminator1.push_back(j->bDiscriminator("combinedInclusiveSecondaryVertexV2BJetTags"));
    TTHJet_bDiscriminator2.push_back(j->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
  } 
}
void TTHJetSelector::SetBranches(){
  if(debug_) std::cout<<"setting branches: calling AddBranch of baseTree"<<std::endl;
  //Kinematics	  
  AddBranch(&TTHJet_pt,                  "TTHJet_pt");
  AddBranch(&TTHJet_ptJesUp,             "TTHJet_ptJesUp");
  AddBranch(&TTHJet_ptJesDown,           "TTHJet_ptJesDown");
  AddBranch(&TTHJet_eta,                 "TTHJet_eta");
  AddBranch(&TTHJet_phi,                 "TTHJet_phi");
  AddBranch(&TTHJet_energy,              "TTHJet_energy");
  //ID
  AddBranch(&TTHJet_bDiscriminator,      "TTHJet_bDiscriminator");
  AddBranch(&TTHJet_bDiscriminator1,     "TTHJet_bDiscriminator1");
  AddBranch(&TTHJet_bDiscriminator2,     "TTHJet_bDiscriminator2");
  if(debug_) std::cout<<"set branches"<<std::endl;
}
void TTHJetSelector::Clear(){
  //Kinematics	  
  TTHJet_pt.clear();
  TTHJet_ptJesUp.clear();
  TTHJet_ptJesDown.clear();
  TTHJet_eta.clear();
  TTHJet_phi.clear();
  TTHJet_energy.clear();
  //ID
  TTHJet_bDiscriminator.clear();
  TTHJet_bDiscriminator1.clear();
  TTHJet_bDiscriminator2.clear();
}
//Get jets before any correction
std::vector<pat::Jet> TTHJetSelector::GetUncorrectedJets(const std::vector<pat::Jet> &inputJets){
  std::vector<pat::Jet> outputJets;
  outputJets.reserve(inputJets.size());
  for(std::vector<pat::Jet>::const_iterator it = inputJets.begin(), ed = inputJets.end(); it != ed; ++it){
    pat::Jet jet = (*it);
    jet.setP4(it->correctedJet(0).p4());
    outputJets.push_back(jet);
  }
  return outputJets;
}
//Remove muons from jet constituents
std::vector<pat::Jet> TTHJetSelector::RemoveMuons(const std::vector<pat::Muon>& Muons, const std::vector<pat::Jet>& UncleanedJets, edm::Handle<reco::VertexCollection> vtx_h){
  reco::VertexCollection::const_iterator firstGoodVertex = vtx_h->end();
  for(reco::VertexCollection::const_iterator it = vtx_h->begin(); it != firstGoodVertex; it++){
    if(it->isFake()) continue;
    if(it->ndof() < 4) continue;
    if(it->position().Rho() > 2) continue;
    if(fabs(it->position().Z()) > 24) continue;
    firstGoodVertex = it;
    break;
  }

  std::vector<pat::Jet> CleanedJets;
  for(std::vector<pat::Jet>::const_iterator iobj1 = UncleanedJets.begin(); iobj1!=UncleanedJets.end(); ++iobj1 ){
    pat::Jet jet = (*iobj1);
    unsigned int nSources1 = jet.numberOfSourceCandidatePtrs();
    bool hasOverlaps = false;
    std::vector<reco::CandidatePtr> overlaps;
    for(std::vector<pat::Muon>::const_iterator iobj2 = Muons.begin(); iobj2!=Muons.end(); ++iobj2 ){
      double SumChargedHadronPt = iobj2->pfIsolationR04().sumChargedHadronPt;
      double SumNeutralEt       = iobj2->pfIsolationR04().sumNeutralHadronEt + iobj2->pfIsolationR04().sumPhotonEt;
      double SumPU              = 0.5*iobj2->pfIsolationR04().sumPUPt;
      double SumNeutralCorrEt   = std::max( 0.0, SumNeutralEt - SumPU );
      double relIso = (SumChargedHadronPt + SumNeutralCorrEt)/iobj2->pt();
      if(iobj2->pt()<20) continue;
      if(fabs(iobj2->eta())>2.4) continue;
      if(firstGoodVertex == vtx_h->end()) continue;
      if(!(iobj2->isTightMuon(*firstGoodVertex))) continue;
      if(relIso>0.12) continue;
      unsigned int nSources2 = iobj2->numberOfSourceCandidatePtrs();
      for( unsigned int i1=0; i1<nSources1; i1++ ){
	reco::CandidatePtr source1 = jet.sourceCandidatePtr(i1);
	if( !(source1.isNonnull() && source1.isAvailable()) ) continue;
	for( unsigned int i2=0; i2<nSources2; i2++ ){
	  reco::CandidatePtr source2 = iobj2->sourceCandidatePtr(i2);
	  if( !(source2.isNonnull() && source2.isAvailable()) ) continue;
	  if( source1==source2 ){
	    hasOverlaps = true;
	    overlaps.push_back(source2);
	  }
	}
      }
    }//end loop over iobj22

    pat::Jet CleanedJet = jet;
    if( hasOverlaps ){
      math::XYZTLorentzVector original = CleanedJet.p4();
      for( int iOverlap=0; iOverlap<int(overlaps.size()); iOverlap++ ){
	const reco::Candidate & cOverlap = *(overlaps[iOverlap]);
	math::XYZTLorentzVector overlaper = cOverlap.p4();
	original -= overlaper;
      }
      CleanedJet.setP4( original );
    }
    CleanedJets.push_back(CleanedJet);
  }	
  return CleanedJets;
}
//Remove electrons from jet constituents
std::vector<pat::Jet> TTHJetSelector::RemoveElectrons(edm::Handle< vector< pat::Electron > > Electrons, const std::vector<pat::Jet>& UncleanedJets, const edm::Event& iEvent){
  edm::Handle<edm::ValueMap<bool> > medium_id_decisions;
  iEvent.getByToken(electronMediumIdMapToken_,medium_id_decisions);

  std::vector<pat::Jet> CleanedJets;
  for(std::vector<pat::Jet>::const_iterator iobj1 = UncleanedJets.begin(); iobj1!=UncleanedJets.end(); ++iobj1 ){
    pat::Jet jet = (*iobj1);
    unsigned int nSources1 = jet.numberOfSourceCandidatePtrs();
    bool hasOverlaps = false;
    std::vector<reco::CandidatePtr> overlaps;
    for(std::vector<pat::Electron>::const_iterator iobj2 = Electrons->begin(); iobj2!=Electrons->end(); ++iobj2 ){
      const Ptr<pat::Electron> elPtr(Electrons, iobj2 - Electrons->begin() );
      bool isPassMedium = (*medium_id_decisions)[ elPtr ];
      if(iobj2->pt()<20) continue;
      if(fabs(iobj2->eta())>2.4) continue;
      if(abs(iobj2->superCluster()->position().eta())>1.4442 && abs(iobj2->superCluster()->position().eta())<1.5660) continue;
      if(!(iobj2->passConversionVeto())) continue; 
      if(!isPassMedium) continue;
      unsigned int nSources2 = iobj2->numberOfSourceCandidatePtrs();
      for( unsigned int i1=0; i1<nSources1; i1++ ){
	reco::CandidatePtr source1 = jet.sourceCandidatePtr(i1);
	if( !(source1.isNonnull() && source1.isAvailable()) ) continue;
	for( unsigned int i2=0; i2<nSources2; i2++ ){
	  reco::CandidatePtr source2 = iobj2->sourceCandidatePtr(i2);
	  if( !(source2.isNonnull() && source2.isAvailable()) ) continue;
	  if( source1==source2 ){
	    hasOverlaps = true;
	    overlaps.push_back(source2);
	  }
	}
      }
    }//end loop over iobj22

    pat::Jet CleanedJet = jet;
    if( hasOverlaps ){
      math::XYZTLorentzVector original = CleanedJet.p4();
      for( int iOverlap=0; iOverlap<int(overlaps.size()); iOverlap++ ){
	const reco::Candidate & cOverlap = *(overlaps[iOverlap]);
	math::XYZTLorentzVector overlaper = cOverlap.p4();
	original -= overlaper;
      }
      CleanedJet.setP4( original );
    }
    CleanedJets.push_back(CleanedJet);
  }	
  return CleanedJets;
}
//Get the corrected jets
std::vector<pat::Jet> TTHJetSelector::GetCorrectedJets(const std::vector<pat::Jet>& inputJets, const edm::Event& event, const edm::EventSetup& setup){
  std::vector<pat::Jet> outputJets;
  const JetCorrector* corrector = JetCorrector::getJetCorrector( "ak4PFchsL1L2L3", setup ); //Get the jet corrector from the event setup
  for( std::vector<pat::Jet>::const_iterator it = inputJets.begin(), ed = inputJets.end(); it != ed; ++it ){
    pat::Jet jet = (*it);
    double scale = corrector->correction(*it, event, setup);
    jet.scaleEnergy( scale );
    outputJets.push_back(jet);
  }
  return outputJets;
}
//Get the jet uncertainties for the systematics
void TTHJetSelector::GetJetUncertainties(pat::Jet jet, JetCorrectionUncertainty *jecUnc_, float &JesUncertainties){
  jecUnc_->setJetEta(jet.eta());
  jecUnc_->setJetPt(jet.pt()); // here you must use the CORRECTED jet pt
  JesUncertainties = jecUnc_->getUncertainty(true);
  //float uncUp = jecUnc_->getUncertainty(true);
  //JesUp = 1 + uncUp;
  //float uncDown = jecUnc_->getUncertainty(false);
  //JesDown = 1 - uncDown;
}
/*
//Alignement for TTHLep analysis
const reco::Vertex &PV = vtx_h->front();
edm::Handle<double> rhoHandle;
iEvent.getByLabel("fixedGridRhoFastjetAll",rhoHandle);
double rho = *rhoHandle;
cout<<"Event "<<iEvent.id().event()<<endl;
cout<<"Muon"<<endl;
//double much1 = 0;
//double much2 = 0;
//int    muc   = 0;
//for(std::vector<pat::Muon>::const_iterator mu = muon_h->begin(); mu!=muon_h->end(); ++mu ){
// if(muc==0) much1 = mu->charge();
// if(muc==1) much2 = mu->charge();
// muc++;
// if(muc==2) break;
//}
if(1){
for(std::vector<pat::Muon>::const_iterator mu = muon_h->begin(); mu!=muon_h->end(); ++mu ){
if(mu->pt()>5 && mu->innerTrack().isNonnull() && mu->globalTrack().isNonnull()){
cout<<mu->pt()<<setw(20)<<mu->eta()<<setw(20)<<mu->phi()<<endl;
double pfIsoCharged = mu->pfIsolationR03().sumChargedHadronPt;
double pfIsoNeutral = mu->pfIsolationR03().sumNeutralHadronEt + mu->pfIsolationR03().sumPhotonEt;
double Eta = abs(mu->eta());
double EffArea = 9999.;
if(abs(Eta) < 0.8)  EffArea = 0.0913;
else if(abs(Eta) < 1.3)  EffArea = 0.0765;
else if(abs(Eta) < 2.0)  EffArea = 0.0546;
else if(abs(Eta) < 2.2)  EffArea = 0.0728;
else EffArea = 0.1177;
double correction = rho*EffArea;
double pfIsoPUSubtracted = std::max( 0.0, pfIsoNeutral - correction);
double result = (pfIsoCharged + pfIsoPUSubtracted)/mu->pt();
cout<<fabs(mu->innerTrack()->dxy(PV.position()))<<setw(20)<<fabs(mu->innerTrack()->dz(PV.position()))<<setw(20)<<result<<setw(20)<<fabs(mu->dB(pat::Muon::PV3D)/mu->edB(pat::Muon::PV3D))<<endl;
cout<<mu->isPFMuon()<<setw(20)<<mu->isGlobalMuon()<<setw(20)<<mu->innerTrack()->ptError()/mu->innerTrack()->pt()<<endl;
cout<<mu->isTrackerMuon()<<setw(20)<<mu->globalTrack()->normalizedChi2()<<setw(20)<<mu->combinedQuality().chi2LocalPosition<<setw(20)<<mu->combinedQuality().trkKink<<setw(20)<<mu->innerTrack()->validFraction()<<setw(20)<<mu->segmentCompatibility()<<endl;
cout<<result-mu->pfIsolationR03().sumChargedHadronPt/mu->pt()<<setw(20)<<mu->pfIsolationR03().sumChargedHadronPt/mu->pt()<<setw(20)<<endl;
double mindr = 999;
pat::Jet mujet;
const JetCorrector* corrector = JetCorrector::getJetCorrector( "ak4PFchsL1L2L3", iSetup );
for(const pat::Jet &j : *jets){
pat::Jet jet = j.correctedJet(0);
double scale = corrector->correction(jet, iEvent, iSetup);
jet.scaleEnergy( scale );
double dr = deltaR(mu->p4(),jet.p4());
if(dr<mindr){
mindr = dr;
mujet = jet;
}
}
cout<<mindr<<setw(20)<<mu->pt()<<setw(20)<<mujet.pt()<<setw(20)<<min(mu->pt()/mujet.pt(), 1.5)<<setw(20)<<max(mujet.bDiscriminator("combinedInclusiveSecondaryVertexV2BJetTags"), float(0.0))<<setw(20)<<fabs(mu->dB(pat::Muon::PV3D)/mu->edB(pat::Muon::PV3D))<<setw(20)<<fabs(mu->innerTrack()->dxy(PV.position()))<<setw(20)<<fabs(mu->innerTrack()->dz(PV.position()))<<setw(20)<<mu->segmentCompatibility()<<endl;
}
}  
}
cout<<"Electron"<<endl;
//double elech1 = 0;
//double elech2 = 0;
//int    elec   = 0;
//for(std::vector<pat::Electron>::const_iterator ele = electron_pat->begin(); ele!=electron_pat->end(); ++ele ){
// if(elec==0) elech1 = ele->charge();
// if(elec==1) elech2 = ele->charge();
// elec++;
// if(elec==2) break;
//}
if(1){
for(std::vector<pat::Electron>::const_iterator ele = electron_pat->begin(); ele!=electron_pat->end(); ++ele ){
if(ele->pt()>10){
cout<<ele->pt()<<setw(20)<<ele->eta()<<setw(20)<<ele->phi()<<endl;
if(ele->gsfTrack().isAvailable()){
double pfIsoCharged = ele->pfIsolationVariables().sumChargedHadronPt;
double pfIsoNeutral = ele->pfIsolationVariables().sumNeutralHadronEt + ele->pfIsolationVariables().sumPhotonEt;
double Eta = abs(ele->eta());
double EffArea = 9999.;
if(abs(Eta) < 0.8)  EffArea = 0.1013;
else if(abs(Eta) < 1.3)  EffArea = 0.0988;
else if(abs(Eta) < 2.0)  EffArea = 0.0572;
else if(abs(Eta) < 2.2)  EffArea = 0.0842;
else EffArea = 0.1530;
double correction = rho*EffArea;
double pfIsoPUSubtracted = std::max( 0.0, pfIsoNeutral - correction);
double result = (pfIsoCharged + pfIsoPUSubtracted)/ele->pt();
cout<<fabs(ele->gsfTrack()->dxy(PV.position()))<<setw(20)<<fabs(ele->gsfTrack()->dz(PV.position()))<<setw(20)<<result<<setw(20)<<fabs(ele->dB(pat::Electron::PV3D)/ele->edB(pat::Electron::PV3D))<<endl;
cout<<ele->gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS)<<setw(20)<<ele->isGsfCtfScPixChargeConsistent()<<setw(20)<<ele->passConversionVeto()<<setw(20)<<endl;
cout<<result-ele->pfIsolationVariables().sumChargedHadronPt/ele->pt()<<setw(20)<<ele->pfIsolationVariables().sumChargedHadronPt/ele->pt()<<setw(20)<<endl;
double mindr = 999;
pat::Jet elejet;
const JetCorrector* corrector = JetCorrector::getJetCorrector( "ak4PFchsL1L2L3", iSetup );
for(const pat::Jet &j : *jets){
pat::Jet jet = j.correctedJet(0);
double scale = corrector->correction(jet, iEvent, iSetup);
jet.scaleEnergy( scale );
double dr = deltaR(ele->p4(),jet.p4());
if(dr<mindr){
mindr = dr;
elejet = jet;
}
}
cout<<mindr<<setw(20)<<ele->pt()<<setw(20)<<elejet.pt()<<setw(20)<<min(double(ele->pt()/elejet.pt()), double(1.5))<<setw(20)<<max(double(elejet.bDiscriminator("combinedInclusiveSecondaryVertexV2BJetTags")), double(0.0))<<setw(20)<<fabs(ele->dB(pat::Electron::PV3D)/ele->edB(pat::Electron::PV3D))<<setw(20)<<fabs(ele->gsfTrack()->dxy(PV.position()))<<setw(20)<<fabs(ele->gsfTrack()->dz(PV.position()))<<endl;
}
}
}
}
*/
