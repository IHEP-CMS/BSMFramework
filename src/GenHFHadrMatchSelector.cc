#include "NtupleMaker/BSM3G_TNT_Maker/interface/GenHFHadrMatchSelector.h"

GenHFHadrMatchSelector::GenHFHadrMatchSelector(std::string name, TTree* tree, bool debug, const pset& iConfig, edm::ConsumesCollector && ic):
  baseTree(name,tree,debug),
  genJetPtMin_(iConfig.getParameter<double>("genJetPtMin")),
  genJetAbsEtaMax_(iConfig.getParameter<double>("genJetAbsEtaMax")),
  genJetsToken_(ic.consumes<reco::GenJetCollection>(iConfig.getParameter<edm::InputTag>("genJets"))),
  genBHadJetIndexToken_(ic.consumes<std::vector<int> >(iConfig.getParameter<edm::InputTag>("genBHadJetIndex"))),
  genBHadFlavourToken_(ic.consumes<std::vector<int> >(iConfig.getParameter<edm::InputTag>("genBHadFlavour"))),
  genBHadFromTopWeakDecayToken_(ic.consumes<std::vector<int> >(iConfig.getParameter<edm::InputTag>("genBHadFromTopWeakDecay"))),
  genBHadPlusMothersToken_(ic.consumes<std::vector<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("genBHadPlusMothers"))),
  genBHadPlusMothersIndicesToken_(ic.consumes<std::vector<std::vector<int> > >(iConfig.getParameter<edm::InputTag>("genBHadPlusMothersIndices"))),
  genBHadIndexToken_(ic.consumes<std::vector<int> >(iConfig.getParameter<edm::InputTag>("genBHadIndex"))),
  genBHadLeptonHadronIndexToken_(ic.consumes<std::vector<int> >(iConfig.getParameter<edm::InputTag>("genBHadLeptonHadronIndex"))),
  genBHadLeptonViaTauToken_(ic.consumes<std::vector<int> >(iConfig.getParameter<edm::InputTag>("genBHadLeptonViaTau"))),
  genCHadJetIndexToken_(ic.consumes<std::vector<int> >(iConfig.getParameter<edm::InputTag>("genCHadJetIndex"))),
  genCHadFlavourToken_(ic.consumes<std::vector<int> >(iConfig.getParameter<edm::InputTag>("genCHadFlavour"))),
  genCHadFromTopWeakDecayToken_(ic.consumes<std::vector<int> >(iConfig.getParameter<edm::InputTag>("genCHadFromTopWeakDecay"))),
  genCHadBHadronIdToken_(ic.consumes<std::vector<int> >(iConfig.getParameter<edm::InputTag>("genCHadBHadronId")))

{
  if(debug) std::cout<<"in GenHFHadrMatchSelector constructor"<<std::endl;
  if(debug) std::cout<<"in GenHFHadrMatchSelector constructor: calling SetBrances()"<<std::endl;
  SetBranches();
}

GenHFHadrMatchSelector::~GenHFHadrMatchSelector(){
  delete tree_;
}

void GenHFHadrMatchSelector::Fill(const edm::Event& iEvent){
  if(debug_)    std::cout<<"getting ttHFCategorization info"<<std::endl;
  // Reading gen jets from the event
  edm::Handle<reco::GenJetCollection> genJets;
  iEvent.getByToken(genJetsToken_, genJets);
    
  // Reading B hadrons related information
  edm::Handle<std::vector<int> > genBHadFlavour;
  iEvent.getByToken(genBHadFlavourToken_, genBHadFlavour);
    
  edm::Handle<std::vector<int> > genBHadJetIndex;
  iEvent.getByToken(genBHadJetIndexToken_, genBHadJetIndex);
    
  edm::Handle<std::vector<int> > genBHadFromTopWeakDecay;
  iEvent.getByToken(genBHadFromTopWeakDecayToken_, genBHadFromTopWeakDecay);
    
  edm::Handle<std::vector<reco::GenParticle> > genBHadPlusMothers;
  iEvent.getByToken(genBHadPlusMothersToken_, genBHadPlusMothers);
    
  edm::Handle<std::vector<std::vector<int> > > genBHadPlusMothersIndices;
  iEvent.getByToken(genBHadPlusMothersIndicesToken_, genBHadPlusMothersIndices);
    
  edm::Handle<std::vector<int> > genBHadIndex;
  iEvent.getByToken(genBHadIndexToken_, genBHadIndex);
    
  edm::Handle<std::vector<int> > genBHadLeptonHadronIndex;
  iEvent.getByToken(genBHadLeptonHadronIndexToken_, genBHadLeptonHadronIndex);
    
  edm::Handle<std::vector<int> > genBHadLeptonViaTau;
  iEvent.getByToken(genBHadLeptonViaTauToken_, genBHadLeptonViaTau);
    
  // Reading C hadrons related information
  edm::Handle<std::vector<int> > genCHadFlavour;
  iEvent.getByToken(genCHadFlavourToken_, genCHadFlavour);
    
  edm::Handle<std::vector<int> > genCHadJetIndex;
  iEvent.getByToken(genCHadJetIndexToken_, genCHadJetIndex);
    
  edm::Handle<std::vector<int> > genCHadFromTopWeakDecay;
  iEvent.getByToken(genCHadFromTopWeakDecayToken_, genCHadFromTopWeakDecay);
    
  edm::Handle<std::vector<int> > genCHadBHadronId;
  iEvent.getByToken(genCHadBHadronIdToken_, genCHadBHadronId);
	
  // Map <jet index, number of specific hadrons in the jet>
  // B jets with b hadrons directly from top quark decay
  std::map<int, int> bJetFromTopIds_all;
  // B jets with b hadrons directly from top quark decay
  std::map<int, int> bJetFromTopIds;
  // B jets with b hadrons after top quark decay
  std::map<int, int> bJetAfterTopIds;
  // B jets with b hadrons before top quark decay chain
  std::map<int, int> bJetBeforeTopIds;
  // C jets with c hadrons before top quark decay chain
  std::map<int, int> cJetBeforeTopIds;
  // C jets with c hadrons after top quark decay
  std::map<int, int> cJetAfterTopIds;
    
  // Counting number of specific hadrons in each b jet
  for(size_t hadronId = 0; hadronId < genBHadIndex->size(); ++hadronId) {
    // Flavour of the hadron's origin
    const int flavour = genBHadFlavour->at(hadronId);
    // Whether hadron radiated before top quark decay
    // const bool fromTopDecay = genBHadFromTopWeakDecay.at(hadronId);
    // Index of a jet associated to the hadron
    const int jetIndex = genBHadJetIndex->at(hadronId);
    // Skipping hadrons which have no associated jet
    if(jetIndex < 0) continue;
    // Jet from direct top quark decay [pdgId(top)=6]
    if(std::abs(flavour) == 6) {
      if(bJetFromTopIds_all.count(jetIndex) < 1) bJetFromTopIds_all[jetIndex] = 1;
      else bJetFromTopIds_all[jetIndex]++;
    }
    // Skipping if jet is not in acceptance
    if(genJets->at(jetIndex).pt() < genJetPtMin_) continue;
    if(std::fabs(genJets->at(jetIndex).eta()) > genJetAbsEtaMax_) continue;
    // Identifying jets with b hadrons not from top quark decay
    // Jet from direct top quark decay [pdgId(top)=6]
    if(std::abs(flavour) == 6) {
      if(bJetFromTopIds.count(jetIndex) < 1) bJetFromTopIds[jetIndex] = 1;
      else bJetFromTopIds[jetIndex]++;
    }
    // Skipping if jet is from top quark decay
    if(std::abs(flavour) == 6) continue;
    // Jet before top quark decay
    // if(!fromTopDecay) {
    if(bJetBeforeTopIds.count(jetIndex) < 1) bJetBeforeTopIds[jetIndex] = 1;
    else bJetBeforeTopIds[jetIndex]++;
    // }
    // // Jet after top quark decay but not directly from top
    // else if(fromTopDecay) {
    //     if(bJetAfterTopIds.count(jetIndex) < 1) bJetAfterTopIds[jetIndex] = 1;
    //     else bJetAfterTopIds[jetIndex]++;
    // }
  }
    
  // Counting number of specific hadrons in each c jet
  for(size_t hadronId = 0; hadronId < genCHadJetIndex->size(); ++hadronId) {
    // Skipping c hadrons that are coming from b hadrons
    if(genCHadBHadronId->at(hadronId) >= 0) continue;
    // Skipping c hadrons coming for W-dcays
    if(abs(genCHadFlavour->at(hadronId))==24) continue;
    // Index of a jet associated to the hadron
    const int jetIndex = genCHadJetIndex->at(hadronId);
    // Whether hadron radiated before top quark decay
    // const bool fromTopDecay = genCHadFromTopWeakDecay.at(hadronId);
    // Skipping hadrons which have no associated jet
    if(jetIndex < 0) continue;
    // Skipping if jet is not in acceptance
    if(genJets->at(jetIndex).pt() < genJetPtMin_) continue;
    if(std::fabs(genJets->at(jetIndex).eta()) > genJetAbsEtaMax_) continue;
    // Jet before top quark decay
    // if(!fromTopDecay) {
    if(cJetBeforeTopIds.count(jetIndex) < 1) cJetBeforeTopIds[jetIndex] = 1;
    else cJetBeforeTopIds[jetIndex]++;
    // }
    // // Jet after top quark decay but not directly from top
    // else if(fromTopDecay) {
    //     if(cJetAfterTopIds.count(jetIndex) < 1) cJetAfterTopIds[jetIndex] = 1;
    //     else cJetAfterTopIds[jetIndex]++;
    // }
  }
    
  // Finding additional b jets (before top decay)
  std::vector<int> additionalBJetIds;
  for(std::map<int, int>::iterator it = bJetBeforeTopIds.begin(); it != bJetBeforeTopIds.end(); ++it) {
    const int jetId = it->first;
    // Skipping the jet if it contains a b hadron directly from top quark decay
    if(bJetFromTopIds.count(jetId) > 0) continue;
    additionalBJetIds.push_back(jetId);
  }
  // Finding pseudo-additional b jets (after top decay)
  std::vector<int> pseudoadditionalBJetIds;
  for(std::map<int, int>::iterator it = bJetAfterTopIds.begin(); it != bJetAfterTopIds.end(); ++it) {
    const int jetId = it->first;
    // Skipping the jet if it contains a b hadron directly from top quark decay
    if(bJetFromTopIds.count(jetId) > 0) continue;
    pseudoadditionalBJetIds.push_back(jetId);
  }
  // Finding additional c jets
  std::vector<int> additionalCJetIds;
  for(std::map<int, int>::iterator it = cJetBeforeTopIds.begin(); it != cJetBeforeTopIds.end(); ++it) {
    const int jetId = it->first;
    if(bJetFromTopIds.count(jetId) > 0) continue;
    additionalCJetIds.push_back(jetId);
  }
  // Finding pseudo-additional c jets (after top decay)
  std::vector<int> pseudoadditionalCJetIds;
  for(std::map<int, int>::iterator it = cJetAfterTopIds.begin(); it != cJetAfterTopIds.end(); ++it) {
    const int jetId = it->first;
    // Skipping the jet if it contains a b hadron directly from top quark decay
    if(bJetFromTopIds.count(jetId) > 0) continue;
    pseudoadditionalCJetIds.push_back(jetId);
  }
    
  // Categorizing event based on number of additional b/c jets 
  // and number of corresponding hadrons in each of them
  // int additionalJetEventId;
  int additionalJetEventId = bJetFromTopIds.size()*100;
  // tt + 1 additional b jet
  if (additionalBJetIds.size() == 1) {
    int nHadronsInJet = bJetBeforeTopIds[additionalBJetIds.at(0)];
    // tt + 1 additional b jet from 1 additional b hadron
    if(nHadronsInJet == 1) additionalJetEventId = 51;
    // tt + 1 additional b jet from >=2 additional b hadrons
    else additionalJetEventId = 52;
  }
  // tt + 2 additional b jets
  else if (additionalBJetIds.size() > 1) {
    int nHadronsInJet1 = bJetBeforeTopIds[additionalBJetIds.at(0)];
    int nHadronsInJet2 = bJetBeforeTopIds[additionalBJetIds.at(1)];
    // tt + 2 additional b jets each from 1 additional b hadron
    if(std::max(nHadronsInJet1, nHadronsInJet2) == 1) additionalJetEventId = 53;
    // tt + 2 additional b jets one of which from >=2 overlapping additional b hadrons
    else if(std::min(nHadronsInJet1, nHadronsInJet2) == 1 && std::max(nHadronsInJet1, nHadronsInJet2) > 1) additionalJetEventId = 54;
    // tt + 2 additional b jets each from >=2 additional b hadrons
    else if(std::min(nHadronsInJet1, nHadronsInJet2) > 1) additionalJetEventId = 55;
  }
  // tt + no additional b jets
  else if(additionalBJetIds.size() == 0) {
    // tt + >=1 pseudo-additional b jet with b hadrons after top quark decay
    if(pseudoadditionalBJetIds.size() > 0) additionalJetEventId = 56;
    // tt + 1 additional c jet
    else if(additionalCJetIds.size() == 1) {
      int nHadronsInJet = cJetBeforeTopIds[additionalCJetIds.at(0)];
      // tt + 1 additional c jet from 1 additional c hadron
      if(nHadronsInJet == 1) additionalJetEventId = 41;
      // tt + 1 additional c jet from >=2 overlapping additional c hadrons
      else additionalJetEventId = 42;
    }
    // tt + >=2 additional c jets
    else if(additionalCJetIds.size() > 1) {
      int nHadronsInJet1 = cJetBeforeTopIds[additionalCJetIds.at(0)];
      int nHadronsInJet2 = cJetBeforeTopIds[additionalCJetIds.at(1)];
      // tt + 2 additional c jets each from 1 additional c hadron
      if(std::max(nHadronsInJet1, nHadronsInJet2) == 1) additionalJetEventId = 43;
      // tt + 2 additional c jets one of which from >=2 overlapping additional c hadrons
      else if(std::min(nHadronsInJet1, nHadronsInJet2) == 1 && std::max(nHadronsInJet1, nHadronsInJet2) > 1) additionalJetEventId = 44;
      // tt + 2 additional c jets each from >=2 additional c hadrons
      else if(std::min(nHadronsInJet1, nHadronsInJet2) > 1) additionalJetEventId = 45;
    }
    // tt + no additional c jets
    else if(additionalCJetIds.size() == 0) {
      // tt + >=1 pseudo-additional c jet with c hadrons after top quark decay
      if(pseudoadditionalCJetIds.size() > 0) additionalJetEventId = 46;
      // tt + light jets
      else additionalJetEventId = 0;
    }
  }
	
  ttHFCategory=additionalJetEventId;
	
  if(debug_)    std::cout<<"got ttHFCategorization info"<<std::endl;
}
void GenHFHadrMatchSelector::SetBranches(){
  if(debug_)    std::cout<<"setting branches: calling AddBranch of ttHFCategorization"<<std::endl;
  AddBranch(&ttHFCategory,"ttHFCategory");
}
