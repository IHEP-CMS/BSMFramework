#ifndef __TTHF_HE_H_ 
#define __TTHF_HE_H_

// system include files
#include <memory>
#include <map>
#include <algorithm>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include <TTree.h>
#include <Math/VectorUtil.h>
#include "FWCore/Framework/interface/ConsumesCollector.h"

#include "baseTree.h"

using namespace std;
using namespace edm;


class GenHFHadrMatchSelector : public baseTree{

 public:
  GenHFHadrMatchSelector(std::string name, TTree* tree, bool debug, const edm::ParameterSet& cfg, edm::ConsumesCollector && iC);
  ~GenHFHadrMatchSelector();
  void Fill(const edm::Event& iEvent);
  void SetBranches();


 private:
  GenHFHadrMatchSelector(){};
  int ttHFCategory;
	
  // Jets configuration
  const double genJetPtMin_=0.;
  const double genJetAbsEtaMax_=0.;
	
  // Input tags
  const edm::EDGetTokenT<reco::GenJetCollection> genJetsToken_;
	
  const edm::EDGetTokenT<std::vector<int> > genBHadJetIndexToken_;
  const edm::EDGetTokenT<std::vector<int> > genBHadFlavourToken_;
  const edm::EDGetTokenT<std::vector<int> > genBHadFromTopWeakDecayToken_;
  const edm::EDGetTokenT<std::vector<reco::GenParticle> > genBHadPlusMothersToken_;
  const edm::EDGetTokenT<std::vector<std::vector<int> > > genBHadPlusMothersIndicesToken_;
  const edm::EDGetTokenT<std::vector<int> > genBHadIndexToken_;
  const edm::EDGetTokenT<std::vector<int> > genBHadLeptonHadronIndexToken_;
  const edm::EDGetTokenT<std::vector<int> > genBHadLeptonViaTauToken_;
	
  const edm::EDGetTokenT<std::vector<int> > genCHadJetIndexToken_;
  const edm::EDGetTokenT<std::vector<int> > genCHadFlavourToken_;
  const edm::EDGetTokenT<std::vector<int> > genCHadFromTopWeakDecayToken_;
  const edm::EDGetTokenT<std::vector<int> > genCHadBHadronIdToken_;
	
};

#endif

