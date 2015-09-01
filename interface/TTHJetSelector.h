#ifndef __TTHJET_MU_H_
#define __TTHJET_MU_H_
/////
//   Include files and namespaces
/////
#include <memory>
#include <iostream>
#include <cmath>
#include <vector>
#include <TBranch.h>
#include <TTree.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <string>
#include <map>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <TRandom3.h>
#include <TBranch.h>                                                                    
#include <TClonesArray.h>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerEventWithRefs.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "baseTree.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
using namespace std;
using namespace pat;
using namespace edm;
/////
//   Class declaration
/////
class TTHJetSelector : public baseTree{
 public:
  TTHJetSelector(std::string name, TTree* tree, bool debug, const edm::ParameterSet& cfg, edm::ConsumesCollector && iC);
  ~TTHJetSelector();
  void Fill(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  void SetBranches();
  void Clear();
 private:
  TTHJetSelector(){};
  /////
  //   Config variables
  /////
  edm::InputTag jetToken_;
  edm::InputTag _muonToken;
  edm::InputTag _patElectronToken;
  edm::InputTag _vertexInputTag;
  edm::EDGetTokenT<edm::ValueMap<bool> > electronMediumIdMapToken_;
  bool _super_TNT;
  /////
  //   IHEP methods/variables
  /////
  std::vector<pat::Jet> GetUncorrectedJets(const std::vector<pat::Jet> &inputJets);
  std::vector<pat::Jet> RemoveMuons(const std::vector<pat::Muon>& Muons, const std::vector<pat::Jet>& UncleanedJets, edm::Handle<reco::VertexCollection> vtx_h);
  std::vector<pat::Jet> RemoveElectrons(edm::Handle< vector< pat::Electron > > Electrons, const std::vector<pat::Jet>& UncleanedJets, const edm::Event& iEvent);
  std::vector<pat::Jet> GetCorrectedJets(const std::vector<pat::Jet>& inputJets, const edm::Event& event, const edm::EventSetup& setup);
  void GetJetUncertainties(pat::Jet jet, JetCorrectionUncertainty *jecUnc, float &JesUncertainties);
  vector<double> TTHJet_pt, TTHJet_eta, TTHJet_phi, TTHJet_energy,  TTHJet_bDiscriminator;
  vector<double> TTHJet_bDiscriminator1,  TTHJet_bDiscriminator2;
  vector<double> TTHJet_ptJesUp, TTHJet_ptJesDown;
};
#endif
