#ifndef __TOPSUBJET_MU_H_
#define __TOPSUBJET_MU_H_
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
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
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
#include "DataFormats/BTauReco/interface/CATopJetTagInfo.h"
using namespace std;
using namespace pat;
using namespace edm;
/////
//   Class declaration
/////
class TopSubJetSelector : public  baseTree{
 public:
  TopSubJetSelector(std::string name, TTree* tree, bool debug, const edm::ParameterSet& cfg);
  ~TopSubJetSelector();
  void Fill(const edm::Event& iEvent);
  void JECInitialization();
  void SetBranches();
  void Clear();
 private:
  TopSubJetSelector(){};
  /////
  //   Config variables
  /////
  edm::InputTag topsubjetToken_;
  edm::InputTag _vertexInputTag;
  edm::FileInPath jecPayloadNamesAK4PFchsMC1_;
  edm::FileInPath jecPayloadNamesAK4PFchsMC2_;
  edm::FileInPath jecPayloadNamesAK4PFchsMC3_;
  edm::FileInPath jecPayloadNamesAK4PFchsMCUnc_;
  edm::FileInPath jecPayloadNamesAK4PFchsDATA1_;
  edm::FileInPath jecPayloadNamesAK4PFchsDATA2_;
  edm::FileInPath jecPayloadNamesAK4PFchsDATA3_;
  edm::FileInPath jecPayloadNamesAK4PFchsDATAUnc_;
  double _Jet_pt_min;
  bool _is_data;
  /////
  //   JEC
  /////
  boost::shared_ptr<FactorizedJetCorrector>   jecAK4PFchsMC_;
  boost::shared_ptr<JetCorrectionUncertainty> jecAK4PFchsMCUnc_;
  boost::shared_ptr<FactorizedJetCorrector>   jecAK4PFchsDATA_;
  boost::shared_ptr<JetCorrectionUncertainty> jecAK4PFchsDATAUnc_;
  /////
  //   IHEP methods/variables
  /////
  vector <double> TopSubjet_pt, TopSubjet_eta, TopSubjet_phi, TopSubjet_energy, TopSubjet_mass, TopSubjet_Btag0, TopSubjet_Btag1, TopSubjet_Btag2;
  //Jet Energy Corrections and Uncertainties
  vector<double> TopSubjet_JesSF, TopSubjet_JesSFup, TopSubjet_JesSFdown;
};
#endif
