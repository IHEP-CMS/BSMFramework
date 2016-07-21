#ifndef __BOOSTEDJET_MU_H_
#define __BOOSTEDJET_MU_H_
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
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "JetMETCorrections/Modules/interface/JetResolution.h"
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
#include "FWCore/Framework/interface/ConsumesCollector.h"
using namespace std;
using namespace pat;
using namespace edm;
/////
//   Class declaration
/////
class BoostedJetSelector : public  baseTree{
 public:
  BoostedJetSelector(std::string name, TTree* tree, bool debug, const edm::ParameterSet& cfg, edm::ConsumesCollector && ic);
  ~BoostedJetSelector();
  void Fill(const edm::Event& iEvent);
  void SetBranches();
  void JECInitialization();
  void Clear();
  void GetJER(pat::Jet jet, float JesSF, float rhoJER, bool AK8PFchs, float &JERScaleFactor, float &JERScaleFactorUP, float &JERScaleFactorDOWN);
 private:
  BoostedJetSelector(){};
  /////
  //   Config variables
  /////
  edm::EDGetTokenT<reco::VertexCollection> vtx_h_;
  edm::EDGetTokenT<pat::JetCollection> fatjets_;
  edm::EDGetTokenT<double> rhopogHandle_;
  edm::EDGetTokenT<double> rhoJERHandle_;
  edm::FileInPath jecPayloadNamesAK8PFchsMC1_;
  edm::FileInPath jecPayloadNamesAK8PFchsMC2_;
  edm::FileInPath jecPayloadNamesAK8PFchsMC3_;
  edm::FileInPath jecPayloadNamesAK8PFchsMCUnc_;
  edm::FileInPath jecPayloadNamesAK8PFchsDATA1_;
  edm::FileInPath jecPayloadNamesAK8PFchsDATA2_;
  edm::FileInPath jecPayloadNamesAK8PFchsDATA3_;
  edm::FileInPath jecPayloadNamesAK8PFchsDATA4_;
  edm::FileInPath jecPayloadNamesAK8PFchsDATAUnc_;
  std::string jerAK8PFchs_;
  std::string jerAK8PFchsSF_;
  std::string jerAK8PFPuppi_;
  std::string jerAK8PFPuppiSF_;
  bool _is_data;
  /////
  //   JEC
  /////
  boost::shared_ptr<FactorizedJetCorrector>   jecAK8PFchsMC_;
  boost::shared_ptr<JetCorrectionUncertainty> jecAK8PFchsMCUnc_;
  boost::shared_ptr<FactorizedJetCorrector>   jecAK8PFchsDATA_;
  boost::shared_ptr<JetCorrectionUncertainty> jecAK8PFchsDATAUnc_;
  /////
  //   IHEP methods/variables
  /////
  vector <double> BoostedJet_pt, BoostedJet_eta, BoostedJet_phi, BoostedJet_energy, BoostedJet_mass, BoostedJet_Uncorr_pt, BoostedJet_pfJetProbabilityBJetTags, BoostedJet_pfCombinedMVAV2BJetTags, BoostedJet_pfCombinedInclusiveSecondaryVertexV2BJetTags, BoostedJet_pfCombinedCvsLJetTags, BoostedJet_pfCombinedCvsBJetTags;
  vector <double> BoostedJet_neutralHadEnergyFraction, BoostedJet_neutralEmEmEnergyFraction, BoostedJet_chargedHadronEnergyFraction;
  vector <double> BoostedJet_chargedEmEnergyFraction, BoostedJet_muonEnergyFraction,BoostedJet_electronEnergy, BoostedJet_photonEnergy;
  vector <int>    BoostedJet_numberOfConstituents, BoostedJet_chargedMultiplicity;
  vector <double> BoostedJet_tau1, BoostedJet_tau2, BoostedJet_tau3;
  vector <double> BoostedJet_softdrop_mass, //BoostedJet_trimmed_mass,
                                                                       BoostedJet_pruned_mass;//, BoostedJet_filtered_mass;
  vector <double> TopTagging_topMass, TopTagging_minMass, TopTagging_wMass;
  vector <int>    TopTagging_nSubJets;
  //Jet Energy Corrections and Uncertainties
  vector<double> BoostedJet_JesSF, BoostedJet_JesSFup, BoostedJet_JesSFdown, BoostedJet_JerSF, BoostedJet_JerSFup, BoostedJet_JerSFdown; 
};
#endif
