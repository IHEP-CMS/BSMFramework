#ifndef __JET_MU_H_
#define __JET_MU_H_
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
#include "FWCore/ParameterSet/interface/FileInPath.h"
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
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "baseTree.h"
using namespace std;
using namespace pat;
using namespace edm;
/////
//   Class declaration
/////
class JetSelector : public  baseTree{
 public:
  JetSelector(std::string name, TTree* tree, bool debug, const edm::ParameterSet& cfg, edm::ConsumesCollector && iC);
  ~JetSelector();
  void Fill(const edm::Event& iEvent);
  void SetBranches();
  void JECInitialization();
  void Clear();
  void GetJER(pat::Jet jet, float JesSF, float rhoJER, bool AK4PFchs, float &JERScaleFactor, float &JERScaleFactorUP, float &JERScaleFactorDOWN);
 private:
  JetSelector(){};
  /////
  //   Config variables
  /////
  edm::EDGetTokenT<reco::VertexCollection> vtx_h_;
  edm::EDGetTokenT<pat::JetCollection> jets_;
  edm::EDGetTokenT<pat::JetCollection> puppijets_;
  edm::EDGetTokenT<double> rhopogHandle_;
  edm::EDGetTokenT<double> rhoJERHandle_;
  edm::FileInPath jecPayloadNamesAK4PFchsMC1_;
  edm::FileInPath jecPayloadNamesAK4PFchsMC2_;
  edm::FileInPath jecPayloadNamesAK4PFchsMC3_;
  edm::FileInPath jecPayloadNamesAK4PFchsMCUnc_;
  edm::FileInPath jecPayloadNamesAK4PFchsDATA1_;
  edm::FileInPath jecPayloadNamesAK4PFchsDATA2_;
  edm::FileInPath jecPayloadNamesAK4PFchsDATA3_;
  edm::FileInPath jecPayloadNamesAK4PFchsDATA4_;
  edm::FileInPath jecPayloadNamesAK4PFchsDATAUnc_;
  edm::FileInPath jecPayloadNamesAK4PFPuppiMC1_;
  edm::FileInPath jecPayloadNamesAK4PFPuppiMC2_;
  edm::FileInPath jecPayloadNamesAK4PFPuppiMC3_;
  edm::FileInPath jecPayloadNamesAK4PFPuppiMCUnc_;
  edm::FileInPath jecPayloadNamesAK4PFPuppiDATA1_;
  edm::FileInPath jecPayloadNamesAK4PFPuppiDATA2_;
  edm::FileInPath jecPayloadNamesAK4PFPuppiDATA3_;
  edm::FileInPath jecPayloadNamesAK4PFPuppiDATA4_;
  edm::FileInPath jecPayloadNamesAK4PFPuppiDATAUnc_;
  std::string jerAK4PFchs_;
  std::string jerAK4PFchsSF_;
  std::string jerAK4PFPuppi_;
  std::string jerAK4PFPuppiSF_;
  double _Jet_pt_min;
  bool _super_TNT;
  bool _PuppiVar;
  bool _is_data;
  /////
  //   JEC
  /////
  boost::shared_ptr<FactorizedJetCorrector>   jecAK4PFchsMC_;
  boost::shared_ptr<JetCorrectionUncertainty> jecAK4PFchsMCUnc_;
  boost::shared_ptr<FactorizedJetCorrector>   jecAK4PFchsDATA_;
  boost::shared_ptr<JetCorrectionUncertainty> jecAK4PFchsDATAUnc_;
  boost::shared_ptr<FactorizedJetCorrector>   jecAK4PFPuppiMC_;
  boost::shared_ptr<JetCorrectionUncertainty> jecAK4PFPuppiMCUnc_;
  boost::shared_ptr<FactorizedJetCorrector>   jecAK4PFPuppiDATA_;
  boost::shared_ptr<JetCorrectionUncertainty> jecAK4PFPuppiDATAUnc_;
  /////
  //   BSM variables
  /////
  ////slimmedJets
  //Kinematics
  vector<double> Jet_pt, Jet_eta, Jet_phi, Jet_energy, Jet_mass, Jet_px, Jet_py, Jet_pz, Jet_Uncorr_pt;
  //ID
  vector<double> Jet_pfCombinedInclusiveSecondaryVertexV2BJetTags, Jet_pfCombinedMVAV2BJetTags, Jet_pfJetProbabilityBJetTags, Jet_pfCombinedCvsLJetTags, Jet_pfCombinedCvsBJetTags, Jet_pileupId, Jet_isPFJet, Jet_isCaloJet;
  //Energy
  vector<double> Jet_neutralHadEnergyFraction, Jet_neutralEmEnergyFraction, Jet_chargedHadronEnergyFraction, Jet_chargedEmEnergyFraction, Jet_muonEnergyFraction, Jet_electronEnergy, Jet_photonEnergy, Jet_emEnergyFraction;
  //Other prop
  vector<double> Jet_numberOfConstituents, Jet_chargedMultiplicity, Jet_vtxMass, Jet_vtxNtracks, Jet_vtx3DVal, Jet_vtx3DSig;
  //Jet Energy Corrections and Uncertainties
  vector<double> Jet_JesSF, Jet_JesSFup, Jet_JesSFdown, Jet_JerSF, Jet_JerSFup, Jet_JerSFdown; 
  //MC 
  vector<double> Jet_partonFlavour, Jet_hadronFlavour;
  ////slimmedJetsPuppi
  //Kinematics
  vector<double> Jet_puppi_pt, Jet_puppi_eta, Jet_puppi_phi, Jet_puppi_energy, Jet_puppi_mass, Jet_puppi_px, Jet_puppi_py, Jet_puppi_pz, Jet_puppi_Uncorr_pt;
  //ID
  vector<double> Jet_puppi_pfCombinedInclusiveSecondaryVertexV2BJetTags, Jet_puppi_pfCombinedMVAV2BJetTags, Jet_puppi_pfJetProbabilityBJetTags, Jet_puppi_pfCombinedCvsLJetTags, Jet_puppi_pfCombinedCvsBJetTags, Jet_puppi_pileupId, Jet_puppi_isPFJet, Jet_puppi_isCaloJet;
  //Energy
  vector<double> Jet_puppi_neutralHadEnergyFraction, Jet_puppi_neutralEmEnergyFraction, Jet_puppi_chargedHadronEnergyFraction, Jet_puppi_chargedEmEnergyFraction, Jet_puppi_muonEnergyFraction, Jet_puppi_electronEnergy, Jet_puppi_photonEnergy, Jet_puppi_emEnergyFraction;
  //Other prop
  vector<double> Jet_puppi_numberOfConstituents, Jet_puppi_chargedMultiplicity, Jet_puppi_vtxMass, Jet_puppi_vtxNtracks, Jet_puppi_vtx3DVal, Jet_puppi_vtx3DSig;
  //Jet Energy Corrections and Uncertainties
  vector<double> Jet_puppi_JesSF, Jet_puppi_JesSFup, Jet_puppi_JesSFdown, Jet_puppi_JerSF, Jet_puppi_JerSFup, Jet_puppi_JerSFdown; 
  //MC 
  vector<double> Jet_puppi_partonFlavour, Jet_puppi_hadronFlavour;
};
#endif
