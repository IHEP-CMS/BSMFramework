#ifndef __TRIGGER_H_ 
#define __TRIGGER_H_
#include <memory>
/////
//   User include files
/////
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "HLTrigger/Muon/interface/HLTMuonIsoFilter.h"
#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Common/interface/TriggerResultsByName.h"
#include "HLTrigger/HLTcore/interface/HLTEventAnalyzerAOD.h"
#include "DataFormats/HLTReco/interface/TriggerRefsCollections.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/PatCandidates/interface/TriggerPath.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/TriggerAlgorithm.h"
#include <TTree.h>
#include <string>
#include <vector>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include "baseTree.h"
#include <TBranch.h>
//namespaces
using namespace std;
using namespace edm;
class TriggerSelector : public baseTree{
 public:
  TriggerSelector(std::string name, TTree* tree, bool debug, const edm::ParameterSet& cfg);
  ~TriggerSelector();
  virtual void startTrigger (edm::EventSetup const& , edm::Run const &);
  void Fill(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  void SetBranches();
  void Clear();
 private:
  TriggerSelector(){};
  /////
  //   Config variables
  /////
  //vector<int> Trigger_decision;
  //vector <string> Trigger_names;
  HLTConfigProvider hltConfig_;
  edm::InputTag triggerBits_;
  double _maxtriggerversion;
  bool   _reduced;

  //TAU
  int HLT_DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg;
  int HLT_DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg;
  int HLT_DoubleTightChargedIsoPFTau40_Trk1_eta2p1_Reg;
  int HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1;
  int HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTau30_eta2p1_CrossL1;
  int HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTau30_eta2p1_CrossL1;
  int HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1;
  int HLT_IsoMu20_eta2p1_MediumChargedIsoPFTau27_eta2p1_CrossL1;
  int HLT_IsoMu20_eta2p1_TightChargedIsoPFTau27_eta2p1_CrossL1;
  int HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_1pr;
  int HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1;
  int HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET100;
  int HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET90;
  int HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr;

  //dineutrino trigger
  int HLT_PFHT500_PFMET100_PFMHT100_IDTight;
  int HLT_PFHT700_PFMET85_PFMHT85_IDTight;
  int HLT_PFHT800_PFMET75_PFMHT75_IDTight;
  int HLT_CaloJet500_NoJetID;
  int HLT_PFHT1050;
  int HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight;
  int HLT_PFMET120_PFMHT120_IDTight;
  int HLT_PFMETNoMu120_PFMHTNoMu120_IDTight;
  int HLT_PFJet500;
  int HLT_AK8PFJet500;
  int HLT_Ele27_eta2p1_WPTight_Gsf;
  int HLT_Ele27_WPTight_Gsf;
  int HLT_Ele25_eta2p1_WPTight_Gsf;
  int HLT_Ele115_CaloIdVT_GsfTrkIdT;
  int HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf;
  int HLT_DoubleEle33_CaloIdL_MW;
  int HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW;
  int HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ;

  int HLT_IsoMu22;
  int HLT_IsoTkMu22;
  int HLT_IsoMu24;
  int HLT_IsoTkMu24;
  int HLT_IsoMu22_eta2p1;
  int HLT_IsoTkMu22_eta2p1;
  int HLT_Mu50;
  int HLT_TkMu50;
  int HLT_DoubleMu33NoFiltersNoVtx;
  int HLT_DoubleMu23NoFiltersNoVtxDisplaced;
  int HLT_Mu30_TkMu11;
  int HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ;
};
#endif
