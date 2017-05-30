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
  bool _is_data;
  /////
  //   TTH
  /////
  //Common
  //Electron
  int HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ;
  //Muon
  int HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ;
  int HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ;
  int HLT_IsoMu20;
  int HLT_IsoMu22;
  int HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL;
  int HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL;
  //TTHbb
  //Electron
  int HLT_Ele115_CaloIdVT_GsfTrkIdT;
  int HLT_Ele27_eta2p1_WP75_Gsf;
  int HLT_Ele27_WP85_Gsf;
  int HLT_Ele27_eta2p1_WPLoose_Gsf;
  int HLT_Ele27_eta2p1_WPTight_Gsf;
  //Muon
  int HLT_Mu45_eta2p1;
  int HLT_Mu50;
  int HLT_TkMu50;
  int HLT_IsoMu17_eta2p1;
  int HLT_IsoMu24_eta2p1;
  int HLT_IsoMu18;
  int HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL;
  int HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL;
  int HLT_IsoMu24;
  int HLT_IsoTkMu24;
  //Cross
  int HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ;
  int HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL;
  int HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL;
  int HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ;
  int HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ;
  //TTHLep
  //Electron
  int HLT_Ele23_WPLoose_Gsf; //Data
  int HLT_Ele23_CaloIdL_TrackIdL_IsoVL; //MC
  int HLT_Ele27_WPTight_Gsf;
  int HLT_Ele25_eta2p1_WPTight_Gsf;
  //Muon
  int HLT_IsoTkMu20;
  int HLT_IsoMu22_eta2p1;
  int HLT_IsoTkMu22_eta2p1;
  int HLT_IsoTkMu22;
  //CrossEle-Mu
  int HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL;
  int HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ;
  int HLT_DiMu9_Ele9_CaloIdL_TrackIdL;
  int HLT_Mu8_DiEle12_CaloIdL_TrackIdL;
  int HLT_TripleMu_12_10_5;
  int HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL;
  //Other
  int HLT_DoubleEle33_CaloIdL_MW;
  int HLT_DoubleEle33_CaloIdL_GsfTrkIdVL;
  //Analysis
  int TTHbb_SL;
  int TTHbb_DL;
  int TTHLep_2Mu;
  int TTHLep_2Ele;
  int TTHLep_MuEle;
  int TTHLep_3L4L;
};
#endif
