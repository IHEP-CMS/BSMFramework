// 
//  Authors:  Andres Florez: Universidad de los Andes, Colombia. 
//  kaur amandeepkalsi: Panjab University, India. 
// 

#ifndef __TRIGGER_H_ 

#define __TRIGGER_H_

#include <memory>

// user include files                                                                      
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
  vector <int> Trigger_decision;
  vector <string> Trigger_names;
  HLTConfigProvider hltConfig_;
  edm::InputTag triggerResultsTag_;
  int triggerSL,triggerDL;
  int HLT_Ele105_CaloIdVT_GsfTrkIdT;
  int HLT_Ele27_eta2p1_WP75_Gsf;
  int HLT_Ele27_WP85_Gsf;
  int HLT_Ele27_eta2p1_WPLoose_Gsf;
  int HLT_Mu50;
  int HLT_IsoMu20;
  int HLT_IsoMu17_eta2p1;
  int HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ;
  int HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL;
  int HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL;
  int HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ;
  int HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ;
  int HLT_IsoMu24_eta2p1;
  int HLT_IsoMu18;
};

#endif

