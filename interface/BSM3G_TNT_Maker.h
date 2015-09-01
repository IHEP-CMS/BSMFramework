//-*- C++ -*-
//
//Package: MiniAOD/BSM3G_TNT_Maker
//Class:   BSM3G_TNT_Maker
// 
//
//Authors: Andres Florez:      Universidad de los Andes, Colombia. 
//         kaur amandeepkalsi: Panjab University, India. 
//
#ifndef  TREE_MAKER_H                                                                                                                          
#define  TREE_MAKER_H
/////
//   User include files
/////
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "NtupleMaker/BSM3G_TNT_Maker/interface/GenParticleSelector.h"
#include "NtupleMaker/BSM3G_TNT_Maker/interface/GenHFHadrMatchSelector.h"
#include "NtupleMaker/BSM3G_TNT_Maker/interface/EventInfoSelector.h"
#include "NtupleMaker/BSM3G_TNT_Maker/interface/TriggerSelector.h"
#include "NtupleMaker/BSM3G_TNT_Maker/interface/PVSelector.h"
#include "NtupleMaker/BSM3G_TNT_Maker/interface/MuonSelector.h"
#include "NtupleMaker/BSM3G_TNT_Maker/interface/ElectronPatSelector.h"
#include "NtupleMaker/BSM3G_TNT_Maker/interface/TauSelector.h"
#include "NtupleMaker/BSM3G_TNT_Maker/interface/JetSelector.h"
#include "NtupleMaker/BSM3G_TNT_Maker/interface/TTHJetSelector.h"
#include "NtupleMaker/BSM3G_TNT_Maker/interface/BoostedJetSelector.h"
#include "NtupleMaker/BSM3G_TNT_Maker/interface/TopSubJetSelector.h"
#include "NtupleMaker/BSM3G_TNT_Maker/interface/BTagReweight.h"
#include "NtupleMaker/BSM3G_TNT_Maker/interface/METSelector.h"
#include "NtupleMaker/BSM3G_TNT_Maker/interface/PhotonSelector.h"
#include "baseTree.h"
#include <TBranch.h>
#include <TTree.h>
#include <TFile.h>
#include <vector>
#ifdef __MAKECINT__
#pragma link C++ class std::vector<std::vector<int> >+;
#pragma link C++ class std::vector<std::vector<std::string> >+;
#pragma link C++ class std::vector<std::vector<TString> >+;
#pragma link C++ class std::vector<std::vector<float> >+;
#pragma link C++ class std::vector<std::vector<double> >+;
#pragma link C++ class std::vector<std::vector<bool> >+;
#pragma extra_include "std::vector";
#endif
/////
//   Class declaration
/////
class BSM3G_TNT_Maker : public edm::EDAnalyzer {
 public:
  explicit BSM3G_TNT_Maker(const edm::ParameterSet&);
  ~BSM3G_TNT_Maker();
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
 private:
  virtual void beginJob() override;
  virtual void beginRun(edm::Run const &, edm::EventSetup const&);
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  /////
  //   Config variables
  /////
  bool _fillgeninfo;
  bool _fillgenHFCategoryinfo;
  bool _filleventinfo;
  bool _filltriggerinfo; 
  bool _fillPVinfo;
  bool _fillmuoninfo;
  bool _fillelectronpatinfo;
  bool _filltauinfo;
  bool _filljetinfo;
  bool _filltthjetinfo;
  bool _fillBoostedJetinfo;
  bool _fillTopSubJetinfo;
  bool _fillBTagReweight;
  bool _fillMETinfo;
  bool _fillphotoninfo;
  /////
  //   BSM methods/variables
  /////
  TFile* file;
  TTree* tree_;
  const size_t MaxN;
  bool debug_;
  GenParticleSelector    *genselector;
  GenHFHadrMatchSelector *genhfselector;
  EventInfoSelector      *eventinfoselector;
  TriggerSelector        *trselector;
  PVSelector             *pvselector;
  MuonSelector           *muselector;
  ElectronPatSelector    *elpatselector;
  TauSelector            *tauselector;
  JetSelector            *jetselector;
  TTHJetSelector         *tthjetselector;
  BoostedJetSelector     *BoostedJetselector;
  TopSubJetSelector      *TopSubJetselector;
  BTagReweight           *btagreweight;
  METSelector            *metselector;
  PhotonSelector         *photonselector;
};
#endif 
