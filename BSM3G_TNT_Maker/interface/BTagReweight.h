#ifndef __BTAGREWEIGHT_HE_H_ 
#define __BTAGREWEIGHT_HE_H_
/////
//   Include files and namespaces
/////
#include <memory>
#include <map>
#include <algorithm>
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
#include "DataFormats/PatCandidates/interface/Jet.h"
#include <TTree.h>
#include <Math/VectorUtil.h>
#include "baseTree.h"
using namespace std;
using namespace edm;
/////
//   Class declaration
/////
class BTagReweight : public baseTree{
 public:
  BTagReweight(std::string name, TTree* tree, bool debug, const edm::ParameterSet& cfg);
  ~BTagReweight();
  void Fill(const edm::Event& iEvent);
  void SetBranches();
 private:
  BTagReweight(){};
  /////
  //   Config variables
  /////
  double bWeight;
  edm::InputTag jetToken_;
  /////
  //   IHEP methods/variables
  /////
  edm::FileInPath BTAGReweightfile1_;
  edm::FileInPath BTAGReweightfile2_;
  void fillCSVhistos(TFile *fileHF, TFile *fileLF);
  double get_csv_wgt( std::vector<double> jetPts, std::vector<double> jetEtas, std::vector<double> jetCSVs, std::vector<int> jetFlavors, int iSys, double &csvWgtHF, double &csvWgtLF, double &csvWgtCF );
  // CSV reweighting
  TH1D* h_csv_wgt_hf[9][6];
  TH1D* c_csv_wgt_hf[9][6];
  TH1D* h_csv_wgt_lf[9][4][3];
};
#endif
