#ifndef __PILEUPREWEIGHT_HE_H_ 
#define __PILEUPREWEIGHT_HE_H_
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
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
using namespace std;
using namespace edm;
/////
//   Class declaration
/////
class PileupReweight : public baseTree{
 public:
  PileupReweight(std::string name, TTree* tree, bool debug, const edm::ParameterSet& cfg);
  ~PileupReweight();
  void Fill(const edm::Event& iEvent);
  void SetBranches();
 private:
  PileupReweight(){};
  /////
  //   Config variables
  /////
  bool _is_data;
  edm::FileInPath PUReweightfile_;
  double PUWeight;
  int nPUMax_;
  std::vector<double> puWeigths_;
};
#endif
