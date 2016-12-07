#ifndef __EVENTINFO_HE_H_ 
#define __EVENTINFO_HE_H_
/////
//   Include files and namespaces
/////
#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"   
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/LorentzVectorFwd.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/PatCandidates/interface/Isolation.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/normalizedPhi.h"
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/Common/interface/RefVector.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "CLHEP/Random/RandGauss.h"
#include "CommonTools/CandUtils/interface/Booster.h"
#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/METReco/interface/GenMETCollection.h"
#include <Math/VectorUtil.h>
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "TrackingTools/TransientTrack/interface/TrackTransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "CLHEP/Units/GlobalPhysicalConstants.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "CommonTools/ParticleFlow/test/PFIsoReaderDemo.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "boost/regex.hpp"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "PhysicsTools/JetMCUtils/interface/JetMCTag.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "GeneratorInterface/LHEInterface/interface/LHEEvent.h"
#include "SimDataFormats/GeneratorProducts/interface/PdfInfo.h"
#include "LHAPDF/LHAPDF.h"
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TProfile.h>
#include <TTree.h>
#include <string>
#include <vector>
#include <map>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include "baseTree.h"
#include <TRandom3.h>
#include <TBranch.h>
#include "TTree.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
using namespace std;
using namespace edm;
using namespace reco;
using namespace pat;
/////
//   Class declaration
/////
class EventInfoSelector : public baseTree{
 public:
  EventInfoSelector(std::string name, TTree* tree, bool debug, const edm::ParameterSet& cfg, edm::ConsumesCollector && iC);
  ~EventInfoSelector();
  void Fill(const edm::Event& iEvent);
  void SetBranches();
 private:
  EventInfoSelector(){};
  edm::EDGetTokenT<GenEventInfoProduct> genEvtInfo_;
  edm::EDGetTokenT<LHEEventProduct> lheEventProduct_;
  edm::EDGetTokenT<LHEEventProduct> lheEventProductSource_;
  edm::EDGetTokenT<double> rhopogHandle_;
  edm::EDGetTokenT<double> rhotthHandle_;
  edm::EDGetTokenT<double> fixedGridRhoFastjetCentralHandle_;
  edm::EDGetTokenT<double> fixedGridRhoFastjetCentralChargedPileUpHandle_;
  edm::EDGetTokenT<double> fixedGridRhoFastjetCentralNeutralHandle_;
  edm::EDGetTokenT<edm::TriggerResults> metFilterBits_;
  void Initialise();
  //Event quantities
  int EVENT_event_, EVENT_run_, EVENT_lumiBlock_;
  double EVENT_genWeight_, EVENT_genHT;
  bool _is_data; 
  double EVENT_rhopog_, EVENT_rhotth_; 
  double EVENT_originalXWGTUP_, EVENT_scalePDF_;
  double EVENT_PDFtthbbWeightUp_, EVENT_PDFtthbbWeightDown_, EVENT_Q2tthbbWeightUp_, EVENT_Q2tthbbWeightDown_;
  vector<double> EVENT_genWeights_;
  vector<double> EVENT_genWeightsSource_;
  double EVENT_fixedGridRhoFastjetCentral, EVENT_fixedGridRhoFastjetCentralChargedPileUp, EVENT_fixedGridRhoFastjetCentralNeutral;
  //Event filters
  int Flag_HBHENoiseFilter;
  int Flag_HBHENoiseIsoFilter;
  int Flag_CSCTightHaloFilter;
  int Flag_CSCTightHaloTrkMuUnvetoFilter;
  int Flag_CSCTightHalo2015Filter;
  int Flag_HcalStripHaloFilter;
  int Flag_hcalLaserEventFilter;
  int Flag_EcalDeadCellTriggerPrimitiveFilter;
  int Flag_EcalDeadCellBoundaryEnergyFilter;
  int Flag_goodVertices;
  int Flag_eeBadScFilter;
  int Flag_ecalLaserCorrFilter;
  int Flag_trkPOGFilters;
  int Flag_chargedHadronTrackResolutionFilter;
  int Flag_muonBadTrackFilter;
  int Flag_trkPOG_manystripclus53X;
  int Flag_trkPOG_toomanystripclus53X;
  int Flag_trkPOG_logErrorTooManyClusters;
  int Flag_METFilters;
  double BJetness_jetpt0;
  double BJetness_numleps;
  double BJetness_npvTrkOVcollTrk;
  double BJetness_avip3dval;
  double BJetness_avip3dsig;
  double BJetness_avsip3dsig;
  double BJetness_avip1dsig;  
  LHAPDF::PDFSet *read_PDFSet;
  std::vector<LHAPDF::PDF*> _systPDFs;
  bool bjetnessproducer_;
  edm::EDGetTokenT<vector<float> > bjetness_;
};
#endif
