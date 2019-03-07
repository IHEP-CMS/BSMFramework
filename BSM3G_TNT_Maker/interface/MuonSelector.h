#ifndef __MUON_MU_H_                                                                                                            
#define __MUON_MU_H_
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
#include <TClonesArray.h>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "baseTree.h"
#include "TLorentzVector.h"
//Track builder infos
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexUpdator.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexTrackCompatibilityEstimator.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexSmoother.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/TrackExtra.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "RecoBTag/BTagTools/interface/SignedTransverseImpactParameter.h"
#include "RecoJets/JetProducers/interface/QGTagger.h"
#include "TMath.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
using namespace std;
using namespace pat;
using namespace edm;
using namespace reco;
namespace SelfVetoPolicyMu{
  enum SelfVetoPolicyMu{
    selfVetoNone=0, selfVetoAll=1, selfVetoFirst=2
  };
}
/////
//   Class declaration
/////
class MuonSelector : public  baseTree{
 public:
  MuonSelector(std::string name, TTree* tree, bool debug, const edm::ParameterSet& cfg, edm::ConsumesCollector && iC);
  ~MuonSelector();
  void Fill(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  void SetBranches();
  void Clear();
  bool isGoodVertex(const reco::Vertex& vtx);
 private:
  MuonSelector(){};
  /////
  //   Config variables
  /////
  edm::EDGetTokenT<reco::VertexCollection> vtx_h_;
  edm::EDGetTokenT<reco::BeamSpot> beamSpot_;
  edm::EDGetTokenT<edm::View<pat::Muon> > muon_h_;
  edm::EDGetTokenT<pat::PackedCandidateCollection> pfToken_;
  edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
  double _Muon_pt_min;
  double _Muon_eta_max;
  int    _vtx_ndof_min;
  int    _vtx_rho_max;
  double _vtx_position_z_max;
  bool   _is_data;
  bool   _reduced;
  /////
  //   BSM 
  /////
  //Variables
  //Kinematics
  vector<double> Muon_pt, Muon_eta, Muon_phi, Muon_energy, Muon_px, Muon_py, Muon_pz, Muon_p, Muon_dB;
  //Charge
  vector<double> Muon_charge;
  //ID
  vector<int>    Muon_soft, Muon_loose, Muon_medium, Muon_tight, Muon_isHighPt, Muon_POGisGood, Muon_pdgId, Muon_pf, Muon_isGlobal, Muon_isTrackerMuon, Muon_tunePBestTrackType;
  //Isolation
  vector<double> Muon_isoR04Charged, Muon_isoR04NeutralHadron, Muon_isoR04Photon, Muon_isoR04PU, Muon_relIsoDeltaBetaR04, Muon_isoR04CharParPt, Muon_isoR03Charged, Muon_isoR03NeutralHadron, Muon_isoR03Photon, Muon_isoR03PU, Muon_relIsoDeltaBetaR03, Muon_isoR03CharParPt, Muon_trackIso, Muon_TrackerIso, Muon_ecalIso, Muon_hcalIso, Muon_caloIso, Muon_isoSum, Muon_pfEcalEnergy;
  //Track related variables 
  vector<double> Muon_chi2, Muon_chi2LocalPosition, Muon_matchedStat, Muon_validHits, Muon_validHitsInner, Muon_TLayers, Muon_ndof, Muon_validFraction, Muon_pixelLayersWithMeasurement, Muon_qualityhighPurity, Muon_trkKink, Muon_segmentCompatibility; 
  //IP
  vector<double> Muon_dz_pv, Muon_dxy_pv, Muon_dz_bs, Muon_dxy_bs, Muon_dzError, Muon_dxyError, Muon_vtx, Muon_vty, Muon_vtz; 
  // MC
  vector<double> Muon_gen_pt, Muon_gen_eta, Muon_gen_phi, Muon_gen_en;
  vector<int> Muon_gen_pdgId, Muon_gen_isPromptFinalState, Muon_gen_isDirectPromptTauDecayProductFinalState;
  //TRIGGER
  vector<int>    Muon_isMatchedToTrigger;
  //Methods
  int MatchingToTrigger(const edm::Event& iEvent, edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects, edm::Handle<edm::TriggerResults> triggerBits, float eta, float phi);
  double get_effarea(double eta);
};
#endif
