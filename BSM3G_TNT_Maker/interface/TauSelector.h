#ifndef __TAU_MU_H_
#define __TAU_MU_H_
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
#include "DataFormats/PatCandidates/interface/Muon.h"
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
//#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "RecoBTag/BTagTools/interface/SignedTransverseImpactParameter.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "TMath.h"
using namespace std;
using namespace pat;
using namespace edm;
using namespace reco;
/////
//   Class declaration
/////
class TauSelector : public baseTree{
 public:
  TauSelector(std::string name, TTree* tree, bool debug, const edm::ParameterSet& cfg, edm::ConsumesCollector && ic);
  ~TauSelector();
  void Fill(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  void SetBranches();
  void Clear();
  bool isGoodVertex(const reco::Vertex& vtxxx);
 private:
  TauSelector(){};
  /////
  //   Config variables
  /////
  edm::EDGetTokenT<reco::VertexCollection> vtx_h_;
  edm::EDGetTokenT<reco::BeamSpot> beamSpot_;
  edm::EDGetTokenT<edm::View<pat::Tau> > taus_;
  edm::EDGetTokenT<pat::PackedCandidateCollection> pfToken_;
  double _Tau_pt_min;
  double _Tau_eta_max;
  int _Tau_vtx_ndof_min;
  int _Tau_vtx_rho_max;
  double _Tau_vtx_position_z_max;
  bool _reduced;
  /////
  //   BSM variables
  /////
  //Kinematic
  vector<double> Tau_pt, Tau_eta, Tau_phi, Tau_energy, Tau_px, Tau_py, Tau_pz, Tau_p;
  vector<double> Tau_leadChargedCandPt, Tau_leadChargedCandEta, Tau_leadChargedCandPhi, Tau_leadChargedCandE;
  vector<double> Tau_leadChargedCandTrack_pt, Tau_leadChargedCandTrack_ptError;
  //Charge
  vector<double> Tau_charge, Tau_leadChargedCandCharge;
  //Decay mode finding
  vector<int> Tau_decayModeFinding, Tau_decayModeFindingOldDMs, Tau_decayModeFindingNewDMs; 
  //Against Muon
  vector<int> Tau_againstMuonLoose3, Tau_againstMuonTight3; 
  //Against Electron
  vector<int> Tau_againstElectronVLooseMVA6, Tau_againstElectronLooseMVA6, Tau_againstElectronMediumMVA6, Tau_againstElectronTightMVA6;
  //Isolation
  vector<int> Tau_byLooseCombinedIsolationDeltaBetaCorr3Hits, Tau_byMediumCombinedIsolationDeltaBetaCorr3Hits, Tau_byTightCombinedIsolationDeltaBetaCorr3Hits, Tau_byVLooseIsolationMVArun2v1DBoldDMwLT, Tau_byLooseIsolationMVArun2v1DBoldDMwLT, Tau_byMediumIsolationMVArun2v1DBoldDMwLT, Tau_byTightIsolationMVArun2v1DBoldDMwLT, Tau_byVTightIsolationMVArun2v1DBoldDMwLT, Tau_byVLooseIsolationMVArun2v1DBnewDMwLT, Tau_byLooseIsolationMVArun2v1DBnewDMwLT, Tau_byMediumIsolationMVArun2v1DBnewDMwLT, Tau_byTightIsolationMVArun2v1DBnewDMwLT, Tau_byVTightIsolationMVArun2v1DBnewDMwLT, Tau_byLooseCombinedIsolationDeltaBetaCorr3HitsdR03, Tau_byMediumCombinedIsolationDeltaBetaCorr3HitsdR03, Tau_byTightCombinedIsolationDeltaBetaCorr3HitsdR03, Tau_byVLooseIsolationMVArun2v1DBdR03oldDMwLT, Tau_byLooseIsolationMVArun2v1DBdR03oldDMwLT, Tau_byMediumIsolationMVArun2v1DBdR03oldDMwLT, Tau_byTightIsolationMVArun2v1DBdR03oldDMwLT, Tau_byVTightIsolationMVArun2v1DBdR03oldDMwLT, Tau_byVVLooseIsolationMVArun2017v2DBoldDMwLT2017, Tau_byVLooseIsolationMVArun2017v2DBoldDMwLT2017, Tau_byLooseIsolationMVArun2017v2DBoldDMwLT2017, Tau_byMediumIsolationMVArun2017v2DBoldDMwLT2017, Tau_byTightIsolationMVArun2017v2DBoldDMwLT2017, Tau_byVTightIsolationMVArun2017v2DBoldDMwLT2017, Tau_byVVTightIsolationMVArun2017v2DBoldDMwLT2017, Tau_byVVLooseIsolationMVArun2017v2DBnewDMwLT2017, Tau_byVLooseIsolationMVArun2017v2DBnewDMwLT2017, Tau_byLooseIsolationMVArun2017v2DBnewDMwLT2017, Tau_byMediumIsolationMVArun2017v2DBnewDMwLT2017, Tau_byTightIsolationMVArun2017v2DBnewDMwLT2017, Tau_byVTightIsolationMVArun2017v2DBnewDMwLT2017, Tau_byVVTightIsolationMVArun2017v2DBnewDMwLT2017, Tau_byVVLooseIsolationMVArun2017v2DBoldDMdR0p3wLT2017, Tau_byVLooseIsolationMVArun2017v2DBoldDMdR0p3wLT2017, Tau_byLooseIsolationMVArun2017v2DBoldDMdR0p3wLT2017, Tau_byMediumIsolationMVArun2017v2DBoldDMdR0p3wLT2017, Tau_byTightIsolationMVArun2017v2DBoldDMdR0p3wLT2017, Tau_byVTightIsolationMVArun2017v2DBoldDMdR0p3wLT2017, Tau_byVVTightIsolationMVArun2017v2DBoldDMdR0p3wLT2017;
  vector<double> Tau_byIsolationMVA3oldDMwLTraw, Tau_chargedIsoPtSum, Tau_neutralIsoPtSum, Tau_puCorrPtSum;
  //Other prop and Track related variables
  vector<double> Tau_nProngs, Tau_leadChargedCandNdof, Tau_leadChargedCandChi2, Tau_leadChargedCandValidHits;
};
#endif
