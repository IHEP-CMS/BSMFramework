#ifndef __TAUBOOSTED_MU_H_
#define __TAUBOOSTED_MU_H_
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
#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"
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
class TauboostedSelector : public baseTree{
 public:
  TauboostedSelector(std::string name, TTree* tree, bool debug, const edm::ParameterSet& cfg, edm::ConsumesCollector && ic);
  ~TauboostedSelector();
  void Fill(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  void SetBranches();
  void Clear();
  bool isGoodVertex(const reco::Vertex& vtxxx);
 private:
  TauboostedSelector(){};
  /////
  //   Config variables
  /////
  edm::EDGetTokenT<reco::VertexCollection> vtx_h_;
  edm::EDGetTokenT<reco::BeamSpot> beamSpot_;
  edm::EDGetTokenT<edm::View<pat::Tau> > tausboosted_;
  edm::EDGetTokenT<pat::PackedCandidateCollection> pfToken_;
  double _Tauboosted_pt_min;
  double _Tauboosted_eta_max;
  int _Tauboosted_vtx_ndof_min;
  int _Tauboosted_vtx_rho_max;
  double _Tauboosted_vtx_position_z_max;
  bool _reduced;
  /////
  //   BSM variables
  /////
  //Kinematic
  vector<double> Tauboosted_pt, Tauboosted_eta, Tauboosted_phi, Tauboosted_energy, Tauboosted_px, Tauboosted_py, Tauboosted_pz, Tauboosted_p;
  vector<double> Tauboosted_leadChargedCandPt, Tauboosted_leadChargedCandEta, Tauboosted_leadChargedCandPhi, Tauboosted_leadChargedCandE;
  vector<double> Tauboosted_leadChargedCandTrack_pt, Tauboosted_leadChargedCandTrack_ptError;
  //Charge
  vector<double> Tauboosted_charge, Tauboosted_leadChargedCandCharge;
  //Decay mode finding
  vector<int> Tauboosted_decayModeFinding, Tauboosted_decayModeFindingOldDMs, Tauboosted_decayModeFindingNewDMs; 
  //Against Muon
  vector<int> Tauboosted_againstMuonLoose3, Tauboosted_againstMuonTight3; 
  //Against Electron
  vector<int> Tauboosted_againstElectronVLooseMVA6, Tauboosted_againstElectronLooseMVA6, Tauboosted_againstElectronMediumMVA6, Tauboosted_againstElectronTightMVA6;
  //Isolation
  vector<int> Tauboosted_byLooseCombinedIsolationDeltaBetaCorr3Hits, Tauboosted_byMediumCombinedIsolationDeltaBetaCorr3Hits, Tauboosted_byTightCombinedIsolationDeltaBetaCorr3Hits, Tauboosted_byVLooseIsolationMVArun2v1DBoldDMwLT, Tauboosted_byLooseIsolationMVArun2v1DBoldDMwLT, Tauboosted_byMediumIsolationMVArun2v1DBoldDMwLT, Tauboosted_byTightIsolationMVArun2v1DBoldDMwLT, Tauboosted_byVTightIsolationMVArun2v1DBoldDMwLT, Tauboosted_byVLooseIsolationMVArun2v1DBnewDMwLT, Tauboosted_byLooseIsolationMVArun2v1DBnewDMwLT, Tauboosted_byMediumIsolationMVArun2v1DBnewDMwLT, Tauboosted_byTightIsolationMVArun2v1DBnewDMwLT, Tauboosted_byVTightIsolationMVArun2v1DBnewDMwLT, Tauboosted_byLooseCombinedIsolationDeltaBetaCorr3HitsdR03, Tauboosted_byMediumCombinedIsolationDeltaBetaCorr3HitsdR03, Tauboosted_byTightCombinedIsolationDeltaBetaCorr3HitsdR03, Tauboosted_byVLooseIsolationMVArun2v1DBdR03oldDMwLT, Tauboosted_byLooseIsolationMVArun2v1DBdR03oldDMwLT, Tauboosted_byMediumIsolationMVArun2v1DBdR03oldDMwLT, Tauboosted_byTightIsolationMVArun2v1DBdR03oldDMwLT, Tauboosted_byVTightIsolationMVArun2v1DBdR03oldDMwLT;
  vector<double> Tauboosted_byIsolationMVA3oldDMwLTraw, Tauboosted_chargedIsoPtSum, Tauboosted_neutralIsoPtSum, Tauboosted_puCorrPtSum;
  //Other prop and Track related variables
  vector<double> Tauboosted_nProngs, Tauboosted_leadChargedCandNdof, Tauboosted_leadChargedCandChi2, Tauboosted_leadChargedCandValidHits;
};
#endif
