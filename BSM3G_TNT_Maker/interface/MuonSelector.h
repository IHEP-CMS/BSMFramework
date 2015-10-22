// 
//Authors: Andres Florez:      Universidad de los Andes, Colombia. 
//         kaur amandeepkalsi: Panjab University, India. 
// 
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
#include "TMath.h"
using namespace std;
using namespace pat;
using namespace edm;
using namespace reco;
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
  edm::InputTag _muonToken;
  edm::InputTag _vertexInputTag;
  edm::InputTag _beamSpot; 
  double _Muon_pt_min;
  double _Muon_eta_max;
  int    _vtx_ndof_min;
  int    _vtx_rho_max;
  double _vtx_position_z_max;
  bool   _super_TNT; //super tiny ntuple?
  edm::EDGetTokenT<pat::PackedCandidateCollection> pfToken_;
  edm::InputTag jetToken_;
  /////
  //   BSM 
  /////
  //Variables
  //Kinematics
  vector<double> Muon_pt,  Muon_eta, Muon_phi, Muon_energy, Muon_p, Muon_dB, Muon_pt_it, Muon_ptErr_it, Muon_pt_bt, Muon_pTErrOVpT_it, Muon_ptErr_bt, Muon_pTErrOVpT_bt, Muon_pt_tunePbt;
  //Charge
  vector<double> Muon_charge;
  //ID
  vector<int>    Muon_soft, Muon_loose, Muon_medium, Muon_tight, Muon_isHighPt, Muon_POGisGood, Muon_pdgId, Muon_pf, Muon_isGlobal, Muon_isTrackerMuon, Muon_tunePBestTrackType;
  //Isolation
  vector<double> Muon_isoR04Charged, Muon_isoR04NeutralHadron, Muon_isoR04Photon, Muon_isoR04PU, Muon_relIsoDeltaBetaR04, Muon_isoR04CharParPt, 
                 Muon_isoR03Charged, Muon_isoR03NeutralHadron, Muon_isoR03Photon, Muon_isoR03PU, Muon_relIsoDeltaBetaR03, Muon_isoR03CharParPt, 
                 Muon_trackIso, Muon_ecalIso, Muon_hcalIso, Muon_isoSum, Muon_pfEcalEnergy;
  //Track related variables 
  vector<double> Muon_chi2, Muon_chi2LocalPosition, Muon_matchedStat, Muon_validHits, Muon_validHitsInner, Muon_TLayers, Muon_ndof, Muon_validFraction, Muon_pixelLayersWithMeasurement, Muon_qualityhighPurity, Muon_trkKink, Muon_segmentCompatibility; 
  //IP
  vector<double> Muon_dz_pv, Muon_dxy_pv, Muon_dz_bt, Muon_dxy_bt, Muon_dz_bs, Muon_dxy_bs, Muon_dzError, Muon_dxyError, Muon_vtx, Muon_vty, Muon_vtz; 
  vector<double> Muon_track_PCAx_bs, Muon_track_PCAy_bs, Muon_track_PCAz_bs, Muon_track_PCAx_pv, Muon_track_PCAy_pv, Muon_track_PCAz_pv, Muon_trackFitErrorMatrix_00, Muon_trackFitErrorMatrix_01, Muon_trackFitErrorMatrix_02, Muon_trackFitErrorMatrix_11, Muon_trackFitErrorMatrix_12, Muon_trackFitErrorMatrix_22;
  /////
  //   TTH
  /////
  //Methods
  void get_muminiIso_info(const pat::PackedCandidateCollection& pcc,double rho, const pat::Muon& mu, double& miniIso, double& miniIsoCh, double& miniIsoNeu, double& miniIsoPUsub);
  void get_chneupu_pcc(const pat::PackedCandidateCollection& pcc,vector<const pat::PackedCandidate *>& pfc_all,vector<const pat::PackedCandidate *>& pfc_ch,vector<const pat::PackedCandidate *>& pfc_neu,vector<const pat::PackedCandidate *>&pfc_pu);
  double get_isosumraw(const std::vector<const pat::PackedCandidate *> & pcc, const pat::Muon& cand, double IsoConeSize, double innerR, double ptTh, int pdgId);
  double get_effarea(double eta);
  //double get_iso_rho(const pat::Muon& mu, double& rho);
  void get_mujet_info(const pat::Muon& mu, const edm::Event& iEvent, const edm::EventSetup& iSetup, double& mujet_mindr, double& mujet_pt, double& muptOVmujetpt, double& mujet_btagdisc, double& jx, double& jy, double& jz, double& muptrel);
  //Variables
  vector<double> Muon_miniIsoRel, Muon_miniIsoCh, Muon_miniIsoNeu, Muon_miniIsoPUsub, Muon_jetdr, Muon_jetpt, Muon_jetptratio, Muon_jetcsv, Muon_ptrel, Muon_IP3Dsig_it;
};
#endif
