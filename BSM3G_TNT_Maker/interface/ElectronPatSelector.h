// 
//Authors:  Andres Florez: Universidad de los Andes, Colombia. 
//          kaur amandeepkalsi: Panjab University, India. 
//
#ifndef __ELECTRON_PAT_H_
#define __ELECTRON_PAT_H_
/////
//   Include files
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
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "EGammaMvaEleEstimatorFWLite.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "Math/VectorUtil.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "DataFormats/Math/interface/deltaR.h"
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
class ElectronPatSelector : public  baseTree{
 public:
  ElectronPatSelector(std::string name, TTree* tree, bool debug, const edm::ParameterSet& cfg, edm::ConsumesCollector && iC);
  ~ElectronPatSelector();
  void Fill(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  void SetBranches();
  void Clear();
  bool isGoodVertex(const reco::Vertex& vtx);
 private:
  ElectronPatSelector(){};
  EGammaMvaEleEstimatorFWLite* mvaID_;
  /////
  //   Config variables
  /////
  edm::InputTag _patElectronToken;
  edm::InputTag _vertexInputTag;
  edm::InputTag _beamSpot;
  double _patElectron_pt_min;
  double _patElectron_eta_max;
  int    _vtx_ndof_min;
  int    _vtx_rho_max;
  double _vtx_position_z_max;
  edm::EDGetTokenT<pat::PackedCandidateCollection> pfToken_;
  edm::InputTag jetToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > electronVetoIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > electronLooseIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > electronMediumIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > electronTightIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > eleMVATrigIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > eleHEEPIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<float> > elemvaValuesMapToken_;
  edm::EDGetTokenT<edm::ValueMap<int> >   elemvaCategoriesMapToken_;
  bool _AJVar;
  bool _tthlepVar;
  bool _is_data;
  /////
  //   BSM 
  /////
  //Variables
  //Kinematics
  vector<double> patElectron_pt, patElectron_eta, patElectron_phi, patElectron_energy, patElectron_Et, patElectron_SCeta, patElectron_inCrack;
  //Charge
  vector<double> patElectron_charge;
  //ID
  vector<int>  passVetoId_, passLooseId_, passMediumId_, passTightId_, passMvatrigId_, passHEEPId_, patElectron_mvaCategory_, patElectron_pdgId, patElectron_isEcalDriven;
  vector<float> patElectron_mvaValue_;
  //Isolation
  vector<double> patElectron_isoChargedHadrons, patElectron_isoNeutralHadrons, patElectron_isoPhotons, patElectron_isoPU, patElectron_relIsoDeltaBeta, patElectron_relIsoRhoEA, patElectron_dr03EcalRecHitSumEt, patElectron_dr03HcalDepth1TowerSumEt, patElectron_isolPtTracks, patElectron_ecalPFClusterIso, patElectron_hcalPFClusterIso;
  //Shape, Track related variables, other prop
  vector<double> patElectron_dEtaIn, patElectron_dPhiIn, 
                 patElectron_full5x5_sigmaIetaIeta, patElectron_full5x5_e2x5Max, patElectron_full5x5_e5x5, patElectron_full5x5_e1x5,
                 patElectron_hOverE, patElectron_ooEmooP, passConversionVeto_, expectedMissingInnerHits, patElectron_gsfTrack_ndof, patElectron_gsfTrack_normChi2;
  //IP
  vector<double> patElectron_gsfTrack_dz_pv, patElectron_gsfTrack_dxy_pv, patElectron_d0, patElectron_gsfTrack_dz_bs, patElectron_gsfTrack_dxy_bs, patElectron_dzError, patElectron_dxyError, patElectron_gsfTrack_vtx, patElectron_gsfTrack_vty, patElectron_gsfTrack_vtz;
  vector<double> patElectron_gsfTrack_PCAx_pv, patElectron_gsfTrack_PCAy_pv, patElectron_gsfTrack_PCAz_pv,
                 patElectron_gsfTrack_PCAx_bs, patElectron_gsfTrack_PCAy_bs, patElectron_gsfTrack_PCAz_bs,  
                 patElectron_gsfTrackFitErrorMatrix_00, patElectron_gsfTrackFitErrorMatrix_01, patElectron_gsfTrackFitErrorMatrix_02, patElectron_gsfTrackFitErrorMatrix_11, patElectron_gsfTrackFitErrorMatrix_12, patElectron_gsfTrackFitErrorMatrix_22;
  //patElectron_relIsoDeltaBeta,patElectron_relIsoRhoEA;
  /////
  //   TTH
  /////
  //Methods
  void get_eleminiIso_info(const pat::PackedCandidateCollection& pcc,double rho, const pat::Electron& cand, double& miniIso, double& miniIsoCh, double& miniIsoNeu, double& miniIsoPUsub);
  void get_chneupu_pcc(const pat::PackedCandidateCollection& pcc,vector<const pat::PackedCandidate *>& pfc_all,vector<const pat::PackedCandidate *>& pfc_ch,vector<const pat::PackedCandidate *>& pfc_neu,vector<const pat::PackedCandidate *>&pfc_pu);
  double get_isosumraw(const std::vector<const pat::PackedCandidate *> & pcc, const pat::Electron& cand, double IsoConeSize, double innerR, double ptTh, int pdgId);
  double get_effarea(double eta);
  void get_elejet_info(edm::View<pat::Electron>::const_iterator& ele, const edm::Event& iEvent, const edm::EventSetup& iSetup, double& elejet_mindr, double& elejet_pt, double& eleptSUelejetpt, double& elejet_btagdisc, double& jx, double& jy, double& jz, double& eleptrel);
  //Variables
  vector<double> patElectron_miniIsoRel, patElectron_miniIsoCh, patElectron_miniIsoNeu, patElectron_miniIsoPUsub, patElectron_jetdr, patElectron_jetpt, patElectron_jetptratio, patElectron_jetcsv, patElectron_ptrel, patElectron_IP3Dsig, patElectron_eleMVASpring15NonTrig25ns, patElectron_eleMVASpring15NonTrig25ns_VL;
  /////
  //   MC
  /////
   vector<double> patElectron_gen_pt, patElectron_gen_eta, patElectron_gen_phi, patElectron_gen_en;
   vector<int>    patElectron_gen_pdgId, patElectron_gen_isPromptFinalState, patElectron_gen_isDirectPromptTauDecayProductFinalState;
};
#endif
