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
  edm::EDGetTokenT<pat::JetCollection> jets_;
  edm::EDGetTokenT<edm::View<pat::Jet>> jetsToken;
  edm::EDGetTokenT<pat::PackedCandidateCollection> pfToken_;
  edm::EDGetTokenT<double> rhoHandle_;
  edm::EDGetTokenT<edm::ValueMap<float> > qgToken_;  
  double _Muon_pt_min;
  double _Muon_eta_max;
  int    _vtx_ndof_min;
  int    _vtx_rho_max;
  double _vtx_position_z_max;
  bool   _super_TNT; //super tiny ntuple?
  bool   _is_data;
  bool _AJVar;
  bool _tthlepVar;
  bool _qglVar;
  /////
  //   BSM 
  /////
  //Variables
  //Kinematics
  vector<double> Muon_pt, Muon_eta, Muon_phi, Muon_energy, Muon_px, Muon_py, Muon_pz, Muon_p, Muon_dB, Muon_pt_it, Muon_ptErr_it, Muon_pt_bt, Muon_pTErrOVpT_it, Muon_ptErr_bt, Muon_pTErrOVpT_bt, Muon_pt_tunePbt;
  //Charge
  vector<double> Muon_charge;
  //ID
  vector<int>    Muon_soft, Muon_loose, Muon_medium, Muon_tight, Muon_isHighPt, Muon_POGisGood, Muon_pdgId, Muon_pf, Muon_isGlobal, Muon_isTrackerMuon, Muon_tunePBestTrackType;
  //Isolation
  vector<double> Muon_isoR04Charged, Muon_isoR04NeutralHadron, Muon_isoR04Photon, Muon_isoR04PU, Muon_relIsoDeltaBetaR04, Muon_isoR04CharParPt, 
                 Muon_isoR03Charged, Muon_isoR03NeutralHadron, Muon_isoR03Photon, Muon_isoR03PU, Muon_relIsoDeltaBetaR03, Muon_isoR03CharParPt, 
                 Muon_trackIso, Muon_TrackerIso, Muon_ecalIso, Muon_hcalIso, Muon_caloIso, Muon_isoSum, Muon_pfEcalEnergy;
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
  double get_isosumraw(const std::vector<const pat::PackedCandidate *> & pcc, const pat::Muon& cand, double IsoConeSize, double innerR, double ptTh, SelfVetoPolicyMu::SelfVetoPolicyMu selfVeto, int pdgId);
  double get_effarea(double eta);
  void get_mujet_info(const pat::Muon& mu, const edm::Event& iEvent, const edm::EventSetup& iSetup,
                      double& mujet_mindr, double& mujet_pt, double& muptOVmujetpt,
                      double& mujet_pfCombinedInclusiveSecondaryVertexV2BJetTags, double& mujet_pfJetProbabilityBJetTags, double& mujet_pfCombinedMVABJetTags, double& mujet_qgl,
                      double& jx, double& jy, double& jz, double& muptrel, double& muparel, double& muparelsub,
                      int& lepjetidx);
  int pvassociation(const pat::Muon& mu, const pat::PackedCandidateCollection& pcc);
  double relEta(const math::XYZVector& dir, const math::XYZVector& track);
  double parEta(const math::XYZVector& dir, const math::XYZVector& track);
  double get_lepWmass(const pat::Muon& mu, const edm::Event& iEvent, int& lepjetidx);
  double get_lepTopmass(const pat::Muon& mu, const edm::Event& iEvent, int& lepjetidx);
  double get_lepWTopmass(const pat::Muon& mu, const edm::Event& iEvent, int& lepjetidx);
  void IP3D2D(TransientTrack ttrk, const reco::Vertex& vtx, GlobalVector gv, double& IP3D_val,double& IP3D_err,double& IP3D_sig, double& sIP3D_val,double& sIP3D_err,double& sIP3D_sig, double& IP2D_val,double& IP2D_err,double& IP2D_sig, double& sIP2D_val,double& sIP2D_err,double& sIP2D_sig);
  void zIP1D(TransientTrack ttrk, const reco::Vertex& vtx, GlobalVector gv, double& IP1D_val,double& IP1D_err,double& IP1D_sig, double& sIP1D_val,double& sIP1D_err,double& sIP1D_sig);
  void lepjetIP(const pat::Jet& jet, const reco::Vertex& vtx, GlobalVector lepjetgv, const TransientTrackBuilder& ttrkbuilder,
                double& lepjetMaxIP3D_val, double& lepjetMaxIP3D_sig, double& lepjetMaxsIP3D_val, double& lepjetMaxsIP3D_sig, double& lepjetMaxIP2D_val, double& lepjetMaxIP2D_sig, double& lepjetMaxsIP2D_val, double& lepjetMaxsIP2D_sig, double& lepjetMaxIP1D_val, double& lepjetMaxIP1D_sig, double& lepjetMaxsIP1D_val, double& lepjetMaxsIP1D_sig,
                double& lepjetAvIP3D_val, double& lepjetAvIP3D_sig, double& lepjetAvsIP3D_val, double& lepjetAvsIP3D_sig, double& lepjetAvIP2D_val, double& lepjetAvIP2D_sig, double& lepjetAvsIP2D_val, double& lepjetAvsIP2D_sig, double& lepjetAvIP1D_val, double& lepjetAvIP1D_sig, double& lepjetAvsIP1D_val, double& lepjetAvsIP1D_sig,
                double& denlepjetAvIP3D_val, double& denlepjetAvIP3D_sig, double& denlepjetAvsIP3D_val, double& denlepjetAvsIP3D_sig, double& denlepjetAvIP2D_val, double& denlepjetAvIP2D_sig, double& denlepjetAvsIP2D_val, double& denlepjetAvsIP2D_sig, double& denlepjetAvIP1D_val, double& denlepjetAvIP1D_sig, double& denlepjetAvsIP1D_val, double& denlepjetAvsIP1D_sig,
                double& Lep_IP3D_val
               );
  void lepjetTrks(const pat::Jet& jet, const reco::Vertex& vtx, double& lepjetchtrks, double& lepjetpvchtrks, double& lepjetnonpvchtrks, double& lepjetndaus);
  void lepjetVtxCompatibility(const pat::Jet& jet, const reco::Vertex& vtx, const TransientTrackBuilder& ttrkbuilder, double& lepjetpvchi2, double& lepjetnumno2tr);
  void get_2trksinfo(vector<TransientTrack> ttrks, double& num2v, double& numno2v);
  bool is_goodtrk(Track trk,const reco::Vertex& vtx);
  //Variables
  vector<double> Muon_miniIsoRel, Muon_miniIsoCh, Muon_miniIsoNeu, Muon_miniIsoPUsub;
  vector<double> Muon_jetdr, Muon_jetpt, Muon_jetptratio, Muon_jetcsv, Muon_ptrel, Muon_parel, Muon_parelsub, Muon_etarel, Muon_etarelsub, Muon_etapar, Muon_etaparsub, Muon_IP3Dsig_it;
  vector<double> Muon_pvass, Muon_ptOVen, Muon_mujet_pfJetProbabilityBJetTag, Muon_mujet_pfCombinedMVABJetTags, Muon_mujet_qgl;
  vector<double> Muon_mumass, Muon_mujet_mass, Muon_mujet_Wmass, Muon_mujet_Topmass, Muon_mujet_WTopmass;
  vector<double> Muon_IP3D_val, Muon_IP3D_err, Muon_IP3D_sig, Muon_IP2D_val, Muon_IP2D_err, Muon_IP2D_sig, Muon_sIP3D_val, Muon_sIP3D_err, Muon_sIP3D_sig, Muon_sIP2D_val, Muon_sIP2D_err, Muon_sIP2D_sig, Muon_IP1D_val, Muon_IP1D_err, Muon_IP1D_sig, Muon_sIP1D_val, Muon_sIP1D_err, Muon_sIP1D_sig;
  vector<double> Muon_lepjetMaxIP3D_val, Muon_lepjetMaxIP3D_sig, Muon_lepjetMaxsIP3D_val, Muon_lepjetMaxsIP3D_sig, Muon_lepjetMaxIP2D_val, Muon_lepjetMaxIP2D_sig, Muon_lepjetMaxsIP2D_val, Muon_lepjetMaxsIP2D_sig, Muon_lepjetMaxIP1D_val, Muon_lepjetMaxIP1D_sig, Muon_lepjetMaxsIP1D_val, Muon_lepjetMaxsIP1D_sig, Muon_lepjetAvIP3D_val, Muon_lepjetAvIP3D_sig, Muon_lepjetAvsIP3D_val, Muon_lepjetAvsIP3D_sig, Muon_lepjetAvIP2D_val, Muon_lepjetAvIP2D_sig, Muon_lepjetAvsIP2D_val, Muon_lepjetAvsIP2D_sig, Muon_lepjetAvIP1D_val, Muon_lepjetAvIP1D_sig, Muon_lepjetAvsIP1D_val, Muon_lepjetAvsIP1D_sig;
  vector<double> Muon_lepjetchtrks, Muon_lepjetpvchtrks, Muon_lepjetnonpvchtrks, Muon_lepjetndaus; 
  vector<double> Muon_lepjetpvchi2, Muon_lepjetnumno2trk;
  /////
  //   MC
  /////
  vector<double> Muon_gen_pt, Muon_gen_eta, Muon_gen_phi, Muon_gen_en;
  vector<int> Muon_gen_pdgId, Muon_gen_isPromptFinalState, Muon_gen_isDirectPromptTauDecayProductFinalState;
};
#endif
