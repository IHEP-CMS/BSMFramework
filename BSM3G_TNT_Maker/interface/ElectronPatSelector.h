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
namespace SelfVetoPolicyEle{
  enum SelfVetoPolicyEle{
    selfVetoNone=0, selfVetoAll=1, selfVetoFirst=2
  };
}
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
  edm::EDGetTokenT<edm::View<pat::Jet>> jetsToken;
  edm::EDGetTokenT<edm::ValueMap<float> > qgToken;
  edm::EDGetTokenT<edm::ValueMap<bool>  > electronVetoIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool>  > electronLooseIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool>  > electronMediumIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool>  > electronTightIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool>  > eleMVATrigIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool>  > eleHEEPIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<float> > elemvaValuesMapToken_nonTrig_;
  edm::EDGetTokenT<edm::ValueMap<int>   > elemvaCategoriesMapToken_nonTrig_;
  edm::EDGetTokenT<edm::ValueMap<float> > elemvaValuesMapToken_Trig_;
  edm::EDGetTokenT<edm::ValueMap<int>   > elemvaCategoriesMapToken_Trig_;
  bool _AJVar;
  bool _tthlepVar;
  bool _is_data;
  /////
  //   BSM 
  /////
  //Variables
  //Kinematics
  vector<double> patElectron_pt, patElectron_eta, patElectron_phi, patElectron_energy, patElectron_px, patElectron_py, patElectron_pz, patElectron_p, patElectron_Et, patElectron_SCeta, patElectron_inCrack;
  //Charge
  vector<double> patElectron_charge;
  //ID
  vector<int>  passVetoId_, passLooseId_, passMediumId_, passTightId_, passHEEPId_, passMvatrigId_, patElectron_mvaCategory_, patElectron_pdgId, patElectron_isEcalDriven;
  vector<float> patElectron_mvaValue_nonTrig_, patElectron_mvaCategory_nonTrig_, patElectron_mvaValue_Trig_, patElectron_mvaCategory_Trig_;
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
  double get_isosumraw(const std::vector<const pat::PackedCandidate *> & pcc, const pat::Electron& cand, double IsoConeSize, double innerR, double ptTh, SelfVetoPolicyEle::SelfVetoPolicyEle selfVeto, int pdgId);
  double get_effarea(double eta);
  void get_elejet_info(edm::View<pat::Electron>::const_iterator& ele, const edm::Event& iEvent, const edm::EventSetup& iSetup,
                       double& elejet_mindr, double& elejet_pt, double& eleptOVelejetpt,
                       double& elejet_pfCombinedInclusiveSecondaryVertexV2BJetTags, double& elejet_pfJetProbabilityBJetTags, double& elejet_pfCombinedMVABJetTags, double& elejet_qgl,                       
                       double& jx, double& jy, double& jz, double& eleptrel,
                       int& lepjetidx);
  int pvassociation(edm::View<pat::Electron>::const_iterator& ele, const pat::PackedCandidateCollection& pcc);
  double relativeEta(const math::XYZVector& vector, const math::XYZVector& axis);
  double get_lepWmass(edm::View<pat::Electron>::const_iterator& ele, const edm::Event& iEvent, int& lepjetidx);
  double get_lepTopmass(edm::View<pat::Electron>::const_iterator& ele, const edm::Event& iEvent, int& lepjetidx);
  double get_lepWTopmass(edm::View<pat::Electron>::const_iterator& ele, const edm::Event& iEvent, int& lepjetidx);
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
  vector<double> patElectron_miniIsoRel, patElectron_miniIsoCh, patElectron_miniIsoNeu, patElectron_miniIsoPUsub;
  vector<double> patElectron_jetdr, patElectron_jetpt, patElectron_jetptratio, patElectron_jetcsv, patElectron_ptrel, patElectron_IP3Dsig, patElectron_eleMVASpring15NonTrig25ns, patElectron_eleMVASpring15NonTrig25ns_VL;
  vector<double> patElectron_pvass, patElectron_etarel, patElectron_ptOVen, patElectron_elejet_pfJetProbabilityBJetTag, patElectron_elejet_pfCombinedMVABJetTags, patElectron_elejet_qgl;
  vector<double> patElectron_elemass, patElectron_elejet_mass, patElectron_elejet_Wmass, patElectron_elejet_Topmass, patElectron_elejet_WTopmass;
  vector<double> patElectron_IP3D_val, patElectron_IP3D_err, patElectron_IP3D_sig, patElectron_IP2D_val, patElectron_IP2D_err, patElectron_IP2D_sig, patElectron_sIP3D_val, patElectron_sIP3D_err, patElectron_sIP3D_sig, patElectron_sIP2D_val, patElectron_sIP2D_err, patElectron_sIP2D_sig, patElectron_IP1D_val, patElectron_IP1D_err, patElectron_IP1D_sig, patElectron_sIP1D_val, patElectron_sIP1D_err, patElectron_sIP1D_sig;
  vector<double> patElectron_lepjetMaxIP3D_val, patElectron_lepjetMaxIP3D_sig, patElectron_lepjetMaxsIP3D_val, patElectron_lepjetMaxsIP3D_sig, patElectron_lepjetMaxIP2D_val, patElectron_lepjetMaxIP2D_sig, patElectron_lepjetMaxsIP2D_val, patElectron_lepjetMaxsIP2D_sig, patElectron_lepjetMaxIP1D_val, patElectron_lepjetMaxIP1D_sig, patElectron_lepjetMaxsIP1D_val, patElectron_lepjetMaxsIP1D_sig, patElectron_lepjetAvIP3D_val, patElectron_lepjetAvIP3D_sig, patElectron_lepjetAvsIP3D_val, patElectron_lepjetAvsIP3D_sig, patElectron_lepjetAvIP2D_val, patElectron_lepjetAvIP2D_sig, patElectron_lepjetAvsIP2D_val, patElectron_lepjetAvsIP2D_sig, patElectron_lepjetAvIP1D_val, patElectron_lepjetAvIP1D_sig, patElectron_lepjetAvsIP1D_val, patElectron_lepjetAvsIP1D_sig;
  vector<double> patElectron_lepjetchtrks, patElectron_lepjetpvchtrks, patElectron_lepjetnonpvchtrks, patElectron_lepjetndaus;
  vector<double> patElectron_lepjetpvchi2, patElectron_lepjetnumno2trk;
  /////
  //   MC
  /////
  vector<double> patElectron_gen_pt, patElectron_gen_eta, patElectron_gen_phi, patElectron_gen_en;
  vector<int>    patElectron_gen_pdgId, patElectron_gen_isPromptFinalState, patElectron_gen_isDirectPromptTauDecayProductFinalState;
};
#endif
