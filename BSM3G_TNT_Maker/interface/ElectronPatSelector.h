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
//#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "RecoBTag/BTagTools/interface/SignedTransverseImpactParameter.h"
#include "TMath.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
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
  edm::EDGetTokenT<reco::VertexCollection> vtx_h_;
  edm::EDGetTokenT<reco::BeamSpot> beamSpot_;
  edm::EDGetTokenT<edm::View<pat::Electron> > electron_pat_;
  edm::EDGetTokenT<pat::PackedCandidateCollection> pfToken_;
  edm::EDGetTokenT<double> rhopogHandle_;
  //edm::EDGetTokenT<double> rhotthHandle_;
  edm::EDGetTokenT<edm::ValueMap<bool>  > electronVetoIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool>  > electronLooseIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool>  > electronMediumIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool>  > electronTightIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool>  > eleMVATrigIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool>  > eleMVAnonTrigIdMap_;
  edm::EDGetTokenT<edm::ValueMap<bool>  > eleMVATrigwp90IdMap_;
  edm::EDGetTokenT<edm::ValueMap<bool>  > eleMVAnonTrigwp90IdMap_;
  edm::EDGetTokenT<edm::ValueMap<bool>  > eleMVATrigwpLooseIdMap_;
  edm::EDGetTokenT<edm::ValueMap<bool>  > eleMVAnonTrigwpLooseIdMap_;
  edm::EDGetTokenT<edm::ValueMap<bool>  > eleHEEPIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool>  > eleMVAHZZwpLooseIdMap_;
  edm::EDGetTokenT<edm::ValueMap<float> > elemvaValuesMapToken_nonTrig_;
  edm::EDGetTokenT<edm::ValueMap<int>   > elemvaCategoriesMapToken_nonTrig_;
  edm::EDGetTokenT<edm::ValueMap<float> > elemvaValuesMapToken_Trig_;
  edm::EDGetTokenT<edm::ValueMap<int>   > elemvaCategoriesMapToken_Trig_;
  edm::EDGetTokenT<edm::ValueMap<float> > elemvaValuesMapToken_HZZ_;
  edm::EDGetTokenT<edm::ValueMap<int>   > elemvaCategoriesMapToken_HZZ_;
  edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
  edm::EDGetTokenT<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>>> ebRecHitsToken_;
  edm::Handle<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>>> _ebRecHits;
  double _patElectron_pt_min;
  double _patElectron_eta_max;
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
  vector<double> patElectron_pt, patElectron_eta, patElectron_phi, patElectron_energy, patElectron_energyCorr, patElectron_energySF, patElectron_px, patElectron_py, patElectron_pz, patElectron_p, patElectron_Et, patElectron_SCeta, patElectron_inCrack;
  //Charge
  vector<double> patElectron_charge;
  //ID
  vector<int>  passVetoId_, passLooseId_, passMediumId_, passTightId_, passHEEPId_, passMvatrigId_;
  //Isolation
  vector<double> patElectron_isoChargedHadrons, patElectron_isoNeutralHadrons, patElectron_isoPhotons, patElectron_isoPU, patElectron_relIsoDeltaBeta, patElectron_relIsoRhoEA, patElectron_dr03EcalRecHitSumEt, patElectron_dr03HcalDepth1TowerSumEt, patElectron_isolPtTracks, patElectron_ecalPFClusterIso, patElectron_hcalPFClusterIso;
  //Shape, Track related variables, other prop
  vector<double> patElectron_dEtaIn, patElectron_dPhiIn, patElectron_full5x5_sigmaIetaIeta, patElectron_full5x5_e2x5Max, patElectron_full5x5_e5x5, patElectron_full5x5_e1x5, patElectron_hOverE, patElectron_ooEmooP, passConversionVeto_, expectedMissingInnerHits, patElectron_gsfTrack_ndof, patElectron_gsfTrack_normChi2; 
  //IP
  vector<double> patElectron_gsfTrack_dz_pv, patElectron_gsfTrack_dxy_pv, patElectron_d0, patElectron_dzError, patElectron_dxyError, patElectron_gsfTrack_vtx, patElectron_gsfTrack_vty, patElectron_gsfTrack_vtz;
  //MC 
  vector<double> patElectron_gen_pt, patElectron_gen_eta, patElectron_gen_phi, patElectron_gen_en;
  vector<int>    patElectron_gen_pdgId, patElectron_gen_isPromptFinalState, patElectron_gen_isDirectPromptTauDecayProductFinalState;
  //TRIGGER
  vector<int>    patElectron_isMatchedToTrigger;
  //Methods
  int MatchingToTrigger(const edm::Event& iEvent, edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects, edm::Handle<edm::TriggerResults> triggerBits, float eta, float phi);

  double get_effarea(double eta);
};
#endif
