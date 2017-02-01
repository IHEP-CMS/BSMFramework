//====================================
//       TTHbb_eventSelector.h
//====================================
//      Joshuha Thomas-Wilsker
//   joshuha.thomas-wilsker@cern.ch
//  Institute Of High Energy Physics
//====================================
// Module performs TTHbb analysis
// selection and fills boolean in tree
// signifiying whether an event passed
// or failed the selection.
//====================================
#ifndef TTHbb_eventSelECTOR_H
#define TTHbb_eventSelECTOR_H

// system include files
#include <memory>
#include <map>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"//Needed for edm::ValueMap<>

// additional include files
#include "baseTree.h"
#include <TBranch.h>
#include <TTree.h>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <ios>

#include <cmath>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TRandom3.h>
#include <TBranch.h>
#include <TClonesArray.h>

#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerEventWithRefs.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"

//#include "DataFormats/Candidate/interface/Candidate.h"
//#include "DataFormats/Math/interface/deltaR.h"


#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "JetMETCorrections/Modules/interface/JetResolution.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"

#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
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

using namespace std;
using namespace edm;
using namespace reco;
using namespace pat;

// Want to inherit functions from baseTree.
class TTHbb_eventSelector : public baseTree{
public:
  // Double ampersand syntax declares an rvalue reference.
  // Unlike single ampersand, can bind to an rvalue like a temporary without having to be const.
  // Moves the otherwise temporary object (ic) into the object being constructed.

  TTHbb_eventSelector(std::string name, TTree* tree, bool debug, const edm::ParameterSet& cfg, edm::ConsumesCollector && iC);
  ~TTHbb_eventSelector();
  void Fill(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  bool isGoodVertex(const reco::Vertex& vtx);
  void JECInitialisation();
  void GetJER(pat::Jet jet, float JesSF, float rhoJER, bool AK4PFchs, float &JERScaleFactor, float &JERScaleFactorUP, float &JERScaleFactorDOWN);
  void SetBranches();
  void FillLeptonVectors(const edm::Event& iEvent, std::vector<const reco::Candidate*>& looseleps, std::vector<const reco::Candidate*>& loose_electrons, std::vector<const reco::Candidate*>& loose_muons, std::vector<const reco::Candidate*>& tightleps, std::vector<const reco::Candidate*>& tight_electrons, std::vector<const reco::Candidate*>& tight_muons);
  void FillJetVectors(const edm::Event& iEvent, std::vector<const reco::Candidate*> looseleps, std::vector<const pat::Jet*>& leading_jets, std::vector<const pat::Jet*>& subleading_jets ,std::vector<const pat::Jet*>& leading_btags, std::vector<const pat::Jet*>& subleading_btags, vector<pair<double,int> >& leading_jet_csv_pos, vector<pair<double,int> >& subleading_jet_csv_pos);
  int dilepton_channel(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  int singlelepton_channel(const edm::Event& iEvent, const edm::EventSetup& iSetup);

  //void Clear();
  //bool ApplySelection();

private:
  //TTHbb_eventSelector(){};//Line may be deprecated.
  edm::EDGetTokenT<edm::ValueMap<bool>  > electronVetoIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool>  > electronLooseIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool>  > electronMediumIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool>  > electronTightIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool>  > eleMVATrigIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool>  > eleMVAnonTrigIdMap_;
  edm::EDGetTokenT<edm::ValueMap<bool>  > eleMVATrigwp90IdMap_;
  edm::EDGetTokenT<edm::ValueMap<bool>  > eleMVAnonTrigwp90IdMap_;
  edm::EDGetTokenT<edm::ValueMap<bool>  > eleHEEPIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<float> > elemvaValuesMapToken_nonTrig_;
  edm::EDGetTokenT<edm::ValueMap<int>   > elemvaCategoriesMapToken_nonTrig_;
  edm::EDGetTokenT<edm::ValueMap<float> > elemvaValuesMapToken_Trig_;
  edm::EDGetTokenT<edm::ValueMap<int>   > elemvaCategoriesMapToken_Trig_;
  edm::EDGetTokenT<edm::View<pat::Electron> > electron_pat_;
  edm::EDGetTokenT<edm::View<pat::Muon> > muon_h_;
  edm::EDGetTokenT<pat::JetCollection> jets_;
  edm::EDGetTokenT<double> rhopogHandle_;
  edm::EDGetTokenT<double> rhoJERHandle_;
  edm::FileInPath jecPayloadNamesAK4PFchsMC1_;
  edm::FileInPath jecPayloadNamesAK4PFchsMC2_;
  edm::FileInPath jecPayloadNamesAK4PFchsMC3_;
  edm::FileInPath jecPayloadNamesAK4PFchsMCUnc_;
  edm::FileInPath jecPayloadNamesAK4PFchsDATA1_;
  edm::FileInPath jecPayloadNamesAK4PFchsDATA2_;
  edm::FileInPath jecPayloadNamesAK4PFchsDATA3_;
  edm::FileInPath jecPayloadNamesAK4PFchsDATA4_;
  edm::FileInPath jecPayloadNamesAK4PFchsDATAUnc_;
  std::string jerAK4PFchs_;
  std::string jerAK4PFchsSF_;
  edm::EDGetTokenT<reco::VertexCollection> vtx_h_;
  int    _vtx_ndof_min;
  int    _vtx_rho_max;
  double _vtx_position_z_max;
  int tmp_evtNum;
  bool _is_data;

  edm::EDGetTokenT<edm::View<reco::GenParticle> > prunedGenToken_;
  boost::shared_ptr<FactorizedJetCorrector>   jecAK4PFchsMC_;
  boost::shared_ptr<JetCorrectionUncertainty> jecAK4PFchsMCUnc_;
  boost::shared_ptr<FactorizedJetCorrector>   jecAK4PFchsDATA_;
  boost::shared_ptr<JetCorrectionUncertainty> jecAK4PFchsDATAUnc_;

  bool is_tight_muon(const pat::Muon& mu,const reco::Vertex& vtx);
  bool is_loose_muon(const pat::Muon& mu,const reco::Vertex& vtx);
  double rel_iso_dbc_mu(const pat::Muon& lepton);
  double rel_iso_dbc_ele(const pat::Electron& lepton, double rhopog);
  double get_effarea(double eta);
  bool is_good_jet(const pat::Jet &j, double rho, double rhoJER, int vtxsize, double channel_jpt_cut);
  bool is_loose_electron(const pat::Electron& ele, double rhopog);//, const reco::Vertex& vtx);
  bool is_tight_electron(const pat::Electron& ele, double rhopog);//, const reco::Vertex& vtx);

  int is_singlelep,is_dilepton;
  int is_SL, is_DL;
  vector<double> evSel_numjet,evSel_jetpt, evSel_jeteta, evSel_jetphi, evSel_jetenergy;

};

#endif
