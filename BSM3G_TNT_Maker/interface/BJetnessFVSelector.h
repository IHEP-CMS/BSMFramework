#ifndef __BJetnessFV_MU_H_
#define __BJetnessFV_MU_H_
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
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
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
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "JetMETCorrections/Modules/interface/JetResolution.h"
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
#include "BSMFramework/BSM3G_TNT_Maker/interface/TTHbb_eventSelector.h"

using namespace std;
using namespace pat;
using namespace edm;
using namespace reco;
/////
//   Class declaration
/////
class BJetnessFVSelector : public  baseTree{
  public:
    BJetnessFVSelector(std::string name, TTree* tree, bool debug, const edm::ParameterSet& cfg, edm::ConsumesCollector && iC);
    ~BJetnessFVSelector();
    void Fill(const edm::Event& iEvent, const edm::EventSetup& iSetup);
    void SetBranches();
    void Clear();
    bool isGoodVertex(const reco::Vertex& vtx);
    void JECInitialization();
    void GetJER(pat::Jet jet, float JesSF, float rhoJER, bool AK4PFchs, float &JERScaleFactor, float &JERScaleFactorUP, float &JERScaleFactorDOWN);
  private:
    BJetnessFVSelector(){};

    TTHbb_eventSelector          *evSel;


    /////
    //   Config variables
    /////
    bool _is_data;
    edm::EDGetTokenT<reco::VertexCollection> vtx_h_;
    edm::EDGetTokenT<edm::View<pat::Electron> > electron_pat_;
    edm::EDGetTokenT<edm::View<pat::Muon> > muon_h_;
    edm::EDGetTokenT<pat::JetCollection> jets_;
    edm::EDGetTokenT<double> rhopogHandle_;
    edm::EDGetTokenT<double> rhoJERHandle_;
    edm::EDGetTokenT<edm::ValueMap<bool>  > eleMVATrigIdMapToken_;
    edm::EDGetTokenT<edm::ValueMap<bool>  > eleMVAnonTrigIdMapToken_;
    int    _vtx_ndof_min;
    int    _vtx_rho_max;
    double _vtx_position_z_max;
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
    boost::shared_ptr<FactorizedJetCorrector>   jecAK4PFchsMC_;
    boost::shared_ptr<JetCorrectionUncertainty> jecAK4PFchsMCUnc_;
    boost::shared_ptr<FactorizedJetCorrector>   jecAK4PFchsDATA_;
    boost::shared_ptr<JetCorrectionUncertainty> jecAK4PFchsDATAUnc_;
    /////
    //   Variables
    /////
    //Evt selection
    int BJetnessFV_isSingleLepton, BJetnessFV_isDoubleLepton;
    //BTag discriminators
    vector<double> BJetnessFV_jetcsv, BJetnessFV_pfJetProbabilityBJetTags, BJetnessFV_pfCombinedMVAV2BJetTags;
    //Num_of_trks
    vector<double> BJetnessFV_num_leps;
    vector<double> BJetnessFV_npvTrkOVcollTrk;
    //ImpactParameter
    vector<double> BJetnessFV_avip3d_val, BJetnessFV_avip3d_sig, BJetnessFV_avsip3d_sig, BJetnessFV_avip1d_sig;
    /////
    //   Methods
    /////
    //Methods to be aligned to the TTHbb selection
    bool is_loose_muon(const pat::Muon& mu, const reco::Vertex& vtx);
    bool is_tight_muon(const pat::Muon& mu, const reco::Vertex& vtx);
    double rel_iso_dbc_mu(const pat::Muon& lepton);
    bool is_loose_electron(const pat::Electron& ele, double rhopog);//, const reco::Vertex& vtx);
    bool is_tight_electron(const pat::Electron& ele, double rhopog);//, const reco::Vertex& vtx);
    double rel_iso_dbc_ele(const pat::Electron& lepton, double rhopog);
    double get_effarea(double eta);
    bool is_good_jet(const pat::Jet &j, double rho, double rhoJER, int vtxsize);
    //Methods for the BJetness variables
    void get_bjetness_vars(
                           vector<pat::Jet> evtjets, const reco::Vertex& vtx, const TransientTrackBuilder& ttrkbuilder, edm::Handle<edm::View<pat::Electron> > electron_pat, edm::Handle<edm::View<pat::Muon> > muon_h,
                           double& bjetnessFV_num_leps, double& bjetnessFV_npvTrkOVcollTrk, double& bjetnessFV_avip3d_val, double& bjetnessFV_avip3d_sig, double& bjetnessFV_avsip3d_sig, double& bjetnessFV_avip1d_sig
                          );
    void get_bjetness_trkinfos(vector<pat::Jet> evtjets, const reco::Vertex& vtx, vector<Track>& jetchtrks, double& bjetness_num_pvtrks, double& bjetness_num_npvtrks, edm::Handle<edm::View<pat::Electron> > electron_pat, edm::Handle<edm::View<pat::Muon> > muon_h, double& bjetness_num_eles, double& bjetness_num_mus, vector<tuple<double, double, double> >& jetsdir);
    bool is_goodtrk(Track trk,const reco::Vertex& vtx);
    bool is_loosePOG_jetmuon(const pat::PackedCandidate &jcand, edm::Handle<edm::View<pat::Muon> > muon_h);
    bool is_softLep_jetelectron(const pat::PackedCandidate &jcand, edm::Handle<edm::View<pat::Electron> > electron_pat, const reco::Vertex& vtx);
    bool is_loosePOGNoIPNoIso_jetelectron(const pat::PackedCandidate &jcand, edm::Handle<edm::View<pat::Electron> > electron_pat, const reco::Vertex& vtx);
    void get_avip3d(vector<Track> trks, const TransientTrackBuilder& ttrkbuilder, reco::Vertex vtx, vector<tuple<double, double, double> >& jetsdir, double& jetchtrks_avip3d_val, double& jetchtrks_avip3d_sig, double& jetchtrks_avsip3d_sig);
    void get_avip1d(vector<Track> trks, const TransientTrackBuilder& ttrkbuilder, reco::Vertex vtx, vector<tuple<double, double, double> >& jetsdir, double& jetchtrks_avip1d_sig);
    vector<TransientTrack> get_ttrks(vector<Track> trks, const TransientTrackBuilder& ttrkbuilder);
};
#endif
