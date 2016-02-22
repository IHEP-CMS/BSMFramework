#ifndef __TauJetness_MU_H_
#define __TauJetness_MU_H_
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
using namespace std;
using namespace pat;
using namespace edm;
using namespace reco;
/////
//   Class declaration
/////
class TauJetnessSelector : public  baseTree{
  public:
    TauJetnessSelector(std::string name, TTree* tree, bool debug, const edm::ParameterSet& cfg, edm::ConsumesCollector && iC);
    ~TauJetnessSelector();
    void Fill(const edm::Event& iEvent, const edm::EventSetup& iSetup);
    void SetBranches();
    void Clear();
  private:
    TauJetnessSelector(){};
    /////
    //   Config variables
    /////
    bool _is_data;
    double _Tau_pt_min;
    double _Tau_eta_max;
    edm::EDGetTokenT<reco::VertexCollection> vtx_h_;
    edm::EDGetTokenT<edm::View<pat::Electron> > electron_pat_;
    edm::EDGetTokenT<edm::View<pat::Muon> > muon_h_;
    edm::EDGetTokenT<edm::View<pat::Tau> > taus_;
    edm::EDGetTokenT<pat::JetCollection> jets_;
    edm::EDGetTokenT<double> rhopogHandle_;
    edm::EDGetTokenT<edm::View<reco::GenParticle> > prunedGenToken_;
    edm::EDGetTokenT<edm::ValueMap<bool> > electronLooseIdMapToken_;
    edm::EDGetTokenT<edm::ValueMap<bool>  > eleMVATrigIdMapToken_;
    edm::EDGetTokenT<edm::ValueMap<float> > elemvaValuesMapToken_Trig_;
    edm::EDGetTokenT<edm::ValueMap<int>   > elemvaCategoriesMapToken_Trig_;
    /////
    //   TTH variables
    /////
    //Gen info
    int TauJetness_ngentau;
    //Kinematics and csv 
    vector<double> TauJetness_numtau;
    vector<double> TauJetness_taupt, TauJetness_taueta, TauJetness_tauphi, TauJetness_tauenergy, TauJetness_charge;
    vector<double> TauJetness_againstMuonTight3, TauJetness_againstElectronMediumMVA6, TauJetness_againstElectronTightMVA6, TauJetness_byTightIsolationMVArun2v1DBnewDMwLT, TauJetness_byVTightIsolationMVArun2v1DBnewDMwLT;
    //Num_of_trks
    vector<double> TauJetness_num_pdgid_leps, TauJetness_num_pdgid_eles, TauJetness_num_pdgid_mus, TauJetness_num_soft_leps, TauJetness_num_soft_eles, TauJetness_num_vetonoipnoiso_leps, TauJetness_num_vetonoipnoiso_eles, TauJetness_num_loosenoipnoiso_leps, TauJetness_num_loosenoipnoiso_eles, TauJetness_num_loose_mus;
    vector<double> TauJetness_numtautrks, TauJetness_numtautrkspv, TauJetness_numtautrksnopv;
    vector<double> TauJetness_npvTrkOVcollTrk, TauJetness_pvTrkOVcollTrk, TauJetness_npvTrkOVpvTrk;
    vector<double> TauJetness_npvPtOVcollPt, TauJetness_pvPtOVcollPt, TauJetness_npvPtOVpvPt;
    //Two_trk_info
    vector<double> TauJetness_avnum2v, TauJetness_avnumno2v, TauJetness_avdca3d2t, TauJetness_avdca3dno2t, TauJetness_avdca3d, TauJetness_avdca2d2t, TauJetness_avdca2dno2t, TauJetness_avdca2d;
    //chi2
    vector<double> TauJetness_chi2;
    //ImpactParameter
    vector<double> TauJetness_avip3d_val, TauJetness_avip3d_sig, TauJetness_avsip3d_val, TauJetness_avsip3d_sig, TauJetness_avip2d_val, TauJetness_avip2d_sig, TauJetness_avsip2d_val, TauJetness_avsip2d_sig, TauJetness_avip1d_val, TauJetness_avip1d_sig, TauJetness_avsip1d_val, TauJetness_avsip1d_sig;
    /////
    //   TTH methods 
    /////
    //bool is_soft_muon(const pat::PackedCandidate &jcand, edm::View<pat::Muon>& );
    bool is_loose_muon(const pat::Muon& mu, const reco::Vertex& vtx);
    double rel_iso_dbc_mu(const pat::Muon& lepton);  
    double rel_iso_dbc_ele(const pat::Electron& lepton, double rhopog);
    double get_effarea(double eta);
    bool is_good_tau(const pat::Jet &j);
    bool is_loosePOG_taumuon(const pat::PackedCandidate &jcand, edm::Handle<edm::View<pat::Muon> > muon_h);
    bool is_softLep_tauelectron(const pat::PackedCandidate &jcand, edm::Handle<edm::View<pat::Electron> > electron_pat, const reco::Vertex& vtx);
    bool is_vetoPOGNoIPNoIso_tauelectron(const pat::PackedCandidate &jcand, edm::Handle<edm::View<pat::Electron> > electron_pat, const reco::Vertex& vtx);
    bool is_loosePOGNoIPNoIso_tauelectron(const pat::PackedCandidate &jcand, edm::Handle<edm::View<pat::Electron> > electron_pat, const reco::Vertex& vtx);
    bool is_loose_electron(const pat::Electron& ele, double rhopog);//, const reco::Vertex& vtx);
    bool is_goodtrk(Track trk,const reco::Vertex& vtx);
    TransientTrack get_ttrk(Track trk, const TransientTrackBuilder& ttrkbuilder);
    vector<TransientTrack> get_ttrks(vector<Track> trks, const TransientTrackBuilder& ttrkbuilder);
    TransientVertex get_tv(vector<Track> trks, const TransientTrackBuilder& ttrkbuilder);
    void get_tautrks(const pat::Jet& tau, const reco::Vertex& vtx, const TransientTrackBuilder& ttrkbuilder, vector<Track>& tauchtrks, vector<Track>& tauchtrkspv, vector<Track>& tauchtrksnpv, vector<tuple<double, double, double> >& tausdir, edm::Handle<edm::View<pat::Electron> > electron_pat, edm::Handle<edm::View<pat::Muon> > muon_h, double& btauness_num_pdgid_eles, double& btauness_num_pdgid_mus, double& btauness_num_soft_eles, double& btauness_num_vetonoipnoiso_eles, double& btauness_num_loosenoipnoiso_eles, double& btauness_num_loose_mus);
    void get_chi2(vector<Track> trks, const TransientTrackBuilder& ttrkbuilder, double& chi2red);
    void get_2trksinfo(vector<Track> trks, const TransientTrackBuilder& ttrkbuilder, double& num2v, double& numno2v, double& dca3d2t, double& dca3dno2t, double& dca2d2t, double& dca2dno2t);
    pair<double,double> dca2trks(Track tkA, Track tkB, const TransientTrackBuilder& ttrkbuilder);
    void get_avip3d(vector<Track> trks, const TransientTrackBuilder& ttrkbuilder, reco::Vertex vtx, vector<tuple<double, double, double> >& tausdir, double& tauchtrks_avip3d_val, double& tauchtrks_avip3d_sig, double& tauchtrks_avsip3d_val, double& tauchtrks_avsip3d_sig);
    void get_avip2d(vector<Track> trks, const TransientTrackBuilder& ttrkbuilder, reco::Vertex vtx, vector<tuple<double, double, double> >& tausdir, double& tauchtrks_avip2d_val, double& tauchtrks_avip2d_sig, double& tauchtrks_avsip2d_val, double& tauchtrks_avsip2d_sig);
    void get_avip1d(vector<Track> trks, const TransientTrackBuilder& ttrkbuilder, reco::Vertex vtx, vector<tuple<double, double, double> >& tausdir, double& tauchtrks_avip1d_val, double& tauchtrks_avip1d_sig, double& tauchtrks_avsip1d_val, double& tauchtrks_avsip1d_sig);
};
#endif
