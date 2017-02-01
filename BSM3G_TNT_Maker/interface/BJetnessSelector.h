#ifndef __BJetness_MU_H_
#define __BJetness_MU_H_
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
#include "FWCore/Framework/interface/EventSetup.h"
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
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "TMath.h"
#include "RecoBTag/BTagTools/interface/SignedTransverseImpactParameter.h"
#include "BSMFramework/BSM3G_TNT_Maker/interface/TTHbb_eventSelector.h"

using namespace std;
using namespace pat;
using namespace edm;
using namespace reco;
/////
//   Class declaration
/////
class BJetnessSelector : public  baseTree{
  public:
    BJetnessSelector(std::string name, TTree* tree, bool debug, const edm::ParameterSet& cfg, edm::ConsumesCollector && iC);
    ~BJetnessSelector();
    void Fill(const edm::Event& iEvent, const edm::EventSetup& iSetup, int& bjetnesssel_filter);
    void SetBranches();
    void Clear();
    bool isGoodVertex(const reco::Vertex& vtx);
    void JECInitialization();
    void GetJER(pat::Jet jet, float JesSF, float rhoJER, bool AK4PFchs, float &JERScaleFactor, float &JERScaleFactorUP, float &JERScaleFactorDOWN);
  private:
    BJetnessSelector(){};

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
    edm::EDGetTokenT<edm::View<reco::GenParticle> > prunedGenToken_;
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
    //   TTH variables
    /////
    //Evt selection
    int BJetness_isSingleLepton, BJetness_isDoubleLepton;
    //Gen info
    int BJetness_ngenbh, BJetness_ngenbt, BJetness_ngenb, BJetness_ngenc;
    vector<double> BJetness_partonFlavour, BJetness_hadronFlavour;
    //Kinematics and csv
    vector<double> BJetness_numjet;
    vector<double> BJetness_jetpt, BJetness_jeteta, BJetness_jetphi, BJetness_jetenergy;
    vector<double> BJetness_jetcsv, BJetness_pfJetProbabilityBJetTags, BJetness_pfCombinedMVAV2BJetTags, BJetness_pfCombinedCvsLJetTags, BJetness_pfCombinedCvsBJetTags;
    vector<double> BJetness_pt, BJetness_eta, BJetness_phi, BJetness_en, BJetness_ptOVen;
    //PFcand info
    vector<double> BJetness_jetschpvass, BJetness_jetschfrompv, BJetness_jetschip3dval, BJetness_jetschip3dsig, BJetness_jetschip2dval, BJetness_jetschip2dsig, BJetness_jetschisgoodtrk, BJetness_jetschtrkpur, BJetness_jetschpt, BJetness_jetschen;
    //Num_of_trks
    vector<double> BJetness_num_pdgid_eles, BJetness_num_soft_eles, BJetness_num_vetonoipnoiso_eles, BJetness_num_loosenoipnoiso_eles, BJetness_num_veto_eles, BJetness_num_loose_eles, BJetness_num_medium_eles, BJetness_num_tight_eles, BJetness_num_mvatrig_eles, BJetness_num_mvanontrig_eles, BJetness_num_mvatrigwp90_eles, BJetness_num_mvanontrigwp90_eles, BJetness_num_heep_eles, BJetness_num_pdgid_mus, BJetness_num_loose_mus, BJetness_num_soft_mus, BJetness_num_medium_mus, BJetness_num_tight_mus, BJetness_num_highpt_mus, BJetness_num_POGisGood_mus;
    vector<double> BJetness_numjettrks, BJetness_numjettrkspv, BJetness_numjettrksnopv;
    vector<double> BJetness_npvTrkOVcollTrk, BJetness_pvTrkOVcollTrk, BJetness_npvTrkOVpvTrk;
    vector<double> BJetness_npvPtOVcollPt, BJetness_pvPtOVcollPt, BJetness_npvPtOVpvPt;
    //Trk prop rel to jet dir
    vector<double> BJetness_avprel, BJetness_avppar, BJetness_avetarel, BJetness_avetapar, BJetness_avdr, BJetness_avpreljetpt, BJetness_avpreljeten, BJetness_avpparjetpt, BJetness_avpparjeten;
    //Two_trk_info
    vector<double> BJetness_avnum2v, BJetness_avnumno2v, BJetness_avdca3d2t, BJetness_avdca3dno2t, BJetness_avdca3d, BJetness_avdca2d2t, BJetness_avdca2dno2t, BJetness_avdca2d;
    //chi2
    vector<double> BJetness_chi2;
    //ImpactParameter
    vector<double> BJetness_avip3d_val, BJetness_avip3d_sig, BJetness_avsip3d_val, BJetness_avsip3d_sig, BJetness_numip3dpos, BJetness_numip3dneg, BJetness_avip2d_val, BJetness_avip2d_sig, BJetness_avsip2d_val, BJetness_avsip2d_sig, BJetness_numip2dpos, BJetness_numip2dneg, BJetness_avip1d_val, BJetness_avip1d_sig, BJetness_avsip1d_val, BJetness_avsip1d_sig;
    /////
    //   TTH methods
    /////
    //bool is_soft_muon(const pat::PackedCandidate &jcand, edm::View<pat::Muon>& );
    bool is_loose_muon(const pat::Muon& mu, const reco::Vertex& vtx);
    bool is_tight_muon(const pat::Muon& mu, const reco::Vertex& vtx);
    double rel_iso_dbc_mu(const pat::Muon& lepton);
    double rel_iso_dbc_ele(const pat::Electron& lepton, double rhopog);
    double get_effarea(double eta);
    bool is_good_jet(const pat::Jet &j, double rho, double rhoJER, int vtxsize);
    bool is_loosePOG_jetmuon(const pat::PackedCandidate &jcand, edm::Handle<edm::View<pat::Muon> > muon_h);
    bool is_softLep_jetelectron(const pat::Electron &lele);
    bool is_vetoPOGNoIPNoIso_jetelectron(const pat::Electron &lele);
    bool is_loosePOGNoIPNoIso_jetelectron(const pat::Electron &lele);
    bool is_loose_electron(const pat::Electron& ele, double rhopog);//, const reco::Vertex& vtx);
    bool is_tight_electron(const pat::Electron& ele, double rhopog);//, const reco::Vertex& vtx);
    bool is_goodtrk(Track trk,const reco::Vertex& vtx);
    TransientTrack get_ttrk(Track trk, const TransientTrackBuilder& ttrkbuilder);
    vector<TransientTrack> get_ttrks(vector<Track> trks, const TransientTrackBuilder& ttrkbuilder);
    TransientVertex get_tv(vector<Track> trks, const TransientTrackBuilder& ttrkbuilder);
    void get_jettrks(const pat::Jet& jet, const reco::Vertex& vtx, const TransientTrackBuilder& ttrkbuilder, vector<Track>& jetchtrks, vector<Track>& jetchtrkspv, vector<Track>& jetchtrksnpv, vector<tuple<double, double, double> >& jetsdir, const edm::Event& iEvent, edm::Handle<edm::View<pat::Electron> > electron_pat, edm::Handle<edm::View<pat::Muon> > muon_h, double& bjetness_num_pdgid_eles, double& bjetness_num_soft_eles, double& bjetness_num_vetonoipnoiso_eles, double& bjetness_num_loosenoipnoiso_eles, double& bjetness_num_veto_eles, double& bjetness_num_loose_eles, double& bjetness_num_medium_eles, double& bjetness_num_tight_eles, double& bjetness_num_mvatrig_eles, double& bjetness_num_mvanontrig_eles, double& bjetness_num_mvatrigwp90_eles, double& bjetness_num_mvanontrigwp90_eles, double& bjetness_num_heep_eles, double& bjetness_num_pdgid_mus, double& bjetness_num_loose_mus, double& bjetness_num_soft_mus, double& bjetness_num_medium_mus, double& bjetness_num_tight_mus, double& bjetness_num_highpt_mus, double& bjetness_num_POGisGood_mus);
    void get_avreljet(vector<Track> trks, vector<tuple<double, double, double> >& jetsdir, double& jetchtrks_avprel, double& jetchtrks_avppar, double& jetchtrks_avetarel, double& jetchtrks_avetapar, double& jetchtrks_avdr);
    void get_chi2(vector<Track> trks, const TransientTrackBuilder& ttrkbuilder, double& chi2red);
    void get_2trksinfo(vector<Track> trks, const TransientTrackBuilder& ttrkbuilder, double& num2v, double& numno2v, double& dca3d2t, double& dca3dno2t, double& dca2d2t, double& dca2dno2t);
    pair<double,double> dca2trks(Track tkA, Track tkB, const TransientTrackBuilder& ttrkbuilder);
    void get_avip3d(vector<Track> trks, const TransientTrackBuilder& ttrkbuilder, reco::Vertex vtx, vector<tuple<double, double, double> >& jetsdir, double& jetchtrks_avip3d_val, double& jetchtrks_avip3d_sig, double& jetchtrks_avsip3d_val, double& jetchtrks_avsip3d_sig, double& jetchtrks_numip3dpos, double& jetchtrks_numip3dneg);
    void get_avip2d(vector<Track> trks, const TransientTrackBuilder& ttrkbuilder, reco::Vertex vtx, vector<tuple<double, double, double> >& jetsdir, double& jetchtrks_avip2d_val, double& jetchtrks_avip2d_sig, double& jetchtrks_avsip2d_val, double& jetchtrks_avsip2d_sig, double& jetchtrks_numip2dpos, double& jetchtrks_numip2dneg);
    void get_avip1d(vector<Track> trks, const TransientTrackBuilder& ttrkbuilder, reco::Vertex vtx, vector<tuple<double, double, double> >& jetsdir, double& jetchtrks_avip1d_val, double& jetchtrks_avip1d_sig, double& jetchtrks_avsip1d_val, double& jetchtrks_avsip1d_sig);
};
#endif
