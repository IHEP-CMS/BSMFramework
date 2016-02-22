#include "BSMFramework/BSM3G_TNT_Maker/interface/TauJetnessSelector.h"
#include "BSMFramework/BSM3G_TNT_Maker/interface/ElectronPatSelector.h"
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
KalmanVertexFitter vtxFitterTau(true);
using namespace std;
using namespace pat;
using namespace edm;
using namespace reco;
TauJetnessSelector::TauJetnessSelector(std::string name, TTree* tree, bool debug, const pset& iConfig, edm::ConsumesCollector && ic):
  baseTree(name,tree,debug),
  electronLooseIdMapToken_(ic.consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electronLooseIdMap"))),
  eleMVATrigIdMapToken_(ic.consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleMVATrigIdMap"))),
  elemvaValuesMapToken_Trig_(ic.consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("elemvaValuesMap_Trig"))),
  elemvaCategoriesMapToken_Trig_(ic.consumes<edm::ValueMap<int> >(iConfig.getParameter<edm::InputTag>("elemvaCategoriesMap_Trig")))
{
  vtx_h_               = ic.consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"));
  electron_pat_        = ic.consumes<edm::View<pat::Electron> >(iConfig.getParameter<edm::InputTag>("patElectrons"));
  muon_h_              = ic.consumes<edm::View<pat::Muon> >(iConfig.getParameter<edm::InputTag>("muons"));
  taus_                = ic.consumes<edm::View<pat::Tau> >(iConfig.getParameter<edm::InputTag>("taus"));
  jets_                = ic.consumes<pat::JetCollection >(iConfig.getParameter<edm::InputTag>("jets"));
  _Tau_pt_min          = iConfig.getParameter<double>("Tau_pt_min");
  _Tau_eta_max         = iConfig.getParameter<double>("Tau_eta_max");
  rhopogHandle_        = ic.consumes<double>(edm::InputTag("fixedGridRhoFastjetAll"));
  _is_data = iConfig.getParameter<bool>("is_data");
  if(!_is_data) prunedGenToken_ = ic.consumes<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("pruned"));
  SetBranches();
}
TauJetnessSelector::~TauJetnessSelector(){
  delete tree_;
}
void TauJetnessSelector::Fill(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  Clear();
  /////
  //   Recall collections
  ///// 
  edm::Handle<reco::VertexCollection> vtx_h;
  iEvent.getByToken(vtx_h_, vtx_h);
  edm::Handle<edm::View<pat::Electron> > electron_pat;
  iEvent.getByToken(electron_pat_, electron_pat);
  edm::Handle<edm::View<pat::Muon> > muon_h;
  iEvent.getByToken(muon_h_, muon_h);
  edm::Handle<edm::View<pat::Tau> > taus;
  iEvent.getByToken(taus_, taus);
  edm::Handle<pat::JetCollection> jets;                                       
  iEvent.getByToken(jets_, jets);                                         
  edm::Handle<double> rhopogHandle;
  iEvent.getByToken(rhopogHandle_,rhopogHandle);
  //double rhopog = *rhopogHandle;
  edm::Handle<edm::ValueMap<bool> > loose_id_decisions;
  edm::Handle<edm::ValueMap<bool>  > mvatrig_id_decisions;
  edm::Handle<edm::ValueMap<float> > elemvaValues_Trig;
  edm::Handle<edm::ValueMap<int> >   elemvaCategories_Trig;
  iEvent.getByToken(electronLooseIdMapToken_,loose_id_decisions);
  iEvent.getByToken(eleMVATrigIdMapToken_,     mvatrig_id_decisions);
  iEvent.getByToken(elemvaValuesMapToken_Trig_,        elemvaValues_Trig);
  iEvent.getByToken(elemvaCategoriesMapToken_Trig_,    elemvaCategories_Trig);
  edm::ESHandle<TransientTrackBuilder> ttrkbuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",ttrkbuilder);
  /////
  //   Require a good vertex 
  ///// 
  if(vtx_h->empty()) return; // skip the event if no PV found
  const reco::Vertex &PV = vtx_h->front();
  /*
  /////
  //   Look for loose muons, electrons to clean taus
  /////
  vector<const reco::Candidate*> looseleps;
  //Muons
  for(const pat::Muon &mu : *muon_h){
    if(!is_loose_muon(mu,PV)) continue;
    looseleps.push_back((const reco::Candidate*)&mu);
  }
  //Electrons
  //for(const pat::Electron &ele : *electron_pat){
  for(edm::View<pat::Electron>::const_iterator ele = electron_pat->begin(); ele != electron_pat->end(); ele++){
    const Ptr<pat::Electron> elPtr(electron_pat, ele - electron_pat->begin() );
    bool isPassMvatrig = (*mvatrig_id_decisions)[ elPtr ];
    //float mvaval_Trig  = (*elemvaValues_Trig)[ elPtr ];
    //float mvacat_Trig  = (*elemvaCategories_Trig)[ elPtr ];
    if(!(is_loose_electron(*ele,rhopog) && isPassMvatrig)) continue;
    //bool matchelemu = false;
    //for(uint gl=0; gl<looseleps.size(); gl++) if(deltaR(looseleps[gl]->p4(),ele.p4())<0.3) matchelemu = true;
    //if(matchelemu) continue;
    const pat::Electron &lele = *ele;
    looseleps.push_back((const reco::Candidate*)&lele);
  }
  */
  /////
  //   Get the good taus of the event
  /////  
  int tau_pos = 0; //This counter helps to order taus
  int tau_num = 0; //This counter accounts for the number of good taus in the events
                   //The definition of good tau in the event must be the same of the tau related analyses
  vector<int> tau_pos_v;
  for(edm::View<pat::Tau>::const_iterator tau = taus->begin(); tau != taus->end(); tau++){
    //Selection 
    if(!(tau->pt()>_Tau_pt_min && fabs(tau->eta())<_Tau_eta_max &&
         tau->tauID("againstMuonLoose3") &&
         tau->tauID("againstElectronLooseMVA6") &&
         tau->tauID("byMediumIsolationMVArun2v1DBnewDMwLT")
      )){tau_pos++; continue;}
    //Tau counters
    tau_pos_v.push_back(tau_pos);
    tau_pos++;
    tau_num++;
    //Tau info
    TauJetness_taupt.push_back(tau->pt());
    TauJetness_taueta.push_back(tau->eta());
    TauJetness_tauphi.push_back(tau->phi());
    TauJetness_tauenergy.push_back(tau->energy());
    TauJetness_charge.push_back(tau->charge());
    TauJetness_againstMuonTight3.push_back(tau->tauID("againstMuonTight3"));
    TauJetness_againstElectronMediumMVA6.push_back(tau->tauID("againstElectronMediumMVA6"));
    TauJetness_againstElectronTightMVA6.push_back(tau->tauID("againstElectronTightMVA6"));
    TauJetness_byTightIsolationMVArun2v1DBnewDMwLT.push_back(tau->tauID("byTightIsolationMVArun2v1DBnewDMwLT"));
    TauJetness_byVTightIsolationMVArun2v1DBnewDMwLT.push_back(tau->tauID("byVTightIsolationMVArun2v1DBnewDMwLT"));
  }
  if(tau_num>=2){
    TauJetness_numtau.push_back(tau_num);
    /////
    //   Gen info
    /////
    if(!_is_data){
      TauJetness_ngentau  = 0;
      Handle<edm::View<reco::GenParticle> > pruned;
      iEvent.getByToken(prunedGenToken_, pruned);
      for(size_t i=0; i<pruned->size(); i++){
	if(abs((*pruned)[i].pdgId())==15){
	  const Candidate * gentau = &(*pruned)[i];
	  if(gentau->mother(0)!=0 && abs(gentau->mother(0)->pdgId())!=15) TauJetness_ngentau++; 
	}   
      }
    }
    //Access TauJetness info 
    vector<Track> tauschtrks; tauschtrks.clear();
    vector<Track> tauschtrkspv; tauschtrkspv.clear();
    vector<Track> tauschtrksnpv; tauschtrksnpv.clear();
    vector<tuple<double, double, double> > tausdir; tausdir.clear();
    double tauness_num_pdgid_eles = 0;
    double tauness_num_pdgid_mus  = 0;
    double tauness_num_soft_eles  = 0;
    double tauness_num_vetonoipnoiso_eles  = 0;
    double tauness_num_loosenoipnoiso_eles = 0;
    double tauness_num_loose_mus  = 0;
    int maxtaunum = 2; //Pass this value from python at some point
                       //The values will have to be chosen depending on the definition of the categories (for now ditau analyses)
    if(tau_num<maxtaunum) maxtaunum = tau_num;
    for(int jn=0; jn<maxtaunum; jn++){
      const pat::Tau & pftau = (*taus)[tau_pos_v[jn]];
      //cout<<jn<<setw(20)<<pftau.pt()<<setw(20)<<pftau.eta()<<setw(20)<<pftau.phi()<<endl;
      //At some point you will have to make a stronger association to find the pat::jet of the tau
      //I believe you can access the tau jet from miniAOD somehow, but not sure how (to be investigated)
      //const PFJetRef&  pfjetref = pftau.pfJetRef();
      //cout<<jn<<setw(20)<<pfjetref->pt()<<setw(20)<<pfjetref->eta()<<setw(20)<<pfjetref->phi()<<endl;
      double dR  = 10;
      int jetpos = 0; //Note it means in case of no association by dR, you will take the first jet as the one associated to the tau (to be investigated)
      int taujetpos = 0;
      for(const pat::Jet &j : *jets){
       if(deltaR(pftau.p4(),j.p4())<dR){
        dR = deltaR(pftau.p4(),j.p4());
        taujetpos = jetpos;
       }
       jetpos++;
      }
      const pat::Jet & j = (*jets)[taujetpos];;
      //cout<<jn<<setw(20)<<j.pt()<<setw(20)<<j.eta()<<setw(20)<<j.phi()<<endl;     
      get_tautrks(j, PV, *ttrkbuilder, tauschtrks, tauschtrkspv, tauschtrksnpv, tausdir,
                  electron_pat, muon_h,
                  tauness_num_pdgid_eles, tauness_num_pdgid_mus, tauness_num_soft_eles, tauness_num_vetonoipnoiso_eles, tauness_num_loosenoipnoiso_eles, tauness_num_loose_mus 
                 );
    }  
    /////
    //   Get info to evaluate the event TauJetness
    /////
    //Num_of_trks
    //cout<<setw(20)<<"numlep"<<setw(20)<<tauness_num_pdgid_eles+tauness_num_pdgid_mus<<setw(20)<<tauness_num_pdgid_eles<<setw(20)<<tauness_num_pdgid_mus<<setw(20)<<tauness_num_soft_eles+tauness_num_loose_mus<<setw(20)<<tauness_num_soft_eles<<setw(20)<<tauness_num_vetonoipnoiso_eles+tauness_num_loose_mus<<setw(20)<<tauness_num_vetonoipnoiso_eles<<setw(20)<<tauness_num_loosenoipnoiso_eles+tauness_num_loose_mus<<setw(20)<<tauness_num_loosenoipnoiso_eles<<setw(20)<<tauness_num_loose_mus<<endl;
    TauJetness_num_pdgid_leps.push_back(tauness_num_pdgid_eles+tauness_num_pdgid_mus);
    TauJetness_num_pdgid_eles.push_back(tauness_num_pdgid_eles);
    TauJetness_num_pdgid_mus.push_back(tauness_num_pdgid_mus);
    TauJetness_num_soft_leps.push_back(tauness_num_soft_eles+tauness_num_loose_mus);
    TauJetness_num_soft_eles.push_back(tauness_num_soft_eles);
    TauJetness_num_vetonoipnoiso_leps.push_back(tauness_num_vetonoipnoiso_eles+tauness_num_loose_mus);
    TauJetness_num_vetonoipnoiso_eles.push_back(tauness_num_vetonoipnoiso_eles);
    TauJetness_num_loosenoipnoiso_leps.push_back(tauness_num_loosenoipnoiso_eles+tauness_num_loose_mus);
    TauJetness_num_loosenoipnoiso_eles.push_back(tauness_num_loosenoipnoiso_eles);
    TauJetness_num_loose_mus.push_back(tauness_num_loose_mus);
    TauJetness_numtautrks.push_back(tauschtrks.size());
    TauJetness_numtautrkspv.push_back(tauschtrkspv.size());
    TauJetness_numtautrksnopv.push_back(tauschtrksnpv.size());
    //cout<<setw(20)<<"Num_of_trks"<<setw(20)<<tauschtrks.size()<<setw(20)<<tauschtrkspv.size()<<setw(20)<<tauschtrksnpv.size()<<endl;
    if(tauschtrks.size()!=0){
      //cout<<setw(20)<<"nonPVTrk/CollTrk"<<setw(20)<<"PVTrk/CollTrk"<<setw(20)<<"nonPVTrk/PVTrk"<<endl;
      //cout<<double(tauschtrksnpv.size())/double(tauschtrks.size())<<setw(20)<<double(tauschtrkspv.size())/double(tauschtrks.size())<<setw(20)<<double(tauschtrksnpv.size())/double(tauschtrkspv.size())<<endl;
      TauJetness_npvTrkOVcollTrk.push_back(double(tauschtrksnpv.size())/double(tauschtrks.size()));
      TauJetness_pvTrkOVcollTrk.push_back(double(tauschtrkspv.size())/double(tauschtrks.size()));
      if(tauschtrkspv.size()!=0) TauJetness_npvTrkOVpvTrk.push_back(double(tauschtrksnpv.size())/double(tauschtrkspv.size()));
      else                       TauJetness_npvTrkOVpvTrk.push_back(-997);
      //cout<<setw(20)<<"nonPVPt/CollPT"<<setw(20)<<"PVPt/CollPT"<<setw(20)<<"nonPVPt/PVPt"<<endl;
      double pttautrks    = 0;
      for(uint jt=0; jt<tauschtrks.size(); jt++) pttautrks += tauschtrks[jt].pt();
      double pttautrkspv    = 0;
      for(uint jt=0; jt<tauschtrkspv.size(); jt++) pttautrkspv += tauschtrkspv[jt].pt();
      double pttautrksnpv    = 0;
      for(uint jt=0; jt<tauschtrksnpv.size(); jt++) pttautrksnpv += tauschtrksnpv[jt].pt();
      //cout<<setw(20)<<pttautrksnpv/pttautrks<<setw(20)<<pttautrkspv/pttautrks<<setw(20)<<pttautrksnpv/pttautrkspv<<endl;
      if(pttautrks!=0) TauJetness_npvPtOVcollPt.push_back(pttautrksnpv/pttautrks);
      else             TauJetness_npvPtOVcollPt.push_back(-997);
      //cout<<"TauJetness_npvTrkOVcollTrk "<<double(tauschtrksnpv.size())/double(tauschtrks.size())<<" "<<double(tauschtrksnpv.size())<<" "<<double(tauschtrks.size())<<endl;
      //cout<<"TauJetness_npvPtOVcollPt "<<pttautrksnpv/pttautrks<<" "<<pttautrksnpv<<" "<<pttautrks<<endl;
      if(pttautrks!=0) TauJetness_pvPtOVcollPt.push_back(pttautrkspv/pttautrks);
      else             TauJetness_pvPtOVcollPt.push_back(-997);
      if(pttautrkspv!=0) TauJetness_npvPtOVpvPt.push_back(pttautrksnpv/pttautrkspv);    
      else               TauJetness_npvPtOVpvPt.push_back(-997);
      //Two_trk_info (we may want to use one function per variable. Probably more clear, but slower!)
      double tauchtrks_num2v     = 0;
      double tauchtrks_numno2v   = 0;
      double tauchtrks_dca3d2t   = 0;
      double tauchtrks_dca3dno2t = 0;
      double tauchtrks_dca2d2t   = 0;
      double tauchtrks_dca2dno2t = 0;
      get_2trksinfo(tauschtrks, *ttrkbuilder, tauchtrks_num2v, tauchtrks_numno2v, tauchtrks_dca3d2t, tauchtrks_dca3dno2t, tauchtrks_dca2d2t, tauchtrks_dca2dno2t);
      if((tauchtrks_num2v+tauchtrks_numno2v)!=0){
        double twotrksinfo_valtemp = 0;
        twotrksinfo_valtemp = tauchtrks_num2v/(tauchtrks_num2v+tauchtrks_numno2v);
        if(twotrksinfo_valtemp==twotrksinfo_valtemp) TauJetness_avnum2v.push_back(twotrksinfo_valtemp);
        else                                         TauJetness_avnum2v.push_back(-996);
        twotrksinfo_valtemp = tauchtrks_numno2v/(tauchtrks_num2v+tauchtrks_numno2v);
        if(twotrksinfo_valtemp==twotrksinfo_valtemp) TauJetness_avnumno2v.push_back(twotrksinfo_valtemp);
        else                                         TauJetness_avnumno2v.push_back(-996);
        twotrksinfo_valtemp = tauchtrks_dca3d2t/(tauchtrks_num2v+tauchtrks_numno2v);
        if(twotrksinfo_valtemp==twotrksinfo_valtemp) TauJetness_avdca3d2t.push_back(twotrksinfo_valtemp);
        else                                         TauJetness_avdca3d2t.push_back(-996);
        twotrksinfo_valtemp = tauchtrks_dca3dno2t/(tauchtrks_num2v+tauchtrks_numno2v);
        if(twotrksinfo_valtemp==twotrksinfo_valtemp) TauJetness_avdca3dno2t.push_back(twotrksinfo_valtemp);
        else                                         TauJetness_avdca3dno2t.push_back(-996);
        twotrksinfo_valtemp = (tauchtrks_dca3d2t+tauchtrks_dca3dno2t)/(tauchtrks_num2v+tauchtrks_numno2v);
        if(twotrksinfo_valtemp==twotrksinfo_valtemp) TauJetness_avdca3d.push_back(twotrksinfo_valtemp);
        else                                         TauJetness_avdca3d.push_back(-996);
        twotrksinfo_valtemp = tauchtrks_dca2d2t/(tauchtrks_num2v+tauchtrks_numno2v);
        if(twotrksinfo_valtemp==twotrksinfo_valtemp) TauJetness_avdca2d2t.push_back(twotrksinfo_valtemp);
        else                                         TauJetness_avdca2d2t.push_back(-996);
        twotrksinfo_valtemp = tauchtrks_dca2dno2t/(tauchtrks_num2v+tauchtrks_numno2v);
        if(twotrksinfo_valtemp==twotrksinfo_valtemp) TauJetness_avdca2dno2t.push_back(twotrksinfo_valtemp);
        else                                         TauJetness_avdca2dno2t.push_back(-996);
        twotrksinfo_valtemp = (tauchtrks_dca2d2t+tauchtrks_dca2dno2t)/(tauchtrks_num2v+tauchtrks_numno2v);
        if(twotrksinfo_valtemp==twotrksinfo_valtemp) TauJetness_avdca2d.push_back(twotrksinfo_valtemp);
        else                                         TauJetness_avdca2d.push_back(-996);
      }else{
        TauJetness_avnum2v.push_back(-997);
        TauJetness_avnumno2v.push_back(-997);
        TauJetness_avdca3d2t.push_back(-997);
        TauJetness_avdca3dno2t.push_back(-997);
        TauJetness_avdca3d.push_back(-997);
        TauJetness_avdca2d2t.push_back(-997);
        TauJetness_avdca2dno2t.push_back(-997);
        TauJetness_avdca2d.push_back(-997);
      }
      //cout<<setw(20)<<"Two_trk_info"<<setw(20)<<(tauchtrks_num2v+tauchtrks_numno2v)<<endl;
      //chi2
      double tauchtrks_chi2 = 997;
      get_chi2(tauschtrks, *ttrkbuilder, tauchtrks_chi2);
      TauJetness_chi2.push_back(tauchtrks_chi2);
      //ImpactParameter  
      double ip_valtemp = 0;
      double tauchtrks_avip3d_val  = 0;  
      double tauchtrks_avip3d_sig  = 0;  
      double tauchtrks_avsip3d_val = 0;  
      double tauchtrks_avsip3d_sig = 0;  
      get_avip3d(tauschtrks, *ttrkbuilder, PV, tausdir, tauchtrks_avip3d_val,tauchtrks_avip3d_sig,tauchtrks_avsip3d_val,tauchtrks_avsip3d_sig);
      ip_valtemp = tauchtrks_avip3d_val/tauschtrks.size();
      if(ip_valtemp==ip_valtemp) TauJetness_avip3d_val.push_back(ip_valtemp);
      else                       TauJetness_avip3d_val.push_back(-996);
      ip_valtemp = tauchtrks_avip3d_sig/tauschtrks.size();
      if(ip_valtemp==ip_valtemp) TauJetness_avip3d_sig.push_back(ip_valtemp);
      else                       TauJetness_avip3d_sig.push_back(-996); 
      ip_valtemp = tauchtrks_avsip3d_val/tauschtrks.size();
      if(ip_valtemp==ip_valtemp) TauJetness_avsip3d_val.push_back(ip_valtemp);
      else                       TauJetness_avsip3d_val.push_back(-996); 
      ip_valtemp = tauchtrks_avsip3d_sig/tauschtrks.size();
      if(ip_valtemp==ip_valtemp) TauJetness_avsip3d_sig.push_back(ip_valtemp);
      else                       TauJetness_avsip3d_sig.push_back(-996);
      //cout<<setw(20)<<"ImpactParameter"<<setw(20)<<tauchtrks_avip3d_val<<setw(20)<<tauchtrks_avip3d_sig<<setw(20)<<tauchtrks_avsip3d_val<<setw(20)<<tauchtrks_avsip3d_sig<<endl;
      double tauchtrks_avip2d_val  = 0;  
      double tauchtrks_avip2d_sig  = 0;  
      double tauchtrks_avsip2d_val = 0;  
      double tauchtrks_avsip2d_sig = 0;  
      get_avip2d(tauschtrks, *ttrkbuilder, PV, tausdir, tauchtrks_avip2d_val,tauchtrks_avip2d_sig,tauchtrks_avsip2d_val,tauchtrks_avsip2d_sig);
      ip_valtemp = tauchtrks_avip2d_val/tauschtrks.size();
      if(ip_valtemp==ip_valtemp) TauJetness_avip2d_val.push_back(ip_valtemp);
      else                       TauJetness_avip2d_val.push_back(-996);
      ip_valtemp = tauchtrks_avip2d_sig/tauschtrks.size();
      if(ip_valtemp==ip_valtemp) TauJetness_avip2d_sig.push_back(ip_valtemp);
      else                       TauJetness_avip2d_sig.push_back(-996);
      ip_valtemp = tauchtrks_avsip2d_val/tauschtrks.size();
      if(ip_valtemp==ip_valtemp) TauJetness_avsip2d_val.push_back(ip_valtemp);
      else                       TauJetness_avsip2d_val.push_back(-996);
      ip_valtemp = tauchtrks_avsip2d_sig/tauschtrks.size();
      if(ip_valtemp==ip_valtemp) TauJetness_avsip2d_sig.push_back(ip_valtemp);
      else                       TauJetness_avsip2d_sig.push_back(-996);
      double tauchtrks_avip1d_val  = 0;  
      double tauchtrks_avip1d_sig  = 0;  
      double tauchtrks_avsip1d_val = 0;  
      double tauchtrks_avsip1d_sig = 0;  
      get_avip1d(tauschtrks, *ttrkbuilder, PV, tausdir, tauchtrks_avip1d_val,tauchtrks_avip1d_sig,tauchtrks_avsip1d_val,tauchtrks_avsip1d_sig);
      ip_valtemp = tauchtrks_avip1d_val/tauschtrks.size();
      if(ip_valtemp==ip_valtemp) TauJetness_avip1d_val.push_back(ip_valtemp);
      else                       TauJetness_avip1d_val.push_back(-996);
      ip_valtemp = tauchtrks_avip1d_sig/tauschtrks.size();
      if(ip_valtemp==ip_valtemp) TauJetness_avip1d_sig.push_back(ip_valtemp);
      else                       TauJetness_avip1d_sig.push_back(-996);
      ip_valtemp = tauchtrks_avsip1d_val/tauschtrks.size();
      if(ip_valtemp==ip_valtemp) TauJetness_avsip1d_val.push_back(ip_valtemp);
      else                       TauJetness_avsip1d_val.push_back(-996);
      ip_valtemp = tauchtrks_avsip1d_sig/tauschtrks.size();
      if(ip_valtemp==ip_valtemp) TauJetness_avsip1d_sig.push_back(ip_valtemp);
      else                       TauJetness_avsip1d_sig.push_back(-996);
    }else{//if(tauschtrks.size()!=0)
      //Num trk pt
      TauJetness_npvTrkOVcollTrk.push_back(-998);
      TauJetness_pvTrkOVcollTrk.push_back(-998);
      TauJetness_npvTrkOVpvTrk.push_back(-998);
      TauJetness_npvPtOVcollPt.push_back(-998);
      TauJetness_pvPtOVcollPt.push_back(-998);
      TauJetness_npvPtOVpvPt.push_back(-998);
      //Two_trk_info (we may want to use one function per variable. Probably more clear, but slower!)
      TauJetness_avnum2v.push_back(-998);
      TauJetness_avnumno2v.push_back(-998);
      TauJetness_avdca3d2t.push_back(-998);
      TauJetness_avdca3dno2t.push_back(-998);
      TauJetness_avdca3d.push_back(-998);
      TauJetness_avdca2d2t.push_back(-998);
      TauJetness_avdca2dno2t.push_back(-998);
      TauJetness_avdca2d.push_back(-998);
      //chi2
      TauJetness_chi2.push_back(-998);
      //ImpactParameter  
      TauJetness_avip3d_val.push_back(-998);
      TauJetness_avip3d_sig.push_back(-998);
      TauJetness_avsip3d_val.push_back(-998);
      TauJetness_avsip3d_sig.push_back(-998);
      TauJetness_avip2d_val.push_back(-998);
      TauJetness_avip2d_sig.push_back(-998);
      TauJetness_avsip2d_val.push_back(-998);
      TauJetness_avsip2d_sig.push_back(-998);
      TauJetness_avip1d_val.push_back(-998);
      TauJetness_avip1d_sig.push_back(-998);
      TauJetness_avsip1d_val.push_back(-998);
      TauJetness_avsip1d_sig.push_back(-998);
    }
 }else{//if(tau_num!=0)
   TauJetness_numtau.push_back(-999);
   //Kinematics and csv 
   TauJetness_taupt.push_back(-999);
   TauJetness_taueta.push_back(-999);
   TauJetness_tauphi.push_back(-999);
   TauJetness_tauenergy.push_back(-999);
   TauJetness_charge.push_back(-999);
   TauJetness_againstMuonTight3.push_back(-999);
   TauJetness_againstElectronMediumMVA6.push_back(-999);
   TauJetness_againstElectronTightMVA6.push_back(-999);
   TauJetness_byTightIsolationMVArun2v1DBnewDMwLT.push_back(-999);
   TauJetness_byVTightIsolationMVArun2v1DBnewDMwLT.push_back(-999);
   //Num_of_trks
   TauJetness_num_pdgid_leps.push_back(-999);
   TauJetness_num_pdgid_eles.push_back(-999);
   TauJetness_num_pdgid_mus.push_back(-999);
   TauJetness_num_soft_leps.push_back(-999);
   TauJetness_num_soft_eles.push_back(-999);
   TauJetness_num_vetonoipnoiso_leps.push_back(-999);
   TauJetness_num_vetonoipnoiso_eles.push_back(-999);
   TauJetness_num_loosenoipnoiso_leps.push_back(-999);
   TauJetness_num_loosenoipnoiso_eles.push_back(-999);
   TauJetness_num_loose_mus.push_back(-999);
   TauJetness_numtautrks.push_back(-999);
   TauJetness_numtautrkspv.push_back(-999);
   TauJetness_numtautrksnopv.push_back(-999);
   TauJetness_npvTrkOVcollTrk.push_back(-999);
   TauJetness_pvTrkOVcollTrk.push_back(-999);
   TauJetness_npvTrkOVpvTrk.push_back(-999);
   TauJetness_npvPtOVcollPt.push_back(-999);
   TauJetness_pvPtOVcollPt.push_back(-999);
   TauJetness_npvPtOVpvPt.push_back(-999);
   //Two_trk_info (we may want to use one function per variable. Probably more clear, but slower!)
   TauJetness_avnum2v.push_back(-999);
   TauJetness_avnumno2v.push_back(-999);
   TauJetness_avdca3d2t.push_back(-999);
   TauJetness_avdca3dno2t.push_back(-999);
   TauJetness_avdca3d.push_back(-999);
   TauJetness_avdca2d2t.push_back(-999);
   TauJetness_avdca2dno2t.push_back(-999);
   TauJetness_avdca2d.push_back(-999);
   //chi2
   TauJetness_chi2.push_back(-999);
   //ImpactParameter  
   TauJetness_avip3d_val.push_back(-999);
   TauJetness_avip3d_sig.push_back(-999);
   TauJetness_avsip3d_val.push_back(-999);
   TauJetness_avsip3d_sig.push_back(-999);
   TauJetness_avip2d_val.push_back(-999);
   TauJetness_avip2d_sig.push_back(-999);
   TauJetness_avsip2d_val.push_back(-999);
   TauJetness_avsip2d_sig.push_back(-999);
   TauJetness_avip1d_val.push_back(-999);
   TauJetness_avip1d_sig.push_back(-999);
   TauJetness_avsip1d_val.push_back(-999);
   TauJetness_avsip1d_sig.push_back(-999);
 }
}
void TauJetnessSelector::SetBranches(){
  if(debug_) std::cout<<"setting branches: calling AddBranch of baseTree"<<std::endl;
  //Gen info
  AddBranch(&TauJetness_ngentau                    ,"TauJetness_ngentau");
  AddBranch(&TauJetness_numtau                   ,"TauJetness_numtau");
  //Kinematics and csv 
  AddBranch(&TauJetness_taupt                    ,"TauJetness_taupt");
  AddBranch(&TauJetness_taueta                   ,"TauJetness_taueta");
  AddBranch(&TauJetness_tauphi                   ,"TauJetness_tauphi");
  AddBranch(&TauJetness_tauenergy                ,"TauJetness_tauenergy");
  AddBranch(&TauJetness_charge                   ,"TauJetness_charge");
  AddBranch(&TauJetness_againstMuonTight3                    ,"TauJetness_againstMuonTight3");
  AddBranch(&TauJetness_againstElectronMediumMVA6            ,"TauJetness_againstElectronMediumMVA6");
  AddBranch(&TauJetness_againstElectronTightMVA6             ,"TauJetness_againstElectronTightMVA6");
  AddBranch(&TauJetness_byTightIsolationMVArun2v1DBnewDMwLT  ,"TauJetness_byTightIsolationMVArun2v1DBnewDMwLT");
  AddBranch(&TauJetness_byVTightIsolationMVArun2v1DBnewDMwLT ,"TauJetness_byVTightIsolationMVArun2v1DBnewDMwLT");
  //Num_of_trks
  AddBranch(&TauJetness_num_pdgid_leps           ,"TauJetness_num_pdgid_leps");
  AddBranch(&TauJetness_num_pdgid_eles           ,"TauJetness_num_pdgid_eles");
  AddBranch(&TauJetness_num_pdgid_mus            ,"TauJetness_num_pdgid_mus");
  AddBranch(&TauJetness_num_soft_leps            ,"TauJetness_num_soft_leps");
  AddBranch(&TauJetness_num_soft_eles            ,"TauJetness_num_soft_eles");
  AddBranch(&TauJetness_num_vetonoipnoiso_leps   ,"TauJetness_num_vetonoipnoiso_leps");
  AddBranch(&TauJetness_num_vetonoipnoiso_eles   ,"TauJetness_num_vetonoipnoiso_eles");
  AddBranch(&TauJetness_num_loosenoipnoiso_leps  ,"TauJetness_num_loosenoipnoiso_leps");
  AddBranch(&TauJetness_num_loosenoipnoiso_eles  ,"TauJetness_num_loosenoipnoiso_eles");
  AddBranch(&TauJetness_num_loose_mus            ,"TauJetness_num_loose_mus");
  AddBranch(&TauJetness_numtautrks               ,"TauJetness_numtautrks");
  AddBranch(&TauJetness_numtautrkspv             ,"TauJetness_numtautrkspv");
  AddBranch(&TauJetness_numtautrksnopv           ,"TauJetness_numtautrksnopv");
  AddBranch(&TauJetness_npvTrkOVcollTrk          ,"TauJetness_npvTrkOVcollTrk");
  AddBranch(&TauJetness_pvTrkOVcollTrk           ,"TauJetness_pvTrkOVcollTrk");
  AddBranch(&TauJetness_npvTrkOVpvTrk            ,"TauJetness_npvTrkOVpvTrk");
  AddBranch(&TauJetness_npvPtOVcollPt            ,"TauJetness_npvPtOVcollPt");
  AddBranch(&TauJetness_pvPtOVcollPt             ,"TauJetness_pvPtOVcollPt");
  AddBranch(&TauJetness_npvPtOVpvPt              ,"TauJetness_npvPtOVpvPt");
  //Two_trk_info
  AddBranch(&TauJetness_avnum2v                  ,"TauJetness_avnum2v");
  AddBranch(&TauJetness_avnumno2v                ,"TauJetness_avnumno2v");
  AddBranch(&TauJetness_avdca3d2t                ,"TauJetness_avdca3d2t");
  AddBranch(&TauJetness_avdca3dno2t              ,"TauJetness_avdca3dno2t");
  AddBranch(&TauJetness_avdca3d                  ,"TauJetness_avdca3d");
  AddBranch(&TauJetness_avdca2d2t                ,"TauJetness_avdca2d2t");
  AddBranch(&TauJetness_avdca2dno2t              ,"TauJetness_avdca2dno2t");
  AddBranch(&TauJetness_avdca2d                  ,"TauJetness_avdca2d");
  //chi2
  AddBranch(&TauJetness_chi2                     ,"TauJetness_chi2");
  //ImpactParameter
  AddBranch(&TauJetness_avip3d_val               ,"TauJetness_avip3d_val");
  AddBranch(&TauJetness_avip3d_sig               ,"TauJetness_avip3d_sig");
  AddBranch(&TauJetness_avsip3d_val              ,"TauJetness_avsip3d_val");
  AddBranch(&TauJetness_avsip3d_sig              ,"TauJetness_avsip3d_sig");
  AddBranch(&TauJetness_avip2d_val               ,"TauJetness_avip2d_val");
  AddBranch(&TauJetness_avip2d_sig               ,"TauJetness_avip2d_sig");
  AddBranch(&TauJetness_avsip2d_val              ,"TauJetness_avsip2d_val");
  AddBranch(&TauJetness_avsip2d_sig              ,"TauJetness_avsip2d_sig");
  AddBranch(&TauJetness_avip1d_val               ,"TauJetness_avip1d_val");
  AddBranch(&TauJetness_avip1d_sig               ,"TauJetness_avip1d_sig");
  AddBranch(&TauJetness_avsip1d_val              ,"TauJetness_avsip1d_val");
  AddBranch(&TauJetness_avsip1d_sig              ,"TauJetness_avsip1d_sig");
  if(debug_) std::cout<<"set branches"<<std::endl;
}
void TauJetnessSelector::Clear(){
  //Gen info
  TauJetness_ngentau  = -9999;
  TauJetness_numtau.clear();
  //Kinematics and csv 
  TauJetness_taupt.clear(); 
  TauJetness_taueta.clear();
  TauJetness_tauphi.clear();
  TauJetness_tauenergy.clear();
  TauJetness_charge.clear();
  TauJetness_againstMuonTight3.clear();
  TauJetness_againstElectronMediumMVA6.clear();
  TauJetness_againstElectronTightMVA6.clear();
  TauJetness_byTightIsolationMVArun2v1DBnewDMwLT.clear();
  TauJetness_byVTightIsolationMVArun2v1DBnewDMwLT.clear();
  //Num_of_trks
  TauJetness_num_pdgid_leps.clear();
  TauJetness_num_pdgid_eles.clear();
  TauJetness_num_pdgid_mus.clear();
  TauJetness_num_soft_leps.clear();
  TauJetness_num_soft_eles.clear();
  TauJetness_num_vetonoipnoiso_leps.clear();
  TauJetness_num_vetonoipnoiso_eles.clear();
  TauJetness_num_loosenoipnoiso_leps.clear();
  TauJetness_num_loosenoipnoiso_eles.clear();
  TauJetness_num_loose_mus.clear();
  TauJetness_numtautrks.clear();
  TauJetness_numtautrkspv.clear();
  TauJetness_numtautrksnopv.clear();
  TauJetness_npvTrkOVcollTrk.clear();
  TauJetness_pvTrkOVcollTrk.clear();
  TauJetness_npvTrkOVpvTrk.clear();
  TauJetness_npvPtOVcollPt.clear();
  TauJetness_pvPtOVcollPt.clear();
  TauJetness_npvPtOVpvPt.clear();
  //Two_trk_info
  TauJetness_avnum2v.clear();
  TauJetness_avnumno2v.clear();
  TauJetness_avdca3d2t.clear();
  TauJetness_avdca3dno2t.clear();
  TauJetness_avdca3d.clear();
  TauJetness_avdca2d2t.clear();
  TauJetness_avdca2dno2t.clear();
  TauJetness_avdca2d.clear();
  //chi2
  TauJetness_chi2.clear();
  //ImpactParameter
  TauJetness_avip3d_val.clear();
  TauJetness_avip3d_sig.clear();
  TauJetness_avsip3d_val.clear();
  TauJetness_avsip3d_sig.clear();
  TauJetness_avip2d_val.clear();
  TauJetness_avip2d_sig.clear();
  TauJetness_avsip2d_val.clear();
  TauJetness_avsip2d_sig.clear();
  TauJetness_avip1d_val.clear();
  TauJetness_avip1d_sig.clear();
  TauJetness_avsip1d_val.clear();
  TauJetness_avsip1d_sig.clear();
}
//Look for loose muon
bool TauJetnessSelector::is_loose_muon(const pat::Muon& mu, const reco::Vertex& vtx){
  bool isloosemu = false;
  if(//mu.muonBestTrack().isNonnull() && mu.globalTrack().isNonnull() && mu.innerTrack().isNonnull() && mu.track().isNonnull() &&
    mu.pt()>15 &&
    TMath::Abs(mu.eta()) < 2.4 &&
    mu.isTightMuon(vtx) &&  
    //mu.isPFMuon() && mu.isGlobalMuon() &&
    //fabs(mu.muonBestTrack()->dxy(vtx.position())) < 0.2 && 
    //fabs(mu.muonBestTrack()->dz(vtx.position()))  < 0.5 && 
    //mu.globalTrack()->normalizedChi2() < 10 &&
    //mu.globalTrack()->hitPattern().numberOfValidMuonHits()  > 0 &&
    //mu.innerTrack()->hitPattern().numberOfValidPixelHits()  > 0 &&
    //mu.track()->hitPattern().trackerLayersWithMeasurement() > 5 &&
    //mu.numberOfMatchedStations() > 1 && 
    rel_iso_dbc_mu(mu) < 0.25
    ) isloosemu = true;
  return isloosemu;
}
double TauJetnessSelector::rel_iso_dbc_mu(const pat::Muon& lepton){
  return((lepton.pfIsolationR04().sumChargedHadronPt + max(lepton.pfIsolationR04().sumNeutralHadronEt + lepton.pfIsolationR04().sumPhotonEt - 0.5 * lepton.pfIsolationR04().sumPUPt,0.0))/lepton.pt()
        );
}
bool TauJetnessSelector::is_loose_electron(const pat::Electron& ele, double rhopog){
  bool isele = false;
  if(ele.pt()>15 && TMath::Abs(ele.eta())<2.4 &&
     !(fabs(ele.superCluster()->position().eta()) > 1.4442 && fabs(ele.superCluster()->position().eta()) < 1.5660)){
    if(fabs(ele.superCluster()->position().eta())<1.4442 
      && (ele.full5x5_sigmaIetaIeta()<0.012)                 
      && (ele.hcalOverEcal()<0.09)                           
      && (ele.ecalPFClusterIso()/ele.pt()<0.37)            
      && (ele.hcalPFClusterIso()/ele.pt()<0.25)            
      && (ele.dr03TkSumPt()/ele.pt()<0.18)                 
      && (fabs(ele.deltaEtaSuperClusterTrackAtVtx())<0.0095)
      && (fabs(ele.deltaPhiSuperClusterTrackAtVtx())<0.065)){  
      isele = true;
    }
    if(fabs(ele.superCluster()->position().eta())>1.5660
      && ele.full5x5_sigmaIetaIeta()<0.033
      && ele.hcalOverEcal()<0.09
      && (ele.ecalPFClusterIso()/ele.pt())<0.45
      && (ele.hcalPFClusterIso()/ele.pt())<0.28
      && (ele.dr03TkSumPt()/ele.pt())<0.18){
      isele = true;
    }
  }
  //ele.passConversionVeto()
  return isele; 
}
//Having the vtx in this function is kept for historical reason 
//and it will be kept until we do not have a final electron loose selection
//bool TauJetnessSelector::is_loose_electron(const pat::Electron& ele, const reco::Vertex& vtx){
//  bool islooseele = false;
//  if(ele.pt()>10 && TMath::Abs(ele.eta())<2.4 &&
//     !(abs(ele.superCluster()->position().eta()) > 1.4442 && abs(ele.superCluster()->position().eta()) < 1.5660) &&
//     ele.passConversionVeto()
//    ) islooseele = true;
//  return islooseele; 
//}
double TauJetnessSelector::rel_iso_dbc_ele(const pat::Electron& lepton, double rhopog){
  double effarea          = get_effarea(lepton.superCluster()->position().eta());  
  double SumChHadPt       = lepton.pfIsolationVariables().sumChargedHadronPt;
  double SumNeuHadEt      = lepton.pfIsolationVariables().sumNeutralHadronEt;
  double SumPhotonEt      = lepton.pfIsolationVariables().sumPhotonEt;
  double SumNeutralCorrEt = std::max( 0.0, SumNeuHadEt+SumPhotonEt - rhopog*effarea );
  return ((SumChHadPt + SumNeutralCorrEt)/lepton.pt());
//  return( (lepton.chargedHadronIso() +
//          std::max(0.0, lepton.neutralHadronIso() + lepton.photonIso() - 0.5*lepton.puChargedHadronIso()))/lepton.pt() );
}
double TauJetnessSelector::get_effarea(double eta){
  double effarea = -1;
  if(abs(eta) < 1.0)        effarea = 0.1752;
  else if(abs(eta) < 1.479) effarea = 0.1862;
  else if(abs(eta) < 2.0)   effarea = 0.1411;
  else if(abs(eta) < 2.2)   effarea = 0.1534;
  else if(abs(eta) < 2.3)   effarea = 0.1903;
  else if(abs(eta) < 2.4)   effarea = 0.2243;
  else                      effarea = 0.2687;
  return effarea;
}
//Look for loose electron
//Require good taus
//This function has to be updated in order to select good tau using the TTHbb definition
bool TauJetnessSelector::is_good_tau(const pat::Jet &j){
  bool isgoodtau = true;
  //Acceptance
  if(j.pt() < 30)       isgoodtau = false; //Please note that this requirement is for the SL channel, while for DL channel we require pT > 20! 
  if(fabs(j.eta())>2.4) isgoodtau = false; 
  //ID requirements
  if(j.neutralHadronEnergyFraction() >= 0.99) isgoodtau = false;
  if(j.chargedEmEnergyFraction()     >= 0.99) isgoodtau = false;
  if(j.neutralEmEnergyFraction()     >= 0.99) isgoodtau = false;
  if(j.numberOfDaughters()           <= 1)    isgoodtau = false;
  if(j.chargedHadronEnergyFraction() <= 0.0)  isgoodtau = false;
  if(j.chargedMultiplicity()         <= 0.0)  isgoodtau = false;
  return isgoodtau;
}
bool TauJetnessSelector::is_loosePOG_taumuon(const pat::PackedCandidate &jcand, edm::Handle<edm::View<pat::Muon> > muon_h){
  bool ismu = false;
  for(const pat::Muon &mu : *muon_h){
    if(deltaR(jcand.p4(),mu.p4())<0.1 && fabs(jcand.pt()-mu.pt())/mu.pt()<0.05){
     if(mu.isLooseMuon()) ismu = true;
     if(ismu) break;
    }
  }  
  return ismu;
}
bool TauJetnessSelector::is_softLep_tauelectron(const pat::PackedCandidate &jcand, edm::Handle<edm::View<pat::Electron> > electron_pat, const reco::Vertex& vtx){
  bool isele = false;
  for(edm::View<pat::Electron>::const_iterator ele = electron_pat->begin(); ele != electron_pat->end(); ele++){
    const pat::Electron &lele = *ele; 
    if(deltaR(jcand.p4(),lele.p4())<0.1 && fabs(jcand.pt()-lele.pt())/lele.pt()<0.05 ){
      const HitPattern &hitPattern = lele.gsfTrack().get()->hitPattern();
      uint32_t hit = hitPattern.getHitPattern(HitPattern::TRACK_HITS, 0);
      bool hitCondition = !(HitPattern::validHitFilter(hit) && ((HitPattern::pixelBarrelHitFilter(hit) && HitPattern::getLayer(hit) < 3) || HitPattern::pixelEndcapHitFilter(hit))); 
      if(!hitCondition && lele.passConversionVeto()) isele = true;
      if(isele) break;
    } 
  } 
  return isele; 
}
bool TauJetnessSelector::is_vetoPOGNoIPNoIso_tauelectron(const pat::PackedCandidate &jcand, edm::Handle<edm::View<pat::Electron> > electron_pat, const reco::Vertex& vtx){
  bool isele = false;
  for(edm::View<pat::Electron>::const_iterator ele = electron_pat->begin(); ele != electron_pat->end(); ele++){
    const pat::Electron &lele = *ele; 
    if(deltaR(jcand.p4(),lele.p4())<0.1 && fabs(jcand.pt()-lele.pt())/lele.pt()<0.05){
      double ooEmooP = 999;
      if(lele.ecalEnergy()==0)                   ooEmooP = 1e30;
      else if(!std::isfinite(lele.ecalEnergy())) ooEmooP = 1e30;
      else                                       ooEmooP = fabs(1.0/lele.ecalEnergy() - lele.eSuperClusterOverP()/lele.ecalEnergy() );
      if(lele.full5x5_sigmaIetaIeta()<0.012 && fabs(lele.deltaEtaSuperClusterTrackAtVtx())<0.0126 && fabs(lele.deltaPhiSuperClusterTrackAtVtx())<0.107
         && lele.hcalOverEcal()<0.186 && ooEmooP<0.239
         //&& fabs((-1) * lele.gsfTrack()->dxy(vtx.position()))<0.0227 && fabs(lele.gsfTrack()->dz(vtx.position()))<0.379
         && lele.gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS)<=2 && lele.passConversionVeto()
        ) isele = true;
      if(isele) break;
    } 
  } 
  return isele; 
}
bool TauJetnessSelector::is_loosePOGNoIPNoIso_tauelectron(const pat::PackedCandidate &jcand, edm::Handle<edm::View<pat::Electron> > electron_pat, const reco::Vertex& vtx){
  bool isele = false;
  for(edm::View<pat::Electron>::const_iterator ele = electron_pat->begin(); ele != electron_pat->end(); ele++){
    const pat::Electron &lele = *ele; 
    if(deltaR(jcand.p4(),lele.p4())<0.1 && fabs(jcand.pt()-lele.pt())/lele.pt()<0.05){
      double ooEmooP = 999;
      if(lele.ecalEnergy()==0)                   ooEmooP = 1e30;
      else if(!std::isfinite(lele.ecalEnergy())) ooEmooP = 1e30;
      else                                       ooEmooP = fabs(1.0/lele.ecalEnergy() - lele.eSuperClusterOverP()/lele.ecalEnergy());
      if(!(fabs(lele.superCluster()->position().eta()) > 1.4442 && fabs(lele.superCluster()->position().eta()) < 1.5660)){
        if(fabs(lele.superCluster()->position().eta())<1.4442
          && (lele.full5x5_sigmaIetaIeta()<0.0103)
          && (lele.hcalOverEcal()<0.104)
          && (fabs(lele.deltaEtaSuperClusterTrackAtVtx())<0.0105)
          && (fabs(lele.deltaPhiSuperClusterTrackAtVtx())<0.115)
          && ooEmooP<0.102
          && lele.gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS)<=2
          && lele.passConversionVeto()   
          ){
            isele = true;
        }
        if(fabs(lele.superCluster()->position().eta())>1.5660
          && (lele.full5x5_sigmaIetaIeta()<0.0301)
          && (lele.hcalOverEcal()<0.0897)
          && (fabs(lele.deltaEtaSuperClusterTrackAtVtx())<0.00814)
          && (fabs(lele.deltaPhiSuperClusterTrackAtVtx())<0.182)
          && ooEmooP<0.126
          && lele.gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS)<=1
          && lele.passConversionVeto()   
          ){
            isele = true;
        }
      } 
      if(isele) break;
    } 
  } 
  return isele; 
}
//Check that the track is a good track
bool TauJetnessSelector::is_goodtrk(Track trk,const reco::Vertex& vtx){
 bool isgoodtrk = false;
 if(trk.pt()>1 &&
   trk.hitPattern().numberOfValidHits()>=8 &&
   trk.hitPattern().numberOfValidPixelHits()>=2 &&
   trk.normalizedChi2()<5 &&
   std::abs(trk.dxy(vtx.position()))<0.2 &&
   std::abs(trk.dz(vtx.position()))<17
   ) isgoodtrk = true;
 return isgoodtrk;
}
//Get transient tracks from track
TransientTrack TauJetnessSelector::get_ttrk(Track trk, const TransientTrackBuilder& ttrkbuilder){
  TransientTrack ttrk;
  ttrk = ttrkbuilder.build(&trk);
  return ttrk;
}
vector<TransientTrack> TauJetnessSelector::get_ttrks(vector<Track> trks, const TransientTrackBuilder& ttrkbuilder){
 vector<TransientTrack> ttrks;
 for(uint tr=0; tr<trks.size(); tr++){
  TransientTrack ttrk = ttrkbuilder.build(&trks[tr]);
  ttrks.push_back(ttrk);
 }
 return ttrks;
}
//Get transient pv from tracks
TransientVertex TauJetnessSelector::get_tv(vector<Track> trks, const TransientTrackBuilder& ttrkbuilder){
  TransientVertex tv;
  vector<TransientTrack> ttrks = get_ttrks(trks,ttrkbuilder);
  if(ttrks.size()>=2)  tv = vtxFitterTau.vertex(ttrks);
  return tv;
}
//Get tau tracks info
void TauJetnessSelector::get_tautrks(const pat::Jet& tau, const reco::Vertex& vtx, const TransientTrackBuilder& ttrkbuilder,
                                   vector<Track>& tauchtrks, vector<Track>& tauchtrkspv, vector<Track>& tauchtrksnpv, vector<tuple<double, double, double> >& tausdir,
                                   edm::Handle<edm::View<pat::Electron> > electron_pat, edm::Handle<edm::View<pat::Muon> > muon_h,
                                   double& tauness_num_pdgid_eles, double& tauness_num_pdgid_mus, double& tauness_num_soft_eles, double& tauness_num_vetonoipnoiso_eles, double& tauness_num_loosenoipnoiso_eles, double& tauness_num_loose_mus
                                  ){
  //Access tau daughters
  vector<CandidatePtr> jdaus(tau.daughterPtrVector());
  sort(jdaus.begin(), jdaus.end(), [](const reco::CandidatePtr &p1, const reco::CandidatePtr &p2) {return p1->pt() > p2->pt();});
  for(uint jd=0; jd<jdaus.size(); jd++){ 
    const pat::PackedCandidate &jcand = dynamic_cast<const pat::PackedCandidate &>(*jdaus[jd]);
    if(deltaR(jcand.p4(),tau.p4())>0.4) continue;
    Track trk = Track(jcand.pseudoTrack());
    bool isgoodtrk = is_goodtrk(trk,vtx);
    //Minimal conditions for a track 
    if(isgoodtrk && jcand.charge()!=0 && jcand.fromPV()>1){
      if(fabs(jcand.pdgId())==13){
        tauness_num_pdgid_mus++;
        if(is_loosePOG_taumuon(jcand,muon_h)) tauness_num_loose_mus++;
      } 
      if(fabs(jcand.pdgId())==11){
        tauness_num_pdgid_eles++;
        if(is_softLep_tauelectron(jcand,electron_pat,vtx)) tauness_num_soft_eles++;      
        if(is_vetoPOGNoIPNoIso_tauelectron(jcand,electron_pat,vtx)) tauness_num_vetonoipnoiso_eles++;      
        if(is_loosePOGNoIPNoIso_tauelectron(jcand,electron_pat,vtx)) tauness_num_loosenoipnoiso_eles++;      
      }
      tauchtrks.push_back(trk);
      //Other conditions on tau daughters
      //Using fromPV method
      if(jcand.fromPV()==pat::PackedCandidate::PVUsedInFit){
        tauchtrkspv.push_back(trk);
      }else{
        tauchtrksnpv.push_back(trk);
      }
      //Fill in the tau direction information to keep synchronisation
      tausdir.push_back(make_tuple(tau.px(),tau.py(),tau.pz()));
    }//Ch trks 
  }//Loop on tau daus
}
//Get chi2 information from trk vertex
void TauJetnessSelector::get_chi2(vector<Track> trks, const TransientTrackBuilder& ttrkbuilder, double& chi2red){
  if(trks.size()>=2){
    TransientVertex trks_tv = get_tv(trks, ttrkbuilder);
    if(trks_tv.isValid()) chi2red = trks_tv.totalChiSquared()/trks_tv.degreesOfFreedom(); 
  }
}
//Get Num of two-trk vertices
void TauJetnessSelector::get_2trksinfo(vector<Track> trks, const TransientTrackBuilder& ttrkbuilder, double& num2v, double& numno2v, double& dca3d2t, double& dca3dno2t, double& dca2d2t, double& dca2dno2t){
  double valtemp = 0;
  for(uint t=0; t<trks.size(); t++){
    for(uint t2=t+1; t2<trks.size(); t2++){
      vector<Track> twotrks;
      twotrks.push_back(trks[t]);
      twotrks.push_back(trks[t2]);
      TransientVertex tv = get_tv(twotrks, ttrkbuilder);
      if(tv.isValid() && TMath::Prob(tv.totalChiSquared(),tv.degreesOfFreedom())>0.05){
        num2v++;
        pair<double,double> dca2trks3d2d = dca2trks(trks[t], trks[t2], ttrkbuilder);
        valtemp = dca2trks3d2d.first;
        if(valtemp==valtemp) dca3d2t += dca2trks3d2d.first;
        valtemp = dca2trks3d2d.second;
        if(valtemp==valtemp) dca2d2t += dca2trks3d2d.second;
      }else{
        numno2v++;
        pair<double,double> dca2trks3d2d = dca2trks(trks[t], trks[t2], ttrkbuilder);
        valtemp = dca2trks3d2d.first;
        if(valtemp==valtemp) dca3dno2t += dca2trks3d2d.first;
        valtemp = dca2trks3d2d.second;
        if(valtemp==valtemp) dca2dno2t += dca2trks3d2d.second;
      }
    }
  }
}
//DCA between two trks
pair<double,double> TauJetnessSelector::dca2trks(Track tkA, Track tkB, const TransientTrackBuilder& ttrkbuilder){
  double dca3d2trks_sig = 0;
  double dca2d2trks_sig = 0;
  TransientTrack ttkA = get_ttrk(tkA, ttrkbuilder);
  TransientTrack ttkB = get_ttrk(tkB, ttrkbuilder);
  if(ttkA.impactPointTSCP().isValid() && ttkB.impactPointTSCP().isValid()){
    //Minimum distance
    FreeTrajectoryState state1 = ttkA.impactPointTSCP().theState();
    FreeTrajectoryState state2 = ttkB.impactPointTSCP().theState();
    TwoTrackMinimumDistance minDist;
    minDist.calculate(state1, state2);
    if(minDist.status()){
      //3D distance
      //const float dist3D = minDist.distance();
      std::pair<GlobalPoint,GlobalPoint> pcas = minDist.points();
      GlobalPoint pca1 = pcas.first;
      GlobalPoint pca2 = pcas.second;
      ROOT::Math::SVector<double, 3> distanceVector(pca1.x()-pca2.x(), pca1.y()-pca2.y(), pca1.z()-pca2.z());
      const float twoTkDist3D = ROOT::Math::Mag(distanceVector);
      //3D err distance
      float mass     = 0.139526;
      float massigma = mass*1e-6;
      float chi2 = 0.0f, ndf = 0.0f;
      KinematicParticleFactoryFromTransientTrack pFactory;
      RefCountedKinematicParticle tkAParticle = pFactory.particle(ttkA, mass, chi2, ndf, massigma);
      RefCountedKinematicParticle tkBParticle = pFactory.particle(ttkB, mass, chi2, ndf, massigma);
      float sig[6];
      sig[0] = tkAParticle->stateAtPoint(pca1).kinematicParametersError().matrix()(0,0);
      sig[1] = tkAParticle->stateAtPoint(pca1).kinematicParametersError().matrix()(0,1);
      sig[2] = tkAParticle->stateAtPoint(pca1).kinematicParametersError().matrix()(1,1);
      sig[3] = tkAParticle->stateAtPoint(pca1).kinematicParametersError().matrix()(0,2);
      sig[4] = tkAParticle->stateAtPoint(pca1).kinematicParametersError().matrix()(1,2);
      sig[5] = tkAParticle->stateAtPoint(pca1).kinematicParametersError().matrix()(2,2);
      ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > pca1Cov(sig, sig+6);
      sig[0] = tkBParticle->stateAtPoint(pca2).kinematicParametersError().matrix()(0,0);
      sig[1] = tkBParticle->stateAtPoint(pca2).kinematicParametersError().matrix()(0,1);
      sig[2] = tkBParticle->stateAtPoint(pca2).kinematicParametersError().matrix()(1,1);
      sig[3] = tkBParticle->stateAtPoint(pca2).kinematicParametersError().matrix()(0,2);
      sig[4] = tkBParticle->stateAtPoint(pca2).kinematicParametersError().matrix()(1,2);
      sig[5] = tkBParticle->stateAtPoint(pca2).kinematicParametersError().matrix()(2,2);
      ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > pca2Cov(sig, sig+6);
      ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > totCov = pca1Cov + pca2Cov;
      double twoTkDist3DEr = TMath::Sqrt(fabs(ROOT::Math::Similarity(totCov, distanceVector))) / twoTkDist3D;
      if(twoTkDist3DEr!=0) dca3d2trks_sig = twoTkDist3D/twoTkDist3DEr;
      else                 dca3d2trks_sig = twoTkDist3D/sqrt(twoTkDist3D);
      //2D distance and err distance
      distanceVector(2) = 0.0;
      double twoTkDist2D   = ROOT::Math::Mag(distanceVector);
      double twoTkDist2DEr = TMath::Sqrt(fabs(ROOT::Math::Similarity(totCov, distanceVector))) / twoTkDist2D;
      if(twoTkDist2DEr!=0) dca2d2trks_sig = twoTkDist2D/twoTkDist2DEr;
      else                 dca2d2trks_sig = twoTkDist2D/sqrt(twoTkDist2D); 
    }//if(minDist.status())
  }//ttkA.impactPointTSCP().isValid() && ttkB.impactPointTSCP().isValid()
  return make_pair(dca3d2trks_sig,dca2d2trks_sig);
}
void TauJetnessSelector::get_avip3d(vector<Track> trks, const TransientTrackBuilder& ttrkbuilder, reco::Vertex vtx, vector<tuple<double, double, double> >& tausdir, double& tauchtrks_avip3d_val, double& tauchtrks_avip3d_sig, double& tauchtrks_avsip3d_val, double& tauchtrks_avsip3d_sig){
  double valtemp = 0;
  vector<TransientTrack> ttrks = get_ttrks(trks,ttrkbuilder);
  for(uint t=0; t<ttrks.size(); t++){
    valtemp = IPTools::absoluteImpactParameter3D(ttrks[t],vtx).second.value();
    if(valtemp==valtemp) tauchtrks_avip3d_val  += valtemp;
    valtemp = IPTools::absoluteImpactParameter3D(ttrks[t],vtx).second.significance();
    if(valtemp==valtemp) tauchtrks_avip3d_sig  += valtemp;
    GlobalVector tausdirgv(get<0>(tausdir[t]),get<1>(tausdir[t]),get<2>(tausdir[t]));
    valtemp = IPTools::signedImpactParameter3D(ttrks[t],tausdirgv,vtx).second.value();
    if(valtemp==valtemp) tauchtrks_avsip3d_val += valtemp;
    valtemp = IPTools::signedImpactParameter3D(ttrks[t],tausdirgv,vtx).second.significance();
    if(valtemp==valtemp) tauchtrks_avsip3d_sig += valtemp;
  }
}
void TauJetnessSelector::get_avip2d(vector<Track> trks, const TransientTrackBuilder& ttrkbuilder, reco::Vertex vtx, vector<tuple<double, double, double> >& tausdir, double& tauchtrks_avip2d_val, double& tauchtrks_avip2d_sig, double& tauchtrks_avsip2d_val, double& tauchtrks_avsip2d_sig){
  double valtemp = 0;
  vector<TransientTrack> ttrks = get_ttrks(trks,ttrkbuilder);
  for(uint t=0; t<ttrks.size(); t++){
    valtemp = IPTools::absoluteTransverseImpactParameter(ttrks[t],vtx).second.value();
    if(valtemp==valtemp) tauchtrks_avip2d_val  += valtemp;
    valtemp = IPTools::absoluteTransverseImpactParameter(ttrks[t],vtx).second.significance();
    if(valtemp==valtemp) tauchtrks_avip2d_sig  += valtemp;
    GlobalVector tausdirgv(get<0>(tausdir[t]),get<1>(tausdir[t]),get<2>(tausdir[t]));
    valtemp = IPTools::signedTransverseImpactParameter(ttrks[t],tausdirgv,vtx).second.value();
    if(valtemp==valtemp) tauchtrks_avsip2d_val += valtemp;
    valtemp = IPTools::signedTransverseImpactParameter(ttrks[t],tausdirgv,vtx).second.significance();
    if(valtemp==valtemp) tauchtrks_avsip2d_sig += valtemp;
  }
}
void TauJetnessSelector::get_avip1d(vector<Track> trks, const TransientTrackBuilder& ttrkbuilder, reco::Vertex vtx, vector<tuple<double, double, double> >& tausdir, double& tauchtrks_avip1d_val, double& tauchtrks_avip1d_sig, double& tauchtrks_avsip1d_val, double& tauchtrks_avsip1d_sig){
  double valtemp = 0;
  vector<TransientTrack> ttrks = get_ttrks(trks,ttrkbuilder);
  SignedTransverseImpactParameter stip;
  for(uint t=0; t<ttrks.size(); t++){
    GlobalVector tausdirgv(get<0>(tausdir[t]),get<1>(tausdir[t]),get<2>(tausdir[t]));
    valtemp = fabs(stip.zImpactParameter(ttrks[t],tausdirgv,vtx).second.value());
    if(valtemp==valtemp) tauchtrks_avip1d_val  += valtemp;
    valtemp = fabs(stip.zImpactParameter(ttrks[t],tausdirgv,vtx).second.significance());
    if(valtemp==valtemp) tauchtrks_avip1d_sig  += valtemp;
    valtemp = stip.zImpactParameter(ttrks[t],tausdirgv,vtx).second.value();
    if(valtemp==valtemp) tauchtrks_avsip1d_val += valtemp;
    valtemp = stip.zImpactParameter(ttrks[t],tausdirgv,vtx).second.significance();
    if(valtemp==valtemp) tauchtrks_avsip1d_sig += valtemp; 
  }
}
//To do
//- Think if categories must be defined considering taus ordered by decreasing csv values
//- Include dR between loose muons and loose electrons
//- Include JER for tau selection
//- Check if electron isolation is included in the ID and how it is defined (as it is now or sa it changed for the muon) 
//- Having the vtx in this function is kept for historical reason 
//  and it will be kept until we do not have a final electron loose selection
//- Currently take events where taucsv!=taucsv, but need to understand what to do with them
//- Make sure the tau selection is exactly the same of the TTHbb analysis
//  Same PV, same tau cleaning from mu,ele
//- At the moment valtemp are not included in the average quantities, but the denominator is still the same
//  it should not change too much if the den is high, but check a better implementation
//- Pass maxtaunum from python according to the analysis categorisation
//- For the dca you are considering only the significance.
//  You may want to include the value as well, but probably after data/MC validation
//- Need to decide what to do precisely with 
//  else                       TauJetness_npvTrkOVpvTrk.push_back(-997);
//  If npv!=0 this is likely a btau evt. Check also the pT analogous variables
//- Fare l'average IP solo per le tracce secondarie o primarie?
//  Anche se l'effetto della tauness si dovrebbe vedere con tutte le tracce
