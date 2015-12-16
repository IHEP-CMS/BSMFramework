#include "BSMFramework/BSM3G_TNT_Maker/interface/BJetnessSelector.h"
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
KalmanVertexFitter vtxFitter(true);
using namespace std;
using namespace pat;
using namespace edm;
using namespace reco;
BJetnessSelector::BJetnessSelector(std::string name, TTree* tree, bool debug, const pset& iConfig, edm::ConsumesCollector && ic):
  baseTree(name,tree,debug),
  electronLooseIdMapToken_(ic.consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electronLooseIdMap")))
  {
    _is_data = iConfig.getParameter<bool>("is_data");
  _muonToken        = iConfig.getParameter<edm::InputTag>("muons");
  _patElectronToken = iConfig.getParameter<edm::InputTag>("patElectrons");
  jetToken_         = iConfig.getParameter<edm::InputTag>("jets");
  _vertexInputTag   = iConfig.getParameter<edm::InputTag>("vertices");
  SetBranches();
}
BJetnessSelector::~BJetnessSelector(){
  delete tree_;
}
void BJetnessSelector::Fill(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  Clear();
  /////
  //   Recall collections
  ///// 
  edm::Handle<reco::VertexCollection> vtx_h;
  iEvent.getByLabel(_vertexInputTag, vtx_h);
  edm::Handle<edm::View<pat::Muon> > muon_h;
  iEvent.getByLabel(_muonToken, muon_h);
  edm::Handle<edm::View<pat::Electron> > electron_pat;
  iEvent.getByLabel(_patElectronToken, electron_pat);
  edm::Handle<edm::ValueMap<bool> > loose_id_decisions;
  iEvent.getByToken(electronLooseIdMapToken_,loose_id_decisions);
  edm::Handle<pat::JetCollection> jets;                                       
  iEvent.getByLabel(jetToken_, jets);                                         
  edm::ESHandle<TransientTrackBuilder> ttrkbuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",ttrkbuilder);
  /////
  //   Require a good vertex 
  ///// 
  if(vtx_h->empty()) return; // skip the event if no PV found
  const reco::Vertex &PV = vtx_h->front();
  /////
  //   Look for loose muons, electrons to clean jets
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
    bool isPassLoose  = (*loose_id_decisions) [ elPtr ];
    if(!(is_loose_electron(*ele) && isPassLoose)) continue;
    //bool matchelemu = false;
    //for(uint gl=0; gl<looseleps.size(); gl++) if(deltaR(looseleps[gl]->p4(),ele.p4())<0.3) matchelemu = true;
    //if(matchelemu) continue;
    const pat::Electron &lele = *ele;
    looseleps.push_back((const reco::Candidate*)&lele);
  }
  /////
  //   Get the good jets of the event
  /////  
  //Iterate to access jet by decreasing b-tagging value
  int jet_pos = 0; //This counter helps to order jets
  int jet_num = 0; //This counter accounts for the number of good jets in the events
                   //The definition of good jet in the event must be the same of the TTHbb analysis
                   //so that jet_num corresponds to the number of jets that define the categories in the TTHbb search
  vector<pair<double,int> > jet_csv_pos;
  for(const pat::Jet &j : *jets){ 
    //Minimum jet selections
    if(!is_good_jet(j)){jet_pos++; continue;}
    bool jetmatchedlepts = false;
    for(uint gl=0; gl<looseleps.size(); gl++) if(deltaR(looseleps[gl]->p4(),j.p4())<0.4) jetmatchedlepts = true;
    if(jetmatchedlepts){jet_pos++; continue;}
    double csvcurrjet = j.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
    jet_csv_pos.push_back(make_pair(csvcurrjet,jet_pos));
    jet_pos++;
    jet_num++;
  }
  sort(jet_csv_pos.rbegin(), jet_csv_pos.rend());
  if(jet_num!=0){
    //cout<<"Num of jet is"<<setw(20)<<savebjetnessevt<<endl;
    BJetness_numjet.push_back(jet_num);
    for(int jn=0; jn<jet_num; jn++){
      const pat::Jet & j = (*jets)[jet_csv_pos[jn].second];
      //Kinematics and csv 
      BJetness_jetpt.push_back(j.pt());
      BJetness_jeteta.push_back(j.eta());
      BJetness_jetphi.push_back(j.phi());
      BJetness_jetenergy.push_back(j.energy());
      double jetcsv = j.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"); 
      if(jetcsv!=jetcsv) jetcsv = -996;
      BJetness_jetcsv.push_back(jetcsv);
      double jetbjetprob = j.bDiscriminator("pfJetProbabilityBJetTags");
      if(jetbjetprob!=jetbjetprob) jetbjetprob = -996; 
      BJetness_pfJetProbabilityBJetTags.push_back(jetbjetprob);
      double jetbjetcombmva = j.bDiscriminator("pfCombinedMVABJetTags");
      if(jetbjetcombmva!=jetbjetcombmva) jetbjetcombmva = -996;
      BJetness_pfCombinedMVABJetTags.push_back(jetbjetcombmva);
      //cout<<"Jet pf"<<setw(20)<<"hf"<<setw(20)<<"pt"<<setw(20)<<"eta"<<setw(20)<<"phi"<<setw(20)<<"energy"<<setw(20)<<"csv"<<endl;
      //cout<<j.partonFlavour()<<setw(20)<<j.hadronFlavour()<<setw(20)<<j.pt()<<setw(20)<<j.eta()<<setw(20)<<j.phi()<<setw(20)<<j.energy()<<setw(20)<<j.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags")<<endl;
      //cout<<"CSV "<<setw(20)<<j<<setw(20)<<j.bDiscriminator("pfJetProbabilityBJetTags")<<setw(20)<<j.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags")<<setw(20)<<j.bDiscriminator("pfCombinedMVABJetTags")<<endl;
      BJetness_partonFlavour.push_back(j.partonFlavour());
      BJetness_hadronFlavour.push_back(j.hadronFlavour());
    }
    /////
    //   Gen info
    /////
    if(!_is_data) {
      BJetness_ngenbh = 0;
      BJetness_ngenbt = 0;
      BJetness_ngenb  = 0;
      BJetness_ngenc  = 0;
      Handle< reco::GenParticleCollection > pruned;
      iEvent.getByLabel("prunedGenParticles", pruned);
      for(size_t i=0; i<pruned->size(); i++){
	if(abs((*pruned)[i].pdgId())==5){// && (*pruned)[i].pt()>15)
	  const Candidate * bquark = &(*pruned)[i];
	  if(bquark->mother(0)!=0 && abs(bquark->mother(0)->pdgId())==25)BJetness_ngenbh++;
	  if(bquark->mother(0)!=0 && abs(bquark->mother(0)->pdgId())==6) BJetness_ngenbt++;
	  if(bquark->mother(0)!=0 && abs(bquark->mother(0)->pdgId())!=5) BJetness_ngenb++; 
	}
	if(abs((*pruned)[i].pdgId())==4){// && (*pruned)[i].pt()>15)
	  const Candidate * cquark = &(*pruned)[i];
	  if(cquark->mother(0)!=0 && abs(cquark->mother(0)->pdgId())!=4) BJetness_ngenc++; 
	}   
      }
      //cout<<setw(20)<<"BJetness_ngenbh"<<setw(20)<<"BJetness_ngenbt"<<setw(20)<<"BJetness_ngenb"<<setw(20)<<"BJetness_ngenc"<<endl;
      //cout<<setw(20)<<BJetness_ngenbh<<setw(20)<<BJetness_ngenbt<<setw(20)<<BJetness_ngenb<<setw(20)<<BJetness_ngenc<<endl;
    }
    //Access info jet by jet 
    vector<Track> jetschtrks; jetschtrks.clear();
    vector<Track> jetschtrkspv; jetschtrkspv.clear();
    vector<Track> jetschtrksnpv; jetschtrksnpv.clear();
    vector<tuple<double, double, double> > jetsdir; jetsdir.clear();
    double bjetness_num_pdgid_eles = 0;
    double bjetness_num_pdgid_mus  = 0;
    double bjetness_num_soft_eles  = 0;
    double bjetness_num_vetonoipnoiso_eles  = 0;
    double bjetness_num_loosenoipnoiso_eles = 0;
    double bjetness_num_loose_mus  = 0;
    int maxjetnum = 6; //Pass this value from python at some point
                       //The values will have to be chosen depending on the definition of the categories
    if(jet_num<maxjetnum)  maxjetnum = jet_num;
    for(int jn=1; jn<maxjetnum; jn++){
      const pat::Jet & j = (*jets)[jet_csv_pos[jn].second];
      get_jettrks(j, PV, *ttrkbuilder, jetschtrks, jetschtrkspv, jetschtrksnpv, jetsdir,
                  electron_pat, muon_h,
                  bjetness_num_pdgid_eles, bjetness_num_pdgid_mus, bjetness_num_soft_eles, bjetness_num_vetonoipnoiso_eles, bjetness_num_loosenoipnoiso_eles, bjetness_num_loose_mus 
                 );
    }  
    /////
    //   Get info to evaluate the event BJetness
    /////
    //Num_of_trks
    //cout<<setw(20)<<"numlep"<<setw(20)<<bjetness_num_pdgid_eles+bjetness_num_pdgid_mus<<setw(20)<<bjetness_num_pdgid_eles<<setw(20)<<bjetness_num_pdgid_mus<<setw(20)<<bjetness_num_soft_eles+bjetness_num_loose_mus<<setw(20)<<bjetness_num_soft_eles<<setw(20)<<bjetness_num_vetonoipnoiso_eles+bjetness_num_loose_mus<<setw(20)<<bjetness_num_vetonoipnoiso_eles<<setw(20)<<bjetness_num_loosenoipnoiso_eles+bjetness_num_loose_mus<<setw(20)<<bjetness_num_loosenoipnoiso_eles<<setw(20)<<bjetness_num_loose_mus<<endl;
    BJetness_num_pdgid_leps.push_back(bjetness_num_pdgid_eles+bjetness_num_pdgid_mus);
    BJetness_num_pdgid_eles.push_back(bjetness_num_pdgid_eles);
    BJetness_num_pdgid_mus.push_back(bjetness_num_pdgid_mus);
    BJetness_num_soft_leps.push_back(bjetness_num_soft_eles+bjetness_num_loose_mus);
    BJetness_num_soft_eles.push_back(bjetness_num_soft_eles);
    BJetness_num_vetonoipnoiso_leps.push_back(bjetness_num_vetonoipnoiso_eles+bjetness_num_loose_mus);
    BJetness_num_vetonoipnoiso_eles.push_back(bjetness_num_vetonoipnoiso_eles);
    BJetness_num_loosenoipnoiso_leps.push_back(bjetness_num_loosenoipnoiso_eles+bjetness_num_loose_mus);
    BJetness_num_loosenoipnoiso_eles.push_back(bjetness_num_loosenoipnoiso_eles);
    BJetness_num_loose_mus.push_back(bjetness_num_loose_mus);
    BJetness_numjettrks.push_back(jetschtrks.size());
    BJetness_numjettrkspv.push_back(jetschtrkspv.size());
    BJetness_numjettrksnopv.push_back(jetschtrksnpv.size());
    //cout<<setw(20)<<"Num_of_trks"<<setw(20)<<jetschtrks.size()<<setw(20)<<jetschtrkspv.size()<<setw(20)<<jetschtrksnpv.size()<<endl;
    if(jetschtrks.size()!=0){
      //cout<<setw(20)<<"nonPVTrk/CollTrk"<<setw(20)<<"PVTrk/CollTrk"<<setw(20)<<"nonPVTrk/PVTrk"<<endl;
      //cout<<double(jetschtrksnpv.size())/double(jetschtrks.size())<<setw(20)<<double(jetschtrkspv.size())/double(jetschtrks.size())<<setw(20)<<double(jetschtrksnpv.size())/double(jetschtrkspv.size())<<endl;
      BJetness_npvTrkOVcollTrk.push_back(double(jetschtrksnpv.size())/double(jetschtrks.size()));
      BJetness_pvTrkOVcollTrk.push_back(double(jetschtrkspv.size())/double(jetschtrks.size()));
      if(jetschtrkspv.size()!=0) BJetness_npvTrkOVpvTrk.push_back(double(jetschtrksnpv.size())/double(jetschtrkspv.size()));
      else                       BJetness_npvTrkOVpvTrk.push_back(-997);
      //cout<<setw(20)<<"nonPVPt/CollPT"<<setw(20)<<"PVPt/CollPT"<<setw(20)<<"nonPVPt/PVPt"<<endl;
      double ptjettrks    = 0;
      for(uint jt=0; jt<jetschtrks.size(); jt++) ptjettrks += jetschtrks[jt].pt();
      double ptjettrkspv    = 0;
      for(uint jt=0; jt<jetschtrkspv.size(); jt++) ptjettrkspv += jetschtrkspv[jt].pt();
      double ptjettrksnpv    = 0;
      for(uint jt=0; jt<jetschtrksnpv.size(); jt++) ptjettrksnpv += jetschtrksnpv[jt].pt();
      //cout<<setw(20)<<ptjettrksnpv/ptjettrks<<setw(20)<<ptjettrkspv/ptjettrks<<setw(20)<<ptjettrksnpv/ptjettrkspv<<endl;
      if(ptjettrks!=0) BJetness_npvPtOVcollPt.push_back(ptjettrksnpv/ptjettrks);
      else             BJetness_npvPtOVcollPt.push_back(-997);
      //cout<<"BJetness_npvTrkOVcollTrk "<<double(jetschtrksnpv.size())/double(jetschtrks.size())<<" "<<double(jetschtrksnpv.size())<<" "<<double(jetschtrks.size())<<endl;
      //cout<<"BJetness_npvPtOVcollPt "<<ptjettrksnpv/ptjettrks<<" "<<ptjettrksnpv<<" "<<ptjettrks<<endl;
      if(ptjettrks!=0) BJetness_pvPtOVcollPt.push_back(ptjettrkspv/ptjettrks);
      else             BJetness_pvPtOVcollPt.push_back(-997);
      if(ptjettrkspv!=0) BJetness_npvPtOVpvPt.push_back(ptjettrksnpv/ptjettrkspv);    
      else               BJetness_npvPtOVpvPt.push_back(-997);
      //Two_trk_info (we may want to use one function per variable. Probably more clear, but slower!)
      double jetchtrks_num2v     = 0;
      double jetchtrks_numno2v   = 0;
      double jetchtrks_dca3d2t   = 0;
      double jetchtrks_dca3dno2t = 0;
      double jetchtrks_dca2d2t   = 0;
      double jetchtrks_dca2dno2t = 0;
      get_2trksinfo(jetschtrks, *ttrkbuilder, jetchtrks_num2v, jetchtrks_numno2v, jetchtrks_dca3d2t, jetchtrks_dca3dno2t, jetchtrks_dca2d2t, jetchtrks_dca2dno2t);
      if((jetchtrks_num2v+jetchtrks_numno2v)!=0){
        double twotrksinfo_valtemp = 0;
        twotrksinfo_valtemp = jetchtrks_num2v/(jetchtrks_num2v+jetchtrks_numno2v);
        if(twotrksinfo_valtemp==twotrksinfo_valtemp) BJetness_avnum2v.push_back(twotrksinfo_valtemp);
        else                                         BJetness_avnum2v.push_back(-996);
        twotrksinfo_valtemp = jetchtrks_numno2v/(jetchtrks_num2v+jetchtrks_numno2v);
        if(twotrksinfo_valtemp==twotrksinfo_valtemp) BJetness_avnumno2v.push_back(twotrksinfo_valtemp);
        else                                         BJetness_avnumno2v.push_back(-996);
        twotrksinfo_valtemp = jetchtrks_dca3d2t/(jetchtrks_num2v+jetchtrks_numno2v);
        if(twotrksinfo_valtemp==twotrksinfo_valtemp) BJetness_avdca3d2t.push_back(twotrksinfo_valtemp);
        else                                         BJetness_avdca3d2t.push_back(-996);
        twotrksinfo_valtemp = jetchtrks_dca3dno2t/(jetchtrks_num2v+jetchtrks_numno2v);
        if(twotrksinfo_valtemp==twotrksinfo_valtemp) BJetness_avdca3dno2t.push_back(twotrksinfo_valtemp);
        else                                         BJetness_avdca3dno2t.push_back(-996);
        twotrksinfo_valtemp = (jetchtrks_dca3d2t+jetchtrks_dca3dno2t)/(jetchtrks_num2v+jetchtrks_numno2v);
        if(twotrksinfo_valtemp==twotrksinfo_valtemp) BJetness_avdca3d.push_back(twotrksinfo_valtemp);
        else                                         BJetness_avdca3d.push_back(-996);
        twotrksinfo_valtemp = jetchtrks_dca2d2t/(jetchtrks_num2v+jetchtrks_numno2v);
        if(twotrksinfo_valtemp==twotrksinfo_valtemp) BJetness_avdca2d2t.push_back(twotrksinfo_valtemp);
        else                                         BJetness_avdca2d2t.push_back(-996);
        twotrksinfo_valtemp = jetchtrks_dca2dno2t/(jetchtrks_num2v+jetchtrks_numno2v);
        if(twotrksinfo_valtemp==twotrksinfo_valtemp) BJetness_avdca2dno2t.push_back(twotrksinfo_valtemp);
        else                                         BJetness_avdca2dno2t.push_back(-996);
        twotrksinfo_valtemp = (jetchtrks_dca2d2t+jetchtrks_dca2dno2t)/(jetchtrks_num2v+jetchtrks_numno2v);
        if(twotrksinfo_valtemp==twotrksinfo_valtemp) BJetness_avdca2d.push_back(twotrksinfo_valtemp);
        else                                         BJetness_avdca2d.push_back(-996);
      }else{
        BJetness_avnum2v.push_back(-997);
        BJetness_avnumno2v.push_back(-997);
        BJetness_avdca3d2t.push_back(-997);
        BJetness_avdca3dno2t.push_back(-997);
        BJetness_avdca3d.push_back(-997);
        BJetness_avdca2d2t.push_back(-997);
        BJetness_avdca2dno2t.push_back(-997);
        BJetness_avdca2d.push_back(-997);
      }
      //cout<<setw(20)<<"Two_trk_info"<<setw(20)<<(jetchtrks_num2v+jetchtrks_numno2v)<<endl;
      //chi2
      double jetchtrks_chi2 = 997;
      get_chi2(jetschtrks, *ttrkbuilder, jetchtrks_chi2);
      BJetness_chi2.push_back(jetchtrks_chi2);
      double ip_valtemp = 0;
      //ImpactParameter  
      double jetchtrks_avip3d_val  = 0;  
      double jetchtrks_avip3d_sig  = 0;  
      double jetchtrks_avsip3d_val = 0;  
      double jetchtrks_avsip3d_sig = 0;  
      get_avip3d(jetschtrks, *ttrkbuilder, PV, jetsdir, jetchtrks_avip3d_val,jetchtrks_avip3d_sig,jetchtrks_avsip3d_val,jetchtrks_avsip3d_sig);
      ip_valtemp = jetchtrks_avip3d_val/jetschtrks.size();
      if(ip_valtemp==ip_valtemp) BJetness_avip3d_val.push_back(ip_valtemp);
      else                       BJetness_avip3d_val.push_back(-996);
      ip_valtemp = jetchtrks_avip3d_sig/jetschtrks.size();
      if(ip_valtemp==ip_valtemp) BJetness_avip3d_sig.push_back(ip_valtemp);
      else                       BJetness_avip3d_sig.push_back(-996); 
      ip_valtemp = jetchtrks_avsip3d_val/jetschtrks.size();
      if(ip_valtemp==ip_valtemp) BJetness_avsip3d_val.push_back(ip_valtemp);
      else                       BJetness_avsip3d_val.push_back(-996); 
      ip_valtemp = jetchtrks_avsip3d_sig/jetschtrks.size();
      if(ip_valtemp==ip_valtemp) BJetness_avsip3d_sig.push_back(ip_valtemp);
      else                       BJetness_avsip3d_sig.push_back(-996);
      //cout<<setw(20)<<"ImpactParameter"<<setw(20)<<jetchtrks_avip3d_val<<setw(20)<<jetchtrks_avip3d_sig<<setw(20)<<jetchtrks_avsip3d_val<<setw(20)<<jetchtrks_avsip3d_sig<<endl;
      double jetchtrks_avip2d_val  = 0;  
      double jetchtrks_avip2d_sig  = 0;  
      double jetchtrks_avsip2d_val = 0;  
      double jetchtrks_avsip2d_sig = 0;  
      get_avip2d(jetschtrks, *ttrkbuilder, PV, jetsdir, jetchtrks_avip2d_val,jetchtrks_avip2d_sig,jetchtrks_avsip2d_val,jetchtrks_avsip2d_sig);
      ip_valtemp = jetchtrks_avip2d_val/jetschtrks.size();
      if(ip_valtemp==ip_valtemp) BJetness_avip2d_val.push_back(ip_valtemp);
      else                       BJetness_avip2d_val.push_back(-996);
      ip_valtemp = jetchtrks_avip2d_sig/jetschtrks.size();
      if(ip_valtemp==ip_valtemp) BJetness_avip2d_sig.push_back(ip_valtemp);
      else                       BJetness_avip2d_sig.push_back(-996);
      ip_valtemp = jetchtrks_avsip2d_val/jetschtrks.size();
      if(ip_valtemp==ip_valtemp) BJetness_avsip2d_val.push_back(ip_valtemp);
      else                       BJetness_avsip2d_val.push_back(-996);
      ip_valtemp = jetchtrks_avsip2d_sig/jetschtrks.size();
      if(ip_valtemp==ip_valtemp) BJetness_avsip2d_sig.push_back(ip_valtemp);
      else                       BJetness_avsip2d_sig.push_back(-996);
      double jetchtrks_avip1d_val  = 0;  
      double jetchtrks_avip1d_sig  = 0;  
      double jetchtrks_avsip1d_val = 0;  
      double jetchtrks_avsip1d_sig = 0;  
      get_avip1d(jetschtrks, *ttrkbuilder, PV, jetsdir, jetchtrks_avip1d_val,jetchtrks_avip1d_sig,jetchtrks_avsip1d_val,jetchtrks_avsip1d_sig);
      ip_valtemp = jetchtrks_avip1d_val/jetschtrks.size();
      if(ip_valtemp==ip_valtemp) BJetness_avip1d_val.push_back(ip_valtemp);
      else                       BJetness_avip1d_val.push_back(-996);
      ip_valtemp = jetchtrks_avip1d_sig/jetschtrks.size();
      if(ip_valtemp==ip_valtemp) BJetness_avip1d_sig.push_back(ip_valtemp);
      else                       BJetness_avip1d_sig.push_back(-996);
      ip_valtemp = jetchtrks_avsip1d_val/jetschtrks.size();
      if(ip_valtemp==ip_valtemp) BJetness_avsip1d_val.push_back(ip_valtemp);
      else                       BJetness_avsip1d_val.push_back(-996);
      ip_valtemp = jetchtrks_avsip1d_sig/jetschtrks.size();
      if(ip_valtemp==ip_valtemp) BJetness_avsip1d_sig.push_back(ip_valtemp);
      else                       BJetness_avsip1d_sig.push_back(-996);
    }else{//if(jetschtrks.size()!=0)
      //Num trk pt
      BJetness_npvTrkOVcollTrk.push_back(-998);
      BJetness_pvTrkOVcollTrk.push_back(-998);
      BJetness_npvTrkOVpvTrk.push_back(-998);
      BJetness_npvPtOVcollPt.push_back(-998);
      BJetness_pvPtOVcollPt.push_back(-998);
      BJetness_npvPtOVpvPt.push_back(-998);
      //Two_trk_info (we may want to use one function per variable. Probably more clear, but slower!)
      BJetness_avnum2v.push_back(-998);
      BJetness_avnumno2v.push_back(-998);
      BJetness_avdca3d2t.push_back(-998);
      BJetness_avdca3dno2t.push_back(-998);
      BJetness_avdca3d.push_back(-998);
      BJetness_avdca2d2t.push_back(-998);
      BJetness_avdca2dno2t.push_back(-998);
      BJetness_avdca2d.push_back(-998);
      //chi2
      BJetness_chi2.push_back(-998);
      //ImpactParameter  
      BJetness_avip3d_val.push_back(-998);
      BJetness_avip3d_sig.push_back(-998);
      BJetness_avsip3d_val.push_back(-998);
      BJetness_avsip3d_sig.push_back(-998);
      BJetness_avip2d_val.push_back(-998);
      BJetness_avip2d_sig.push_back(-998);
      BJetness_avsip2d_val.push_back(-998);
      BJetness_avsip2d_sig.push_back(-998);
      BJetness_avip1d_val.push_back(-998);
      BJetness_avip1d_sig.push_back(-998);
      BJetness_avsip1d_val.push_back(-998);
      BJetness_avsip1d_sig.push_back(-998);
    }
 }else{//if(jet_num!=0)
   BJetness_partonFlavour.push_back(-999);
   BJetness_hadronFlavour.push_back(-999);
   BJetness_numjet.push_back(-999);
   //Kinematics and csv 
   BJetness_jetpt.push_back(-999);
   BJetness_jeteta.push_back(-999);
   BJetness_jetphi.push_back(-999);
   BJetness_jetenergy.push_back(-999);
   BJetness_jetcsv.push_back(-999);
   BJetness_pfJetProbabilityBJetTags.push_back(-999);
   BJetness_pfCombinedMVABJetTags.push_back(-999);
   //Num_of_trks
   BJetness_num_pdgid_leps.push_back(-999);
   BJetness_num_pdgid_eles.push_back(-999);
   BJetness_num_pdgid_mus.push_back(-999);
   BJetness_num_soft_leps.push_back(-999);
   BJetness_num_soft_eles.push_back(-999);
   BJetness_num_vetonoipnoiso_leps.push_back(-999);
   BJetness_num_vetonoipnoiso_eles.push_back(-999);
   BJetness_num_loosenoipnoiso_leps.push_back(-999);
   BJetness_num_loosenoipnoiso_eles.push_back(-999);
   BJetness_num_loose_mus.push_back(-999);
   BJetness_numjettrks.push_back(-999);
   BJetness_numjettrkspv.push_back(-999);
   BJetness_numjettrksnopv.push_back(-999);
   BJetness_npvTrkOVcollTrk.push_back(-999);
   BJetness_pvTrkOVcollTrk.push_back(-999);
   BJetness_npvTrkOVpvTrk.push_back(-999);
   BJetness_npvPtOVcollPt.push_back(-999);
   BJetness_pvPtOVcollPt.push_back(-999);
   BJetness_npvPtOVpvPt.push_back(-999);
   //Two_trk_info (we may want to use one function per variable. Probably more clear, but slower!)
   BJetness_avnum2v.push_back(-999);
   BJetness_avnumno2v.push_back(-999);
   BJetness_avdca3d2t.push_back(-999);
   BJetness_avdca3dno2t.push_back(-999);
   BJetness_avdca3d.push_back(-999);
   BJetness_avdca2d2t.push_back(-999);
   BJetness_avdca2dno2t.push_back(-999);
   BJetness_avdca2d.push_back(-999);
   //chi2
   BJetness_chi2.push_back(-999);
   //ImpactParameter  
   BJetness_avip3d_val.push_back(-999);
   BJetness_avip3d_sig.push_back(-999);
   BJetness_avsip3d_val.push_back(-999);
   BJetness_avsip3d_sig.push_back(-999);
   BJetness_avip2d_val.push_back(-999);
   BJetness_avip2d_sig.push_back(-999);
   BJetness_avsip2d_val.push_back(-999);
   BJetness_avsip2d_sig.push_back(-999);
   BJetness_avip1d_val.push_back(-999);
   BJetness_avip1d_sig.push_back(-999);
   BJetness_avsip1d_val.push_back(-999);
   BJetness_avsip1d_sig.push_back(-999);
 }
}
void BJetnessSelector::SetBranches(){
  if(debug_) std::cout<<"setting branches: calling AddBranch of baseTree"<<std::endl;
  //Gen info
  AddBranch(&BJetness_ngenbh                   ,"BJetness_ngenbh");
  AddBranch(&BJetness_ngenbt                   ,"BJetness_ngenbt");
  AddBranch(&BJetness_ngenb                    ,"BJetness_ngenb");
  AddBranch(&BJetness_ngenc                    ,"BJetness_ngenc");
  AddBranch(&BJetness_partonFlavour            ,"BJetness_partonFlavour");
  AddBranch(&BJetness_hadronFlavour            ,"BJetness_hadronFlavour");
  AddBranch(&BJetness_numjet                   ,"BJetness_numjet");
  //Kinematics and csv 
  AddBranch(&BJetness_jetpt                    ,"BJetness_jetpt");
  AddBranch(&BJetness_jeteta                   ,"BJetness_jeteta");
  AddBranch(&BJetness_jetphi                   ,"BJetness_jetphi");
  AddBranch(&BJetness_jetenergy                ,"BJetness_jetenergy");
  AddBranch(&BJetness_jetcsv                   ,"BJetness_jetcsv");
  AddBranch(&BJetness_pfJetProbabilityBJetTags ,"BJetness_pfJetProbabilityBJetTags");
  AddBranch(&BJetness_pfCombinedMVABJetTags    ,"BJetness_pfCombinedMVABJetTags");
  //Num_of_trks
  AddBranch(&BJetness_num_pdgid_leps           ,"BJetness_num_pdgid_leps");
  AddBranch(&BJetness_num_pdgid_eles           ,"BJetness_num_pdgid_eles");
  AddBranch(&BJetness_num_pdgid_mus            ,"BJetness_num_pdgid_mus");
  AddBranch(&BJetness_num_soft_leps            ,"BJetness_num_soft_leps");
  AddBranch(&BJetness_num_soft_eles            ,"BJetness_num_soft_eles");
  AddBranch(&BJetness_num_vetonoipnoiso_leps   ,"BJetness_num_vetonoipnoiso_leps");
  AddBranch(&BJetness_num_vetonoipnoiso_eles   ,"BJetness_num_vetonoipnoiso_eles");
  AddBranch(&BJetness_num_loosenoipnoiso_leps  ,"BJetness_num_loosenoipnoiso_leps");
  AddBranch(&BJetness_num_loosenoipnoiso_eles  ,"BJetness_num_loosenoipnoiso_eles");
  AddBranch(&BJetness_num_loose_mus            ,"BJetness_num_loose_mus");
  AddBranch(&BJetness_numjettrks               ,"BJetness_numjettrks");
  AddBranch(&BJetness_numjettrkspv             ,"BJetness_numjettrkspv");
  AddBranch(&BJetness_numjettrksnopv           ,"BJetness_numjettrksnopv");
  AddBranch(&BJetness_npvTrkOVcollTrk          ,"BJetness_npvTrkOVcollTrk");
  AddBranch(&BJetness_pvTrkOVcollTrk           ,"BJetness_pvTrkOVcollTrk");
  AddBranch(&BJetness_npvTrkOVpvTrk            ,"BJetness_npvTrkOVpvTrk");
  AddBranch(&BJetness_npvPtOVcollPt            ,"BJetness_npvPtOVcollPt");
  AddBranch(&BJetness_pvPtOVcollPt             ,"BJetness_pvPtOVcollPt");
  AddBranch(&BJetness_npvPtOVpvPt              ,"BJetness_npvPtOVpvPt");
  //Two_trk_info
  AddBranch(&BJetness_avnum2v                  ,"BJetness_avnum2v");
  AddBranch(&BJetness_avnumno2v                ,"BJetness_avnumno2v");
  AddBranch(&BJetness_avdca3d2t                ,"BJetness_avdca3d2t");
  AddBranch(&BJetness_avdca3dno2t              ,"BJetness_avdca3dno2t");
  AddBranch(&BJetness_avdca3d                  ,"BJetness_avdca3d");
  AddBranch(&BJetness_avdca2d2t                ,"BJetness_avdca2d2t");
  AddBranch(&BJetness_avdca2dno2t              ,"BJetness_avdca2dno2t");
  AddBranch(&BJetness_avdca2d                  ,"BJetness_avdca2d");
  //chi2
  AddBranch(&BJetness_chi2                     ,"BJetness_chi2");
  //ImpactParameter
  AddBranch(&BJetness_avip3d_val               ,"BJetness_avip3d_val");
  AddBranch(&BJetness_avip3d_sig               ,"BJetness_avip3d_sig");
  AddBranch(&BJetness_avsip3d_val              ,"BJetness_avsip3d_val");
  AddBranch(&BJetness_avsip3d_sig              ,"BJetness_avsip3d_sig");
  AddBranch(&BJetness_avip2d_val               ,"BJetness_avip2d_val");
  AddBranch(&BJetness_avip2d_sig               ,"BJetness_avip2d_sig");
  AddBranch(&BJetness_avsip2d_val              ,"BJetness_avsip2d_val");
  AddBranch(&BJetness_avsip2d_sig              ,"BJetness_avsip2d_sig");
  AddBranch(&BJetness_avip1d_val               ,"BJetness_avip1d_val");
  AddBranch(&BJetness_avip1d_sig               ,"BJetness_avip1d_sig");
  AddBranch(&BJetness_avsip1d_val              ,"BJetness_avsip1d_val");
  AddBranch(&BJetness_avsip1d_sig              ,"BJetness_avsip1d_sig");
  if(debug_) std::cout<<"set branches"<<std::endl;
}
void BJetnessSelector::Clear(){
  //Gen info
  BJetness_ngenbh = -9999;
  BJetness_ngenbt = -9999;
  BJetness_ngenb  = -9999;
  BJetness_ngenc  = -9999;
  BJetness_partonFlavour.clear();
  BJetness_hadronFlavour.clear();
  BJetness_numjet.clear();
  //Kinematics and csv 
  BJetness_jetpt.clear(); 
  BJetness_jeteta.clear();
  BJetness_jetphi.clear();
  BJetness_jetenergy.clear();
  BJetness_jetcsv.clear(); 
  BJetness_pfJetProbabilityBJetTags.clear();
  BJetness_pfCombinedMVABJetTags.clear();
  //Num_of_trks
  BJetness_num_pdgid_leps.clear();
  BJetness_num_pdgid_eles.clear();
  BJetness_num_pdgid_mus.clear();
  BJetness_num_soft_leps.clear();
  BJetness_num_soft_eles.clear();
  BJetness_num_vetonoipnoiso_leps.clear();
  BJetness_num_vetonoipnoiso_eles.clear();
  BJetness_num_loosenoipnoiso_leps.clear();
  BJetness_num_loosenoipnoiso_eles.clear();
  BJetness_num_loose_mus.clear();
  BJetness_numjettrks.clear();
  BJetness_numjettrkspv.clear();
  BJetness_numjettrksnopv.clear();
  BJetness_npvTrkOVcollTrk.clear();
  BJetness_pvTrkOVcollTrk.clear();
  BJetness_npvTrkOVpvTrk.clear();
  BJetness_npvPtOVcollPt.clear();
  BJetness_pvPtOVcollPt.clear();
  BJetness_npvPtOVpvPt.clear();
  //Two_trk_info
  BJetness_avnum2v.clear();
  BJetness_avnumno2v.clear();
  BJetness_avdca3d2t.clear();
  BJetness_avdca3dno2t.clear();
  BJetness_avdca3d.clear();
  BJetness_avdca2d2t.clear();
  BJetness_avdca2dno2t.clear();
  BJetness_avdca2d.clear();
  //chi2
  BJetness_chi2.clear();
  //ImpactParameter
  BJetness_avip3d_val.clear();
  BJetness_avip3d_sig.clear();
  BJetness_avsip3d_val.clear();
  BJetness_avsip3d_sig.clear();
  BJetness_avip2d_val.clear();
  BJetness_avip2d_sig.clear();
  BJetness_avsip2d_val.clear();
  BJetness_avsip2d_sig.clear();
  BJetness_avip1d_val.clear();
  BJetness_avip1d_sig.clear();
  BJetness_avsip1d_val.clear();
  BJetness_avsip1d_sig.clear();
}
//Look for loose muon
bool BJetnessSelector::is_loose_muon(const pat::Muon& mu, const reco::Vertex& vtx){
  bool isloosemu = false;
  if(mu.muonBestTrack().isNonnull() && mu.globalTrack().isNonnull() && mu.innerTrack().isNonnull() && mu.track().isNonnull() &&
    mu.pt()>10 &&
    TMath::Abs(mu.eta()) < 2.4 &&
    mu.isTightMuon(vtx) &&  
    mu.isPFMuon() && mu.isGlobalMuon() &&
    fabs(mu.muonBestTrack()->dxy(vtx.position())) < 0.2 && 
    fabs(mu.muonBestTrack()->dz(vtx.position()))  < 0.5 && 
    mu.globalTrack()->normalizedChi2() < 10 &&
    mu.globalTrack()->hitPattern().numberOfValidMuonHits()  > 0 &&
    mu.innerTrack()->hitPattern().numberOfValidPixelHits()  > 0 &&
    mu.track()->hitPattern().trackerLayersWithMeasurement() > 5 &&
    mu.numberOfMatchedStations() > 1 && 
    rel_iso_dbc_mu(mu) < 0.2
    ) isloosemu = true;
  return isloosemu;
}
double BJetnessSelector::rel_iso_dbc_mu(const pat::Muon& lepton){
  return((lepton.pfIsolationR04().sumChargedHadronPt + max(lepton.pfIsolationR04().sumNeutralHadronEt + lepton.pfIsolationR04().sumPhotonEt - 0.5 * lepton.pfIsolationR04().sumPUPt,0.0))/lepton.pt()
        );
}
bool BJetnessSelector::is_loose_electron(const pat::Electron& ele){
  bool isele = false;
  if(ele.pt()>10 && TMath::Abs(ele.eta())<2.4 &&
     !(abs(ele.superCluster()->position().eta()) > 1.4442 && abs(ele.superCluster()->position().eta()) < 1.5660) &&
     ele.passConversionVeto()
    ) isele = true;
  return isele; 
}
//Having the vtx in this function is kept for historical reason 
//and it will be kept until we do not have a final electron loose selection
//bool BJetnessSelector::is_loose_electron(const pat::Electron& ele, const reco::Vertex& vtx){
//  bool islooseele = false;
//  if(ele.pt()>10 && TMath::Abs(ele.eta())<2.4 &&
//     !(abs(ele.superCluster()->position().eta()) > 1.4442 && abs(ele.superCluster()->position().eta()) < 1.5660) &&
//     ele.passConversionVeto()
//    ) islooseele = true;
//  return islooseele; 
//}
double BJetnessSelector::rel_iso_dbc_ele(const pat::Electron& lepton){
  return( (lepton.chargedHadronIso() +
          std::max(0.0, lepton.neutralHadronIso() + lepton.photonIso() - 0.5*lepton.puChargedHadronIso()))/lepton.pt() );
}
//Look for loose electron
//Require good jets
//This function has to be updated in order to select good jet using the TTHbb definition
bool BJetnessSelector::is_good_jet(const pat::Jet &j){
  bool isgoodjet = true;
  //Acceptance
  if(j.pt() < 30)       isgoodjet = false;
  if(fabs(j.eta())>2.4) isgoodjet = false; 
  //ID requirements
  if(j.neutralHadronEnergyFraction() >= 0.99) isgoodjet = false;
  if(j.chargedEmEnergyFraction()     >= 0.99) isgoodjet = false;
  if(j.neutralEmEnergyFraction()     >= 0.99) isgoodjet = false;
  if(j.numberOfDaughters()           <= 1)    isgoodjet = false;
  if(j.chargedHadronEnergyFraction() <= 0.0)  isgoodjet = false;
  if(j.chargedMultiplicity()         <= 0.0)  isgoodjet = false;
  return isgoodjet;
}
bool BJetnessSelector::is_loosePOG_jetmuon(const pat::PackedCandidate &jcand, edm::Handle<edm::View<pat::Muon> > muon_h){
  bool ismu = false;
  for(const pat::Muon &mu : *muon_h){
    if(deltaR(jcand.p4(),mu.p4())<0.1 && fabs(jcand.pt()-mu.pt())/mu.pt()<0.05){
     if(mu.isLooseMuon()) ismu = true;
     if(ismu) break;
    }
  }  
  return ismu;
}
bool BJetnessSelector::is_softLep_jetelectron(const pat::PackedCandidate &jcand, edm::Handle<edm::View<pat::Electron> > electron_pat, const reco::Vertex& vtx){
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
bool BJetnessSelector::is_vetoPOGNoIPNoIso_jetelectron(const pat::PackedCandidate &jcand, edm::Handle<edm::View<pat::Electron> > electron_pat, const reco::Vertex& vtx){
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
bool BJetnessSelector::is_loosePOGNoIPNoIso_jetelectron(const pat::PackedCandidate &jcand, edm::Handle<edm::View<pat::Electron> > electron_pat, const reco::Vertex& vtx){
  bool isele = false;
  for(edm::View<pat::Electron>::const_iterator ele = electron_pat->begin(); ele != electron_pat->end(); ele++){
    const pat::Electron &lele = *ele; 
    if(deltaR(jcand.p4(),lele.p4())<0.1 && fabs(jcand.pt()-lele.pt())/lele.pt()<0.05) {
      double ooEmooP = 999;
      if(lele.ecalEnergy()==0)                   ooEmooP = 1e30;
      else if(!std::isfinite(lele.ecalEnergy())) ooEmooP = 1e30;
      else                                       ooEmooP = fabs(1.0/lele.ecalEnergy() - lele.eSuperClusterOverP()/lele.ecalEnergy() );
      if(lele.full5x5_sigmaIetaIeta()<0.0105 && fabs(lele.deltaEtaSuperClusterTrackAtVtx())<0.00976 && fabs(lele.deltaPhiSuperClusterTrackAtVtx())<0.0929
         && lele.hcalOverEcal()<0.0765 && ooEmooP<0.184
         //&& fabs((-1) * lele.gsfTrack()->dxy(vtx.position()))<0.0227 && fabs(lele.gsfTrack()->dz(vtx.position()))<0.379
         && lele.gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS)<=2 && lele.passConversionVeto()
        ) isele = true;
      if(isele) break;
    } 
  } 
  return isele; 
}
//Check that the track is a good track
bool BJetnessSelector::is_goodtrk(Track trk,const reco::Vertex& vtx){
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
TransientTrack BJetnessSelector::get_ttrk(Track trk, const TransientTrackBuilder& ttrkbuilder){
  TransientTrack ttrk;
  ttrk = ttrkbuilder.build(&trk);
  return ttrk;
}
vector<TransientTrack> BJetnessSelector::get_ttrks(vector<Track> trks, const TransientTrackBuilder& ttrkbuilder){
 vector<TransientTrack> ttrks;
 for(uint tr=0; tr<trks.size(); tr++){
  TransientTrack ttrk = ttrkbuilder.build(&trks[tr]);
  ttrks.push_back(ttrk);
 }
 return ttrks;
}
//Get transient pv from tracks
TransientVertex BJetnessSelector::get_tv(vector<Track> trks, const TransientTrackBuilder& ttrkbuilder){
  TransientVertex tv;
  vector<TransientTrack> ttrks = get_ttrks(trks,ttrkbuilder);
  if(ttrks.size()>=2)  tv = vtxFitter.vertex(ttrks);
  return tv;
}
//Get jet tracks info
void BJetnessSelector::get_jettrks(const pat::Jet& jet, const reco::Vertex& vtx, const TransientTrackBuilder& ttrkbuilder,
                                   vector<Track>& jetchtrks, vector<Track>& jetchtrkspv, vector<Track>& jetchtrksnpv, vector<tuple<double, double, double> >& jetsdir,
                                   edm::Handle<edm::View<pat::Electron> > electron_pat, edm::Handle<edm::View<pat::Muon> > muon_h,
                                   double& bjetness_num_pdgid_eles, double& bjetness_num_pdgid_mus, double& bjetness_num_soft_eles, double& bjetness_num_vetonoipnoiso_eles, double& bjetness_num_loosenoipnoiso_eles, double& bjetness_num_loose_mus
                                  ){
  //Access jet daughters
  vector<CandidatePtr> jdaus(jet.daughterPtrVector());
  sort(jdaus.begin(), jdaus.end(), [](const reco::CandidatePtr &p1, const reco::CandidatePtr &p2) {return p1->pt() > p2->pt();});
  for(uint jd=0; jd<jdaus.size(); jd++){ 
    const pat::PackedCandidate &jcand = dynamic_cast<const pat::PackedCandidate &>(*jdaus[jd]);
    if(deltaR(jcand.p4(),jet.p4())>0.4) continue;
    Track trk = Track(jcand.pseudoTrack());
    bool isgoodtrk = is_goodtrk(trk,vtx);
    //Minimal conditions for a track 
    if(isgoodtrk && jcand.charge()!=0 && jcand.fromPV()>1){
      if(fabs(jcand.pdgId())==13){
        bjetness_num_pdgid_mus++;
        if(is_loosePOG_jetmuon(jcand,muon_h)) bjetness_num_loose_mus++;
      } 
      if(fabs(jcand.pdgId())==11){
        bjetness_num_pdgid_eles++;
        if(is_softLep_jetelectron(jcand,electron_pat,vtx)) bjetness_num_soft_eles++;      
        if(is_vetoPOGNoIPNoIso_jetelectron(jcand,electron_pat,vtx)) bjetness_num_vetonoipnoiso_eles++;      
        if(is_loosePOGNoIPNoIso_jetelectron(jcand,electron_pat,vtx)) bjetness_num_loosenoipnoiso_eles++;      
      }
      jetchtrks.push_back(trk);
      //Other conditions on jet daughters
      //Using fromPV method
      if(jcand.fromPV()==pat::PackedCandidate::PVUsedInFit){
        jetchtrkspv.push_back(trk);
      }else{
        jetchtrksnpv.push_back(trk);
      }
     //Fill in the jet direction information to keep synchronisation
     jetsdir.push_back(make_tuple(jet.px(),jet.py(),jet.pz()));
    }//Ch trks 
  }//Loop on jet daus
}
//Get chi2 information from trk vertex
void BJetnessSelector::get_chi2(vector<Track> trks, const TransientTrackBuilder& ttrkbuilder, double& chi2red){
  if(trks.size()>=2){
    TransientVertex trks_tv = get_tv(trks, ttrkbuilder);
    if(trks_tv.isValid()) chi2red = trks_tv.totalChiSquared()/trks_tv.degreesOfFreedom(); 
  }
}
//Get Num of two-trk vertices
void BJetnessSelector::get_2trksinfo(vector<Track> trks, const TransientTrackBuilder& ttrkbuilder, double& num2v, double& numno2v, double& dca3d2t, double& dca3dno2t, double& dca2d2t, double& dca2dno2t){
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
pair<double,double> BJetnessSelector::dca2trks(Track tkA, Track tkB, const TransientTrackBuilder& ttrkbuilder){
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
void BJetnessSelector::get_avip3d(vector<Track> trks, const TransientTrackBuilder& ttrkbuilder, reco::Vertex vtx, vector<tuple<double, double, double> >& jetsdir, double& jetchtrks_avip3d_val, double& jetchtrks_avip3d_sig, double& jetchtrks_avsip3d_val, double& jetchtrks_avsip3d_sig){
  double valtemp = 0;
  vector<TransientTrack> ttrks = get_ttrks(trks,ttrkbuilder);
  for(uint t=0; t<ttrks.size(); t++){
    valtemp = IPTools::absoluteImpactParameter3D(ttrks[t],vtx).second.value();
    if(valtemp==valtemp) jetchtrks_avip3d_val  += valtemp;
    valtemp = IPTools::absoluteImpactParameter3D(ttrks[t],vtx).second.significance();
    if(valtemp==valtemp) jetchtrks_avip3d_sig  += valtemp;
    GlobalVector jetsdirgv(get<0>(jetsdir[t]),get<1>(jetsdir[t]),get<2>(jetsdir[t]));
    valtemp = IPTools::signedImpactParameter3D(ttrks[t],jetsdirgv,vtx).second.value();
    if(valtemp==valtemp) jetchtrks_avsip3d_val += valtemp;
    valtemp = IPTools::signedImpactParameter3D(ttrks[t],jetsdirgv,vtx).second.significance();
    if(valtemp==valtemp) jetchtrks_avsip3d_sig += valtemp;
  }
}
void BJetnessSelector::get_avip2d(vector<Track> trks, const TransientTrackBuilder& ttrkbuilder, reco::Vertex vtx, vector<tuple<double, double, double> >& jetsdir, double& jetchtrks_avip2d_val, double& jetchtrks_avip2d_sig, double& jetchtrks_avsip2d_val, double& jetchtrks_avsip2d_sig){
  double valtemp = 0;
  vector<TransientTrack> ttrks = get_ttrks(trks,ttrkbuilder);
  for(uint t=0; t<ttrks.size(); t++){
    valtemp = IPTools::absoluteTransverseImpactParameter(ttrks[t],vtx).second.value();
    if(valtemp==valtemp) jetchtrks_avip2d_val  += valtemp;
    valtemp = IPTools::absoluteTransverseImpactParameter(ttrks[t],vtx).second.significance();
    if(valtemp==valtemp) jetchtrks_avip2d_sig  += valtemp;
    GlobalVector jetsdirgv(get<0>(jetsdir[t]),get<1>(jetsdir[t]),get<2>(jetsdir[t]));
    valtemp = IPTools::signedTransverseImpactParameter(ttrks[t],jetsdirgv,vtx).second.value();
    if(valtemp==valtemp) jetchtrks_avsip2d_val += valtemp;
    valtemp = IPTools::signedTransverseImpactParameter(ttrks[t],jetsdirgv,vtx).second.significance();
    if(valtemp==valtemp) jetchtrks_avsip2d_sig += valtemp;
  }
}
void BJetnessSelector::get_avip1d(vector<Track> trks, const TransientTrackBuilder& ttrkbuilder, reco::Vertex vtx, vector<tuple<double, double, double> >& jetsdir, double& jetchtrks_avip1d_val, double& jetchtrks_avip1d_sig, double& jetchtrks_avsip1d_val, double& jetchtrks_avsip1d_sig){
  double valtemp = 0;
  vector<TransientTrack> ttrks = get_ttrks(trks,ttrkbuilder);
  SignedTransverseImpactParameter stip;
  for(uint t=0; t<ttrks.size(); t++){
    GlobalVector jetsdirgv(get<0>(jetsdir[t]),get<1>(jetsdir[t]),get<2>(jetsdir[t]));
    valtemp = fabs(stip.zImpactParameter(ttrks[t],jetsdirgv,vtx).second.value());
    if(valtemp==valtemp) jetchtrks_avip1d_val  += valtemp;
    valtemp = fabs(stip.zImpactParameter(ttrks[t],jetsdirgv,vtx).second.significance());
    if(valtemp==valtemp) jetchtrks_avip1d_sig  += valtemp;
    valtemp = stip.zImpactParameter(ttrks[t],jetsdirgv,vtx).second.value();
    if(valtemp==valtemp) jetchtrks_avsip1d_val += valtemp;
    valtemp = stip.zImpactParameter(ttrks[t],jetsdirgv,vtx).second.significance();
    if(valtemp==valtemp) jetchtrks_avsip1d_sig += valtemp; 
  }
}
//To do
//- Think if categories must be defined considering jets ordered by decreasing csv values
//- Include dR between loose muons and loose electrons
//- Include JER for jet selection
//- Check if electron isolation is included in the ID and how it is defined (as it is now or sa it changed for the muon) 
//- Having the vtx in this function is kept for historical reason 
//  and it will be kept until we do not have a final electron loose selection
//- Currently take events where jetcsv!=jetcsv, but need to understand what to do with them
//- Make sure the jet selection is exactly the same of the TTHbb analysis
//  Same PV, same jet cleaning from mu,ele
//- At the moment valtemp are not included in the average quantities, but the denominator is still the same
//  it should not change too much if the den is high, but check a better implementation
//- Pass maxjetnum from python according to the analysis categorisation
//- For the dca you are considering only the significance.
//  You may want to include the value as well, but probably after data/MC validation
//- Need to decide what to do precisely with 
//  else                       BJetness_npvTrkOVpvTrk.push_back(-997);
//  If npv!=0 this is likely a bjet evt. Check also the pT analogous variables
//- Fare l'average IP solo per le tracce secondarie o primarie?
//  Anche se l'effetto della bjetness si dovrebbe vedere con tutte le tracce
