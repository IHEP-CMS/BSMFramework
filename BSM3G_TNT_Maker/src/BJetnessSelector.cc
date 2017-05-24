#include "BSMFramework/BSM3G_TNT_Maker/interface/BJetnessSelector.h"
#include "BSMFramework/BSM3G_TNT_Maker/interface/ElectronPatSelector.h"

KalmanVertexFitter vtxFitter(true);
using namespace std;
using namespace pat;
using namespace edm;
using namespace reco;
BJetnessSelector::BJetnessSelector(std::string name, TTree* tree, bool debug, const pset& iConfig, edm::ConsumesCollector && ic):
  baseTree(name,tree,debug),
  electronVetoIdMapToken_(ic.consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electronVetoIdMap"))),
  electronLooseIdMapToken_(ic.consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electronLooseIdMap"))),
  electronMediumIdMapToken_(ic.consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electronMediumIdMap"))),
  electronTightIdMapToken_(ic.consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electronTightIdMap"))),
  eleMVATrigIdMapToken_(ic.consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleMVATrigIdMap"))),
  eleMVAnonTrigIdMap_(ic.consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleMVAnonTrigIdMap"))),
  eleMVATrigwp90IdMap_(ic.consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleMVATrigwp90IdMap"))),
  eleMVAnonTrigwp90IdMap_(ic.consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleMVAnonTrigwp90IdMap"))),
  eleHEEPIdMapToken_(ic.consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleHEEPIdMap"))),
  elemvaValuesMapToken_nonTrig_(ic.consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("elemvaValuesMap_nonTrig"))),
  elemvaCategoriesMapToken_nonTrig_(ic.consumes<edm::ValueMap<int> >(iConfig.getParameter<edm::InputTag>("elemvaCategoriesMap_nonTrig"))),
  elemvaValuesMapToken_Trig_(ic.consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("elemvaValuesMap_Trig"))),
  elemvaCategoriesMapToken_Trig_(ic.consumes<edm::ValueMap<int> >(iConfig.getParameter<edm::InputTag>("elemvaCategoriesMap_Trig")))
{
  vtx_h_              = ic.consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"));
  electron_pat_       = ic.consumes<edm::View<pat::Electron> >(iConfig.getParameter<edm::InputTag>("patElectrons"));
  muon_h_             = ic.consumes<edm::View<pat::Muon> >(iConfig.getParameter<edm::InputTag>("muons"));
  jets_               = ic.consumes<pat::JetCollection >(iConfig.getParameter<edm::InputTag>("jets"));
  rhopogHandle_       = ic.consumes<double>(edm::InputTag("fixedGridRhoFastjetAll"));
  rhoJERHandle_       = ic.consumes<double>(edm::InputTag("fixedGridRhoAll"));
  _vtx_ndof_min       = iConfig.getParameter<int>("vtx_ndof_min");
  _vtx_rho_max        = iConfig.getParameter<int>("vtx_rho_max");
  _vtx_position_z_max = iConfig.getParameter<double>("vtx_position_z_max");
  _is_data = iConfig.getParameter<bool>("is_data");
  if(!_is_data) prunedGenToken_ = ic.consumes<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("pruned"));
  jecPayloadNamesAK4PFchsMC1_   = iConfig.getParameter<edm::FileInPath>("jecPayloadNamesAK4PFchsMC1");
  jecPayloadNamesAK4PFchsMC2_   = iConfig.getParameter<edm::FileInPath>("jecPayloadNamesAK4PFchsMC2");
  jecPayloadNamesAK4PFchsMC3_   = iConfig.getParameter<edm::FileInPath>("jecPayloadNamesAK4PFchsMC3");
  jecPayloadNamesAK4PFchsMCUnc_ = iConfig.getParameter<edm::FileInPath>("jecPayloadNamesAK4PFchsMCUnc");
  jecPayloadNamesAK4PFchsDATA1_   = iConfig.getParameter<edm::FileInPath>("jecPayloadNamesAK4PFchsDATA1");
  jecPayloadNamesAK4PFchsDATA2_   = iConfig.getParameter<edm::FileInPath>("jecPayloadNamesAK4PFchsDATA2");
  jecPayloadNamesAK4PFchsDATA3_   = iConfig.getParameter<edm::FileInPath>("jecPayloadNamesAK4PFchsDATA3");
  jecPayloadNamesAK4PFchsDATA4_   = iConfig.getParameter<edm::FileInPath>("jecPayloadNamesAK4PFchsDATA4");
  jecPayloadNamesAK4PFchsDATAUnc_ = iConfig.getParameter<edm::FileInPath>("jecPayloadNamesAK4PFchsDATAUnc");
  jerAK4PFchs_     = iConfig.getParameter<edm::FileInPath>("jerAK4PFchs").fullPath();
  jerAK4PFchsSF_   = iConfig.getParameter<edm::FileInPath>("jerAK4PFchsSF").fullPath();
  JECInitialization();
  SetBranches();
  // Helper class of helper class. Registration to framework required.
  // consumesCollector may be static - one per helper class.
  // Registration with framework needs to be done using consumesCollector object
  // that is already instantiated for class.
  // 'ic' here is 'edm::ConsumesCollector && ic' - which is an lvalue.
  // It has a name and you can take its address. Its just the type of that lvalue is rvalue reference to key.
  // You need to pass in an rvalue - use std::move.
  evSel = new TTHbb_eventSelector("miniAOD", tree, debug, iConfig, std::move(ic));
}
BJetnessSelector::~BJetnessSelector(){
  delete tree_;
}
void BJetnessSelector::Fill(const edm::Event& iEvent, const edm::EventSetup& iSetup, int& bjetnesssel_filter){
  Clear();
  if(debug_) std::cout << ("BJetnessSelector:Fill") << std::endl;
  int tmp_evtNum     = iEvent.id().event();
  if(debug_) std::cout << ("BJetnessSelector:Fill") << ": Event Number = " << tmp_evtNum << std::endl;
  /////
  //   Recall collections
  /////
  edm::Handle<reco::VertexCollection> vtx_h;
  iEvent.getByToken(vtx_h_, vtx_h);
  edm::Handle<edm::View<pat::Electron> > electron_pat;
  iEvent.getByToken(electron_pat_, electron_pat);
  edm::Handle<edm::View<pat::Muon> > muon_h;
  iEvent.getByToken(muon_h_, muon_h);
  edm::Handle<pat::JetCollection> jets;
  iEvent.getByToken(jets_, jets);
  edm::Handle<double> rhopogHandle;
  iEvent.getByToken(rhopogHandle_,rhopogHandle);
  double rhopog = *rhopogHandle;
  edm::Handle<double> rhoJERHandle;
  iEvent.getByToken(rhoJERHandle_,rhoJERHandle);
  double rhoJER = *rhoJERHandle;
  edm::Handle<edm::ValueMap<bool>  > veto_id_decisions;
  edm::Handle<edm::ValueMap<bool>  > loose_id_decisions;
  edm::Handle<edm::ValueMap<bool>  > medium_id_decisions;
  edm::Handle<edm::ValueMap<bool>  > tight_id_decisions;
  edm::Handle<edm::ValueMap<bool>  > heep_id_decisions;
  edm::Handle<edm::ValueMap<bool>  > mvatrig_id_decisions;
  edm::Handle<edm::ValueMap<bool>  > mvanontrig_id_decisions;
  edm::Handle<edm::ValueMap<bool>  > mvatrigwp90_id_decisions;
  edm::Handle<edm::ValueMap<bool>  > mvanontrigwp90_id_decisions;
  edm::Handle<edm::ValueMap<float> > elemvaValues_nonTrig;
  edm::Handle<edm::ValueMap<int> >   elemvaCategories_nonTrig;
  edm::Handle<edm::ValueMap<float> > elemvaValues_Trig;
  edm::Handle<edm::ValueMap<int> >   elemvaCategories_Trig;
  iEvent.getByToken(electronVetoIdMapToken_,   veto_id_decisions);
  iEvent.getByToken(electronLooseIdMapToken_,  loose_id_decisions);
  iEvent.getByToken(electronMediumIdMapToken_, medium_id_decisions);
  iEvent.getByToken(electronTightIdMapToken_,  tight_id_decisions);
  iEvent.getByToken(eleMVATrigIdMapToken_,     mvatrig_id_decisions);
  iEvent.getByToken(eleMVAnonTrigIdMap_,       mvanontrig_id_decisions);
  iEvent.getByToken(eleMVATrigwp90IdMap_,      mvatrigwp90_id_decisions);
  iEvent.getByToken(eleMVAnonTrigwp90IdMap_,   mvanontrigwp90_id_decisions);
  iEvent.getByToken(eleHEEPIdMapToken_,        heep_id_decisions);
  iEvent.getByToken(elemvaValuesMapToken_nonTrig_,     elemvaValues_nonTrig);
  iEvent.getByToken(elemvaCategoriesMapToken_nonTrig_, elemvaCategories_nonTrig);
  iEvent.getByToken(elemvaValuesMapToken_Trig_,        elemvaValues_Trig);
  iEvent.getByToken(elemvaCategoriesMapToken_Trig_,    elemvaCategories_Trig);
  edm::ESHandle<TransientTrackBuilder> ttrkbuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",ttrkbuilder);

  bool isSL_ = evSel->singlelepton_channel(iEvent,iSetup);
  bool isDL_ = evSel->dilepton_channel(iEvent,iSetup);
  if(debug_) std::cout << ("BJetnessSelector::Fill") << ": isSL_ " << isSL_ << std::endl;
  if(debug_) std::cout << ("BJetnessSelector::Fill") << ": isDL_ " << isDL_ << std::endl;

  //Running only Single lepton selection for BJetness variables atm.
  if (!isSL_ && !isDL_) return;
  BJetness_isSingleLepton = isSL_;
  BJetness_isDoubleLepton = isDL_;

  if(debug_) std::cout << ("BJetnessSelector::Fill") << " Passed SL Selection " << std::endl;

  const reco::Vertex &PV = vtx_h->front();

  if(debug_) std::cout << ("BJetnessSelector::Fill") << " Passed SL Selection " << std::endl;

  /////
  //   Get the good jets of the event
  /////
  //Iterate to access jet by decreasing b-tagging value
  //int jet_pos = 0; //This counter helps to order jets
  int jet_num = 0; //This counter accounts for the number of good jets in the events
                   //The definition of good jet in the event must be the same of the TTHbb analysis
                   //so that jet_num corresponds to the number of jets that define the categories in the TTHbb search
  int jetb_num = 0;
  int DL_jet_num = 0;
  int DL_jetb_num = 0;
  vector<pair<double,int> > leading_jet_csv_pos;
  vector<pair<double,int> > subleading_jet_csv_pos;

  std::vector<const reco::Candidate*> looseleps;
  std::vector<const reco::Candidate*> loose_electrons;
  std::vector<const reco::Candidate*> loose_muons;
  std::vector<const reco::Candidate*> tightleps;
  std::vector<const reco::Candidate*> tight_electrons;
  std::vector<const reco::Candidate*> tight_muons;
  std::vector<const pat::Jet*> subleading_jets;
  std::vector<const pat::Jet*> leading_jets;
  std::vector<const pat::Jet*> subleading_btags;
  std::vector<const pat::Jet*> leading_btags;

  if(debug_) std::cout << ("BJetnessSelector::Fill") << " FillLeptonVectors " << std::endl;
  evSel->FillLeptonVectors(iEvent, looseleps, loose_electrons, loose_muons, tightleps, tight_electrons, tight_muons);
  if(debug_) std::cout << ("BJetnessSelector::Fill") << " FillJetVectors " << std::endl;
  evSel->FillJetVectors(iEvent, looseleps, leading_jets, subleading_jets, leading_btags, subleading_btags, leading_jet_csv_pos, subleading_jet_csv_pos);

  jet_num = leading_jets.size();
  jetb_num = leading_btags.size();

  DL_jet_num = leading_jets.size();
  DL_jetb_num = subleading_btags.size();

  if(debug_) std::cout << ("BJetnessSelector:Fill") << "DL_jet_num: " << DL_jet_num << " DL_jetb_num: " << DL_jetb_num << std::endl;

  //=======================
  //   Event Selection
  //=======================
  if(debug_) std::cout << ("BJetnessSelector:Fill") << ": TTH selection" << std::endl;
  if(debug_) std::cout << ("BJetnessSelector:Fill") << ": tightleps.size(): " << tightleps.size() << std::endl;
  if(debug_) std::cout << ("BJetnessSelector:Fill") << ": looseleps.size() : " << looseleps.size() << std::endl;
  if(debug_) std::cout << ("BJetnessSelector:Fill") << ": jet_num: " << jet_num << std::endl;
  if(debug_) std::cout << ("BJetnessSelector:Fill") << ": jetb_num : " << jetb_num << std::endl;
  if(debug_) std::cout << ("BJetnessSelector:Fill") << ": leading_jet_csv_pos.size(): " << leading_jet_csv_pos.size() << std::endl;

  //if(tightleps.size()==1 && looseleps.size()==1 && jet_num>=4 && jetb_num>=2){
    //BJetness_isSingleLepton = 1;
  //  bjetnesssel_filter      = 1;
  //}
  //if(bjetnesssel_filter==0) return; //This line must be comments if bjetnessselfilter is false in the config file

  //if(tightleps.size()==1 && looseleps.size()==2 && DL_jet_num>=2 && DL_jetb_num>=1) BJetness_isDoubleLepton = 1;

  //===========================
  //     BJetness Variables
  //===========================
  // NOTE: Currently only calculating BJetness variables for SL selection.
  //       Hence, only using leading_jet_csv_pos containing the csv and index of the leading jet in the jet* container.

  if(debug_) std::cout << ("BJetnessSelector::Fill") << " Calculate BJetness variables . . . " << std::endl;
  sort(leading_jet_csv_pos.rbegin(), leading_jet_csv_pos.rend());
  if(jet_num!=0){
    //cout<<"Num of jet is"<<setw(20)<<jet_num<<" "<<jetb_num<<endl;
    BJetness_numjet.push_back(jet_num);
    for(int jn=0; jn<jet_num; jn++){
      if(debug_) std::cout << ("BJetnessSelector:Fill") << ": jet index = " << jn << std::endl;

      const pat::Jet & j = (*jets)[leading_jet_csv_pos[jn].second];
      //Kinematics and csv
      //Jet Energy Corrections and Uncertainties
      int vtxsize = vtx_h->size();
      double corrAK4PFchs     = 1;
      if(debug_) std::cout << ("BJetnessSelector:Fill") << ": JECs" << std::endl;
      reco::Candidate::LorentzVector uncorrJetAK4PFchs = j.correctedP4(0);
      if(!_is_data){
        if(debug_) std::cout << ("BJetnessSelector:Fill") << ": MC JECs" << std::endl;
        jecAK4PFchsMC_->setJetEta( uncorrJetAK4PFchs.eta()    );
        jecAK4PFchsMC_->setJetPt ( uncorrJetAK4PFchs.pt()     );
        jecAK4PFchsMC_->setJetE  ( uncorrJetAK4PFchs.energy() );
        jecAK4PFchsMC_->setRho	( rhopog  );
        jecAK4PFchsMC_->setNPV	( vtxsize  );
        jecAK4PFchsMC_->setJetA  ( j.jetArea()	     );
        corrAK4PFchs = jecAK4PFchsMC_->getCorrection();
      } else {
        if(debug_) std::cout << ("BJetnessSelector:Fill") << ": DATA JECs" << std::endl;
        jecAK4PFchsDATA_->setJetEta( uncorrJetAK4PFchs.eta()    );
        jecAK4PFchsDATA_->setJetPt ( uncorrJetAK4PFchs.pt()     );
        jecAK4PFchsDATA_->setJetE  ( uncorrJetAK4PFchs.energy() );
        jecAK4PFchsDATA_->setRho	( rhopog  );
        jecAK4PFchsDATA_->setNPV	( vtxsize  );
        jecAK4PFchsDATA_->setJetA  ( j.jetArea()	     );
        corrAK4PFchs = jecAK4PFchsDATA_->getCorrection();
      }
      if(debug_) std::cout << ("BJetnessSelector:Fill") << ": JERs"<< std::endl;
      float JERScaleFactor     = 1;
      float JERScaleFactorUP   = 1;
      float JERScaleFactorDOWN = 1;
      if(!_is_data) GetJER(j, corrAK4PFchs, rhoJER, true, JERScaleFactor, JERScaleFactorUP, JERScaleFactorDOWN);
      double jetpt = (j.correctedJet("Uncorrected").pt()*corrAK4PFchs*JERScaleFactor);
      //The code related to "Jet Energy Corrections and Uncertainties" should go in a function
      if(debug_) std::cout << ("BJetnessSelector::Fill") << " Set variables from jets " << std::endl;
      BJetness_jetpt.push_back(jetpt);
      BJetness_jeteta.push_back(j.eta());
      BJetness_jetphi.push_back(j.phi());
      BJetness_jetenergy.push_back(j.correctedJet("Uncorrected").energy()*corrAK4PFchs*JERScaleFactor); //Note that this is not correct by JES/JER (expect if JES in globalTag are ok)
      double jetcsv = j.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
      if(jetcsv!=jetcsv) jetcsv = -996;
      BJetness_jetcsv.push_back(jetcsv);
      double jetbjetprob = j.bDiscriminator("pfJetProbabilityBJetTags");
      if(jetbjetprob!=jetbjetprob) jetbjetprob = -996;
      BJetness_pfJetProbabilityBJetTags.push_back(jetbjetprob);
      double jetbjetcombmva = j.bDiscriminator("pfCombinedMVAV2BJetTags");
      if(jetbjetcombmva!=jetbjetcombmva) jetbjetcombmva = -996;
      BJetness_pfCombinedMVAV2BJetTags.push_back(jetbjetcombmva);
      double jetpfCombinedCvsLJetTags = j.bDiscriminator("pfCombinedCvsLJetTags");
      if(jetpfCombinedCvsLJetTags!=jetpfCombinedCvsLJetTags) jetpfCombinedCvsLJetTags = -996;
      BJetness_pfCombinedCvsLJetTags.push_back(jetpfCombinedCvsLJetTags);
      double jetpfCombinedCvsBJetTags = j.bDiscriminator("pfCombinedCvsBJetTags");
      if(jetpfCombinedCvsBJetTags!=jetpfCombinedCvsBJetTags) jetpfCombinedCvsBJetTags = -996;
      BJetness_pfCombinedCvsBJetTags.push_back(jetpfCombinedCvsBJetTags);
      //cout<<"Jet pf"<<setw(20)<<"hf"<<setw(20)<<"pt"<<setw(20)<<"eta"<<setw(20)<<"phi"<<setw(20)<<"energy"<<setw(20)<<"csv"<<endl;
      //cout<<j.partonFlavour()<<setw(20)<<j.hadronFlavour()<<setw(20)<<j.pt()<<setw(20)<<j.eta()<<setw(20)<<j.phi()<<setw(20)<<j.energy()<<setw(20)<<j.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags")<<endl;
      //cout<<"CSV "<<setw(20)<<j<<setw(20)<<j.bDiscriminator("pfJetProbabilityBJetTags")<<setw(20)<<j.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags")<<setw(20)<<j.bDiscriminator("pfCombinedMVAV2BJetTags")<<endl;
      BJetness_partonFlavour.push_back(j.partonFlavour());
      BJetness_hadronFlavour.push_back(j.hadronFlavour());
    }
    /////
    //   Gen info
    /////
    if(debug_) std::cout << ("BJetnessSelector::Fill") << " Gen info " << std::endl;

    if(!_is_data){
      BJetness_ngenbh = 0;
      BJetness_ngenbt = 0;
      BJetness_ngenb  = 0;
      BJetness_ngenc  = 0;
      Handle<edm::View<reco::GenParticle> > pruned;
      iEvent.getByToken(prunedGenToken_, pruned);
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
    double bjetness_pt  = 0;
    double bjetness_eta = 0;
    double bjetness_phi = 0;
    double bjetness_en  = 0;
    vector<Track> jetschtrks; jetschtrks.clear();
    vector<Track> jetschtrkspv; jetschtrkspv.clear();
    vector<Track> jetschtrksnpv; jetschtrksnpv.clear();
    vector<tuple<double, double, double> > jetsdir; jetsdir.clear();
    double bjetness_num_pdgid_eles          = 0;
    double bjetness_num_soft_eles           = 0;
    double bjetness_num_vetonoipnoiso_eles  = 0;
    double bjetness_num_loosenoipnoiso_eles = 0;
    double bjetness_num_veto_eles           = 0;
    double bjetness_num_loose_eles          = 0;
    double bjetness_num_medium_eles         = 0;
    double bjetness_num_tight_eles          = 0;
    double bjetness_num_mvatrig_eles        = 0;
    double bjetness_num_mvanontrig_eles     = 0;
    double bjetness_num_mvatrigwp90_eles    = 0;
    double bjetness_num_mvanontrigwp90_eles = 0;
    double bjetness_num_heep_eles           = 0;
    double bjetness_num_pdgid_mus           = 0;
    double bjetness_num_loose_mus           = 0;
    double bjetness_num_soft_mus            = 0;
    double bjetness_num_medium_mus          = 0;
    double bjetness_num_tight_mus           = 0;
    double bjetness_num_highpt_mus          = 0;
    double bjetness_num_POGisGood_mus       = 0;
    int maxjetnum = 6; //Pass this value from python at some point
                       //The values will have to be chosen depending on the definition of the categories
    if(jet_num<maxjetnum)  maxjetnum = jet_num;
    for(int jn=1; jn<maxjetnum; jn++){
      const pat::Jet & j = (*jets)[leading_jet_csv_pos[jn].second];
      bjetness_pt  += j.pt();
      bjetness_eta += j.eta();
      bjetness_phi += j.phi();
      bjetness_en  += j.energy();
      // Uses track info to build bjetness variables (passed into method by ref.)
      get_jettrks(j, PV, *ttrkbuilder,
                  jetschtrks, jetschtrkspv, jetschtrksnpv, jetsdir,
                  iEvent, electron_pat, muon_h,
                  bjetness_num_pdgid_eles, bjetness_num_soft_eles, bjetness_num_vetonoipnoiso_eles, bjetness_num_loosenoipnoiso_eles, bjetness_num_veto_eles, bjetness_num_loose_eles, bjetness_num_medium_eles, bjetness_num_tight_eles, bjetness_num_mvatrig_eles, bjetness_num_mvanontrig_eles, bjetness_num_mvatrigwp90_eles, bjetness_num_mvanontrigwp90_eles, bjetness_num_heep_eles, bjetness_num_pdgid_mus, bjetness_num_loose_mus, bjetness_num_soft_mus, bjetness_num_medium_mus, bjetness_num_tight_mus, bjetness_num_highpt_mus, bjetness_num_POGisGood_mus
                 );
    }
    //PFCand info
    /////
    //   Get info to evaluate the event BJetness
    /////
    BJetness_pt.push_back(bjetness_pt);
    BJetness_eta.push_back(bjetness_eta);
    BJetness_phi.push_back(bjetness_phi);
    BJetness_en.push_back(bjetness_en);
    BJetness_ptOVen.push_back(bjetness_pt/bjetness_en);
    BJetness_num_pdgid_eles.push_back(bjetness_num_pdgid_eles);
    BJetness_num_soft_eles.push_back(bjetness_num_soft_eles);
    BJetness_num_vetonoipnoiso_eles.push_back(bjetness_num_vetonoipnoiso_eles);
    BJetness_num_loosenoipnoiso_eles.push_back(bjetness_num_loosenoipnoiso_eles);
    BJetness_num_veto_eles.push_back(bjetness_num_veto_eles);
    BJetness_num_loose_eles.push_back(bjetness_num_loose_eles);
    BJetness_num_medium_eles.push_back(bjetness_num_medium_eles);
    BJetness_num_tight_eles.push_back(bjetness_num_tight_eles);
    BJetness_num_mvatrig_eles.push_back(bjetness_num_mvatrig_eles);
    BJetness_num_mvanontrig_eles.push_back(bjetness_num_mvanontrig_eles);
    BJetness_num_mvatrigwp90_eles.push_back(bjetness_num_mvatrigwp90_eles);
    BJetness_num_mvanontrigwp90_eles.push_back(bjetness_num_mvanontrigwp90_eles);
    BJetness_num_heep_eles.push_back(bjetness_num_heep_eles);
    BJetness_num_pdgid_mus.push_back(bjetness_num_pdgid_mus);
    BJetness_num_loose_mus.push_back(bjetness_num_loose_mus);
    BJetness_num_soft_mus.push_back(bjetness_num_soft_mus);
    BJetness_num_medium_mus.push_back(bjetness_num_medium_mus);
    BJetness_num_tight_mus.push_back(bjetness_num_tight_mus);
    BJetness_num_highpt_mus.push_back(bjetness_num_highpt_mus);
    BJetness_num_POGisGood_mus.push_back(bjetness_num_POGisGood_mus);
    BJetness_numjettrks.push_back(jetschtrks.size());
    BJetness_numjettrkspv.push_back(jetschtrkspv.size());
    BJetness_numjettrksnopv.push_back(jetschtrksnpv.size());
    //cout<<setw(20)<<"Num_of_trks"<<setw(20)<<jetschtrks.size()<<setw(20)<<jetschtrkspv.size()<<setw(20)<<jetschtrksnpv.size()<<endl;
    if(jetschtrks.size()!=0){
      //Num of trks and their pT/En
      BJetness_npvTrkOVcollTrk.push_back(double(jetschtrksnpv.size())/double(jetschtrks.size()));
      BJetness_pvTrkOVcollTrk.push_back(double(jetschtrkspv.size())/double(jetschtrks.size()));
      if(jetschtrkspv.size()!=0) BJetness_npvTrkOVpvTrk.push_back(double(jetschtrksnpv.size())/double(jetschtrkspv.size()));
      else                       BJetness_npvTrkOVpvTrk.push_back(-997);
      double ptjettrks    = 0;
      for(uint jt=0; jt<jetschtrks.size(); jt++) ptjettrks += jetschtrks[jt].pt();
      double ptjettrkspv    = 0;
      for(uint jt=0; jt<jetschtrkspv.size(); jt++) ptjettrkspv += jetschtrkspv[jt].pt();
      double ptjettrksnpv    = 0;
      for(uint jt=0; jt<jetschtrksnpv.size(); jt++) ptjettrksnpv += jetschtrksnpv[jt].pt();
      if(ptjettrks!=0) BJetness_npvPtOVcollPt.push_back(ptjettrksnpv/ptjettrks);
      else             BJetness_npvPtOVcollPt.push_back(-997);
      if(ptjettrks!=0) BJetness_pvPtOVcollPt.push_back(ptjettrkspv/ptjettrks);
      else             BJetness_pvPtOVcollPt.push_back(-997);
      if(ptjettrkspv!=0) BJetness_npvPtOVpvPt.push_back(ptjettrksnpv/ptjettrkspv);
      else               BJetness_npvPtOVpvPt.push_back(-997);
      //Trk prop rel to jet dir
      double jetchtrks_avprel   = 0;
      double jetchtrks_avppar   = 0;
      double jetchtrks_avetarel = 0;
      double jetchtrks_avetapar = 0;
      double jetchtrks_avdr     = 0;
      get_avreljet(jetschtrks,jetsdir,jetchtrks_avprel,jetchtrks_avppar,jetchtrks_avetarel,jetchtrks_avetapar,jetchtrks_avdr);
      BJetness_avprel.push_back(jetchtrks_avprel/double(jetschtrks.size()));
      BJetness_avppar.push_back(jetchtrks_avppar/double(jetschtrks.size()));
      BJetness_avetarel.push_back(jetchtrks_avetarel/double(jetschtrks.size()));
      BJetness_avetapar.push_back(jetchtrks_avetapar/double(jetschtrks.size()));
      BJetness_avdr.push_back(jetchtrks_avdr/double(jetschtrks.size()));
      BJetness_avpreljetpt.push_back(jetchtrks_avprel/bjetness_pt);
      BJetness_avpreljeten.push_back(jetchtrks_avprel/bjetness_en);
      BJetness_avpparjetpt.push_back(jetchtrks_avppar/bjetness_pt);
      BJetness_avpparjeten.push_back(jetchtrks_avppar/bjetness_en);
      //cout<<jetchtrks_avprel/double(jetschtrks.size())<<" "<<jetchtrks_avppar/double(jetschtrks.size())<<" "<<jetchtrks_avetarel/double(jetschtrks.size())<<" "<<jetchtrks_avetapar/double(jetschtrks.size())<<" "<<jetchtrks_avdr/double(jetschtrks.size())<<endl;
      //Two_trk_info (we may want to use one function per variable. Probably more clear, but slower!)
      double jetchtrks_num2v     = 0;
      double jetchtrks_numno2v   = 0;
      double jetchtrks_dca3d2t   = 0;
      double jetchtrks_dca3dno2t = 0;
      double jetchtrks_dca2d2t   = 0;
      double jetchtrks_dca2dno2t = 0;
      //Not working from 80X on?
      //get_2trksinfo(jetschtrks, *ttrkbuilder, jetchtrks_num2v, jetchtrks_numno2v, jetchtrks_dca3d2t, jetchtrks_dca3dno2t, jetchtrks_dca2d2t, jetchtrks_dca2dno2t);
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
      //chi2
      double jetchtrks_chi2 = 997;
      //Not working from 80X on?
      //get_chi2(jetschtrks, *ttrkbuilder, jetchtrks_chi2);
      BJetness_chi2.push_back(jetchtrks_chi2);
      //ImpactParameter
      double ip_valtemp = 0;
      double jetchtrks_avip3d_val  = 0;
      double jetchtrks_avip3d_sig  = 0;
      double jetchtrks_avsip3d_val = 0;
      double jetchtrks_avsip3d_sig = 0;
      double jetchtrks_numip3dpos  = 0;
      double jetchtrks_numip3dneg  = 0;
      get_avip3d(jetschtrks, *ttrkbuilder, PV, jetsdir, jetchtrks_avip3d_val,jetchtrks_avip3d_sig,jetchtrks_avsip3d_val,jetchtrks_avsip3d_sig,jetchtrks_numip3dpos,jetchtrks_numip3dneg);
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
      BJetness_numip3dpos.push_back(jetchtrks_numip3dpos);
      BJetness_numip3dneg.push_back(jetchtrks_numip3dneg);
      double jetchtrks_avip2d_val  = 0;
      double jetchtrks_avip2d_sig  = 0;
      double jetchtrks_avsip2d_val = 0;
      double jetchtrks_avsip2d_sig = 0;
      double jetchtrks_numip2dpos  = 0;
      double jetchtrks_numip2dneg  = 0;
      get_avip2d(jetschtrks, *ttrkbuilder, PV, jetsdir, jetchtrks_avip2d_val,jetchtrks_avip2d_sig,jetchtrks_avsip2d_val,jetchtrks_avsip2d_sig,jetchtrks_numip2dpos,jetchtrks_numip2dneg);
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
      BJetness_numip2dpos.push_back(jetchtrks_numip2dpos);
      BJetness_numip2dneg.push_back(jetchtrks_numip2dneg);
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
      //Trk prop rel to jet dir
      BJetness_avprel.push_back(-998);
      BJetness_avppar.push_back(-998);
      BJetness_avetarel.push_back(-998);
      BJetness_avetapar.push_back(-998);
      BJetness_avdr.push_back(-998);
      BJetness_avpreljetpt.push_back(-998);
      BJetness_avpreljeten.push_back(-998);
      BJetness_avpparjetpt.push_back(-998);
      BJetness_avpparjeten.push_back(-998);
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
      BJetness_numip3dpos.push_back(-998);
      BJetness_numip3dneg.push_back(-998);
      BJetness_avip2d_val.push_back(-998);
      BJetness_avip2d_sig.push_back(-998);
      BJetness_avsip2d_val.push_back(-998);
      BJetness_avsip2d_sig.push_back(-998);
      BJetness_numip2dpos.push_back(-998);
      BJetness_numip2dneg.push_back(-998);
      BJetness_avip1d_val.push_back(-998);
      BJetness_avip1d_sig.push_back(-998);
      BJetness_avsip1d_val.push_back(-998);
      BJetness_avsip1d_sig.push_back(-998);
    }
 }else{
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
   BJetness_pfCombinedMVAV2BJetTags.push_back(-999);
   BJetness_pfCombinedCvsLJetTags.push_back(-999);
   BJetness_pfCombinedCvsBJetTags.push_back(-999);
   BJetness_pt.push_back(-999);
   BJetness_eta.push_back(-999);
   BJetness_phi.push_back(-999);
   BJetness_en.push_back(-999);
   BJetness_ptOVen.push_back(-999);
   //PFcand info
   BJetness_jetschpvass.push_back(-999);
   BJetness_jetschfrompv.push_back(-999);
   BJetness_jetschip3dval.push_back(-999);
   BJetness_jetschip3dsig.push_back(-999);
   BJetness_jetschip2dval.push_back(-999);
   BJetness_jetschip2dsig.push_back(-999);
   BJetness_jetschisgoodtrk.push_back(-999);
   BJetness_jetschtrkpur.push_back(-999);
   BJetness_jetschpt.push_back(-999);
   BJetness_jetscheta.push_back(-999);
   BJetness_jetschen.push_back(-999);
   //Num_of_trks
   BJetness_num_pdgid_eles.push_back(-999);
   BJetness_num_soft_eles.push_back(-999);
   BJetness_num_vetonoipnoiso_eles.push_back(-999);
   BJetness_num_loosenoipnoiso_eles.push_back(-999);
   BJetness_num_veto_eles.push_back(-999);
   BJetness_num_loose_eles.push_back(-999);
   BJetness_num_medium_eles.push_back(-999);
   BJetness_num_tight_eles.push_back(-999);
   BJetness_num_mvatrig_eles.push_back(-999);
   BJetness_num_mvanontrig_eles.push_back(-999);
   BJetness_num_mvatrigwp90_eles.push_back(-999);
   BJetness_num_mvanontrigwp90_eles.push_back(-999);
   BJetness_num_heep_eles.push_back(-999);
   BJetness_num_pdgid_mus.push_back(-999);
   BJetness_num_loose_mus.push_back(-999);
   BJetness_num_soft_mus.push_back(-999);
   BJetness_num_medium_mus.push_back(-999);
   BJetness_num_tight_mus.push_back(-999);
   BJetness_num_highpt_mus.push_back(-999);
   BJetness_num_POGisGood_mus.push_back(-999);
   BJetness_numjettrks.push_back(-999);
   BJetness_numjettrkspv.push_back(-999);
   BJetness_numjettrksnopv.push_back(-999);
   BJetness_npvTrkOVcollTrk.push_back(-999);
   BJetness_pvTrkOVcollTrk.push_back(-999);
   BJetness_npvTrkOVpvTrk.push_back(-999);
   BJetness_npvPtOVcollPt.push_back(-999);
   BJetness_pvPtOVcollPt.push_back(-999);
   BJetness_npvPtOVpvPt.push_back(-999);
   //Trk prop rel to jet dir
   BJetness_avprel.push_back(-999);
   BJetness_avppar.push_back(-999);
   BJetness_avetarel.push_back(-999);
   BJetness_avetapar.push_back(-999);
   BJetness_avdr.push_back(-999);
   BJetness_avpreljetpt.push_back(-999);
   BJetness_avpreljeten.push_back(-999);
   BJetness_avpparjetpt.push_back(-999);
   BJetness_avpparjeten.push_back(-999);
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
   BJetness_numip3dpos.push_back(-999);
   BJetness_numip3dneg.push_back(-999);
   BJetness_avip2d_val.push_back(-999);
   BJetness_avip2d_sig.push_back(-999);
   BJetness_avsip2d_val.push_back(-999);
   BJetness_avsip2d_sig.push_back(-999);
   BJetness_numip2dpos.push_back(-999);
   BJetness_numip2dneg.push_back(-999);
   BJetness_avip1d_val.push_back(-999);
   BJetness_avip1d_sig.push_back(-999);
   BJetness_avsip1d_val.push_back(-999);
   BJetness_avsip1d_sig.push_back(-999);
 }
 if(debug_) std::cout << ("BJetnessSelector:Fill") << " End of Fill " << std::endl;
 return;

}
void BJetnessSelector::SetBranches(){
  if(debug_) std::cout<<"setting branches: calling AddBranch of baseTree"<<std::endl;
  //Evt selection
  AddBranch(&BJetness_isSingleLepton           ,"BJetness_isSingleLepton");
  AddBranch(&BJetness_isDoubleLepton           ,"BJetness_isDoubleLepton");
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
  AddBranch(&BJetness_pfCombinedMVAV2BJetTags  ,"BJetness_pfCombinedMVAV2BJetTags");
  AddBranch(&BJetness_pfCombinedCvsLJetTags    ,"BJetness_pfCombinedCvsLJetTags");
  AddBranch(&BJetness_pfCombinedCvsBJetTags    ,"BJetness_pfCombinedCvsBJetTags");
  AddBranch(&BJetness_pt                       ,"BJetness_pt");
  AddBranch(&BJetness_eta                      ,"BJetness_eta");
  AddBranch(&BJetness_phi                      ,"BJetness_phi");
  AddBranch(&BJetness_en                       ,"BJetness_en");
  AddBranch(&BJetness_ptOVen                   ,"BJetness_ptOVen");
  //PFcand info
  AddBranch(&BJetness_jetschpvass              ,"BJetness_jetschpvass");
  AddBranch(&BJetness_jetschfrompv             ,"BJetness_jetschfrompv");
  AddBranch(&BJetness_jetschip3dval            ,"BJetness_jetschip3dval");
  AddBranch(&BJetness_jetschip3dsig            ,"BJetness_jetschip3dsig");
  AddBranch(&BJetness_jetschip2dval            ,"BJetness_jetschip2dval");
  AddBranch(&BJetness_jetschip2dsig            ,"BJetness_jetschip2dsig");
  AddBranch(&BJetness_jetschisgoodtrk          ,"BJetness_jetschisgoodtrk");
  AddBranch(&BJetness_jetschtrkpur             ,"BJetness_jetschtrkpur");
  AddBranch(&BJetness_jetschpt                 ,"BJetness_jetschpt");
  AddBranch(&BJetness_jetscheta                ,"BJetness_jetscheta");
  AddBranch(&BJetness_jetschen                 ,"BJetness_jetschen");
  //Num_of_trks
  AddBranch(&BJetness_num_pdgid_eles           ,"BJetness_num_pdgid_eles");
  AddBranch(&BJetness_num_soft_eles            ,"BJetness_num_soft_eles");
  AddBranch(&BJetness_num_vetonoipnoiso_eles   ,"BJetness_num_vetonoipnoiso_eles");
  AddBranch(&BJetness_num_loosenoipnoiso_eles  ,"BJetness_num_loosenoipnoiso_eles");
  AddBranch(&BJetness_num_veto_eles            ,"BJetness_num_veto_eles");
  AddBranch(&BJetness_num_loose_eles           ,"BJetness_num_loose_eles");
  AddBranch(&BJetness_num_medium_eles          ,"BJetness_num_medium_eles");
  AddBranch(&BJetness_num_tight_eles           ,"BJetness_num_tight_eles");
  AddBranch(&BJetness_num_mvatrig_eles         ,"BJetness_num_mvatrig_eles");
  AddBranch(&BJetness_num_mvanontrig_eles      ,"BJetness_num_mvanontrig_eles");
  AddBranch(&BJetness_num_mvatrigwp90_eles     ,"BJetness_num_mvatrigwp90_eles");
  AddBranch(&BJetness_num_mvanontrigwp90_eles  ,"BJetness_num_mvanontrigwp90_eles");
  AddBranch(&BJetness_num_heep_eles            ,"BJetness_num_heep_eles");
  AddBranch(&BJetness_num_pdgid_mus            ,"BJetness_num_pdgid_mus");
  AddBranch(&BJetness_num_loose_mus            ,"BJetness_num_loose_mus");
  AddBranch(&BJetness_num_soft_mus             ,"BJetness_num_soft_mus");
  AddBranch(&BJetness_num_medium_mus           ,"BJetness_num_medium_mus");
  AddBranch(&BJetness_num_tight_mus            ,"BJetness_num_tight_mus");
  AddBranch(&BJetness_num_highpt_mus           ,"BJetness_num_highpt_mus");
  AddBranch(&BJetness_num_POGisGood_mus        ,"BJetness_num_POGisGood_mus");
  AddBranch(&BJetness_numjettrks               ,"BJetness_numjettrks");
  AddBranch(&BJetness_numjettrkspv             ,"BJetness_numjettrkspv");
  AddBranch(&BJetness_numjettrksnopv           ,"BJetness_numjettrksnopv");
  AddBranch(&BJetness_npvTrkOVcollTrk          ,"BJetness_npvTrkOVcollTrk");
  AddBranch(&BJetness_pvTrkOVcollTrk           ,"BJetness_pvTrkOVcollTrk");
  AddBranch(&BJetness_npvTrkOVpvTrk            ,"BJetness_npvTrkOVpvTrk");
  AddBranch(&BJetness_npvPtOVcollPt            ,"BJetness_npvPtOVcollPt");
  AddBranch(&BJetness_pvPtOVcollPt             ,"BJetness_pvPtOVcollPt");
  AddBranch(&BJetness_npvPtOVpvPt              ,"BJetness_npvPtOVpvPt");
  //Trk prop rel to jet dir
  AddBranch(&BJetness_avprel                   ,"BJetness_avprel");
  AddBranch(&BJetness_avppar                   ,"BJetness_avppar");
  AddBranch(&BJetness_avetarel                 ,"BJetness_avetarel");
  AddBranch(&BJetness_avetapar                 ,"BJetness_avetapar");
  AddBranch(&BJetness_avdr                     ,"BJetness_avdr");
  AddBranch(&BJetness_avpreljetpt              ,"BJetness_avpreljetpt");
  AddBranch(&BJetness_avpreljeten              ,"BJetness_avpreljeten");
  AddBranch(&BJetness_avpparjetpt              ,"BJetness_avpparjetpt");
  AddBranch(&BJetness_avpparjeten              ,"BJetness_avpparjeten");
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
  AddBranch(&BJetness_numip3dpos               ,"BJetness_numip3dpos");
  AddBranch(&BJetness_numip3dneg               ,"BJetness_numip3dneg");
  AddBranch(&BJetness_avip2d_val               ,"BJetness_avip2d_val");
  AddBranch(&BJetness_avip2d_sig               ,"BJetness_avip2d_sig");
  AddBranch(&BJetness_avsip2d_val              ,"BJetness_avsip2d_val");
  AddBranch(&BJetness_avsip2d_sig              ,"BJetness_avsip2d_sig");
  AddBranch(&BJetness_numip2dpos               ,"BJetness_numip2dpos");
  AddBranch(&BJetness_numip2dneg               ,"BJetness_numip2dneg");
  AddBranch(&BJetness_avip1d_val               ,"BJetness_avip1d_val");
  AddBranch(&BJetness_avip1d_sig               ,"BJetness_avip1d_sig");
  AddBranch(&BJetness_avsip1d_val              ,"BJetness_avsip1d_val");
  AddBranch(&BJetness_avsip1d_sig              ,"BJetness_avsip1d_sig");
  if(debug_) std::cout<<"set branches"<<std::endl;
}
void BJetnessSelector::Clear(){
  //Evt selection
  BJetness_isSingleLepton = 0;
  BJetness_isDoubleLepton = 0;
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
  BJetness_pfCombinedMVAV2BJetTags.clear();
  BJetness_pfCombinedCvsLJetTags.clear();
  BJetness_pfCombinedCvsBJetTags.clear();
  BJetness_pt.clear();
  BJetness_eta.clear();
  BJetness_phi.clear();
  BJetness_en.clear();
  BJetness_ptOVen.clear();
  //PFCand info
  BJetness_jetschpvass.clear();
  BJetness_jetschfrompv.clear();
  BJetness_jetschip3dval.clear();
  BJetness_jetschip3dsig.clear();
  BJetness_jetschip2dval.clear();
  BJetness_jetschip2dsig.clear();
  BJetness_jetschisgoodtrk.clear();
  BJetness_jetschtrkpur.clear();
  BJetness_jetschpt.clear();
  BJetness_jetscheta.clear();
  BJetness_jetschen.clear();
  //Num_of_trks
  BJetness_num_pdgid_eles.clear();
  BJetness_num_soft_eles.clear();
  BJetness_num_vetonoipnoiso_eles.clear();
  BJetness_num_loosenoipnoiso_eles.clear();
  BJetness_num_veto_eles.clear();
  BJetness_num_loose_eles.clear();
  BJetness_num_medium_eles.clear();
  BJetness_num_tight_eles.clear();
  BJetness_num_mvatrig_eles.clear();
  BJetness_num_mvanontrig_eles.clear();
  BJetness_num_mvatrigwp90_eles.clear();
  BJetness_num_mvanontrigwp90_eles.clear();
  BJetness_num_heep_eles.clear();
  BJetness_num_pdgid_mus.clear();
  BJetness_num_loose_mus.clear();
  BJetness_num_soft_mus.clear();
  BJetness_num_medium_mus.clear();
  BJetness_num_tight_mus.clear();
  BJetness_num_highpt_mus.clear();
  BJetness_num_POGisGood_mus.clear();
  BJetness_numjettrks.clear();
  BJetness_numjettrkspv.clear();
  BJetness_numjettrksnopv.clear();
  BJetness_npvTrkOVcollTrk.clear();
  BJetness_pvTrkOVcollTrk.clear();
  BJetness_npvTrkOVpvTrk.clear();
  BJetness_npvPtOVcollPt.clear();
  BJetness_pvPtOVcollPt.clear();
  BJetness_npvPtOVpvPt.clear();
  //Trk prop rel to jet dir
  BJetness_avprel.clear();
  BJetness_avppar.clear();
  BJetness_avetarel.clear();
  BJetness_avetapar.clear();
  BJetness_avdr.clear();
  BJetness_avpreljetpt.clear();
  BJetness_avpreljeten.clear();
  BJetness_avpparjetpt.clear();
  BJetness_avpparjeten.clear();
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
  BJetness_numip3dpos.clear();
  BJetness_numip3dneg.clear();
  BJetness_avip2d_val.clear();
  BJetness_avip2d_sig.clear();
  BJetness_avsip2d_val.clear();
  BJetness_avsip2d_sig.clear();
  BJetness_numip2dpos.clear();
  BJetness_numip2dneg.clear();
  BJetness_avip1d_val.clear();
  BJetness_avip1d_sig.clear();
  BJetness_avsip1d_val.clear();
  BJetness_avsip1d_sig.clear();
}
//Ask for good vertices
bool BJetnessSelector::isGoodVertex(const reco::Vertex& vtx){
  if(vtx.isFake())                                   return false;
  if(vtx.ndof()<_vtx_ndof_min)                       return false;
  if(vtx.position().Rho()>_vtx_rho_max)              return false;
  if(fabs(vtx.position().Z()) > _vtx_position_z_max) return false;
  return true;
}

void BJetnessSelector::JECInitialization(){
  //AK4chs - MC: Get the factorized jet corrector parameters.
  std::vector<std::string> jecPayloadNamesAK4PFchsMC_;
  jecPayloadNamesAK4PFchsMC_.push_back(jecPayloadNamesAK4PFchsMC1_.fullPath());
  jecPayloadNamesAK4PFchsMC_.push_back(jecPayloadNamesAK4PFchsMC2_.fullPath());
  jecPayloadNamesAK4PFchsMC_.push_back(jecPayloadNamesAK4PFchsMC3_.fullPath());
  std::vector<JetCorrectorParameters> vParAK4PFchsMC;
  for ( std::vector<std::string>::const_iterator payloadBegin = jecPayloadNamesAK4PFchsMC_.begin(),
          payloadEnd = jecPayloadNamesAK4PFchsMC_.end(), ipayload = payloadBegin; ipayload != payloadEnd; ++ipayload ) {
    JetCorrectorParameters pars(*ipayload);
    vParAK4PFchsMC.push_back(pars);
  }
  jecAK4PFchsMC_    = boost::shared_ptr<FactorizedJetCorrector>  ( new FactorizedJetCorrector(vParAK4PFchsMC) );
  jecAK4PFchsMCUnc_ = boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(jecPayloadNamesAK4PFchsMCUnc_.fullPath()) );
  //AK4chs - DATA: Get the factorized jet corrector parameters.
  std::vector<std::string> jecPayloadNamesAK4PFchsDATA_;
  jecPayloadNamesAK4PFchsDATA_.push_back(jecPayloadNamesAK4PFchsDATA1_.fullPath());
  jecPayloadNamesAK4PFchsDATA_.push_back(jecPayloadNamesAK4PFchsDATA2_.fullPath());
  jecPayloadNamesAK4PFchsDATA_.push_back(jecPayloadNamesAK4PFchsDATA3_.fullPath());
  jecPayloadNamesAK4PFchsDATA_.push_back(jecPayloadNamesAK4PFchsDATA4_.fullPath());
  std::vector<JetCorrectorParameters> vParAK4PFchsDATA;
  for ( std::vector<std::string>::const_iterator payloadBegin = jecPayloadNamesAK4PFchsDATA_.begin(),
          payloadEnd = jecPayloadNamesAK4PFchsDATA_.end(), ipayload = payloadBegin; ipayload != payloadEnd; ++ipayload ) {
    JetCorrectorParameters pars(*ipayload);
    vParAK4PFchsDATA.push_back(pars);
  }
  jecAK4PFchsDATA_    = boost::shared_ptr<FactorizedJetCorrector>  ( new FactorizedJetCorrector(vParAK4PFchsDATA) );
  jecAK4PFchsDATAUnc_ = boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(jecPayloadNamesAK4PFchsDATAUnc_.fullPath()) );
}
void BJetnessSelector::GetJER(pat::Jet jet, float JesSF, float rhoJER, bool AK4PFchs, float &JERScaleFactor, float &JERScaleFactorUP, float &JERScaleFactorDOWN){
  if(!jet.genJet()) return;
  double jetEta=fabs(jet.eta());
  double cFactorJER = 1.0;
  double cFactorJERdown = 1.0;
  double cFactorJERup = 1.0;
  //https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution#JER_Scaling_factors_and_Uncertai
  if( jetEta<0.5 ){
    cFactorJER = 1.109;
    cFactorJERdown = 1.109-0.008;
    cFactorJERup   = 1.109+0.008;
  } else if( jetEta<0.8 ){
    cFactorJER = 1.138;
    cFactorJERdown = 1.138-0.013;
    cFactorJERup   = 1.138+0.013;
  } else if( jetEta<1.1 ){
    cFactorJER = 1.114;
    cFactorJERdown = 1.114-0.013;
    cFactorJERup   = 1.114+0.013;
  } else if( jetEta<1.3 ){
    cFactorJER = 1.123;
    cFactorJERdown = 1.123-0.024;
    cFactorJERup   = 1.123+0.024;
  } else if( jetEta<1.7 ){
    cFactorJER = 1.084;
    cFactorJERdown = 1.084-0.011;
    cFactorJERup   = 1.084+0.011;
  } else if( jetEta<1.9 ){
    cFactorJER = 1.082;
    cFactorJERdown = 1.082-0.035;
    cFactorJERup   = 1.082+0.035;
  } else if( jetEta<2.1 ){
    cFactorJER = 1.140;
    cFactorJERdown = 1.140-0.047;
    cFactorJERup   = 1.140+0.047;
  } else if( jetEta<2.3 ){
    cFactorJER = 1.067;
    cFactorJERdown = 1.067-0.053;
    cFactorJERup   = 1.067+0.053;
  } else if( jetEta<2.5 ){
    cFactorJER = 1.177;
    cFactorJERdown = 1.177-0.041;
    cFactorJERup   = 1.177+0.041;
  } else if( jetEta<2.8 ){
    cFactorJER = 1.364;
    cFactorJERdown = 1.364-0.039;
    cFactorJERup   = 1.364+0.039;
  } else if( jetEta<3.0 ){
    cFactorJER = 1.857;
    cFactorJERdown = 1.857-0.071;
    cFactorJERup   = 1.857+0.071;
  } else if( jetEta<3.2 ){
    cFactorJER = 0.328;
    cFactorJERdown = 0.328-0.022;
    cFactorJERup   = 0.328+0.022;
  } else if( jetEta<5.0 ){
    cFactorJER = 1.160;
    cFactorJERdown = 1.160-0.029;
    cFactorJERup   = 1.160+0.029;
  }
  //double recoJetPt = jet.pt();//(jet.correctedJet("Uncorrected").pt())*JesSF;
  double recoJetPt = (jet.correctedJet("Uncorrected").pt())*JesSF;
  double genJetPt  = jet.genJet()->pt();
  double diffPt    = recoJetPt - genJetPt;
  JME::JetResolution resolution;
  JME::JetResolutionScaleFactor res_sf;
  if(AK4PFchs){
    resolution = JME::JetResolution(jerAK4PFchs_);
    res_sf = JME::JetResolutionScaleFactor(jerAK4PFchsSF_);
  } else {
    //resolution = JME::JetResolution(jerAK4PFPuppi_);
    //res_sf = JME::JetResolutionScaleFactor(jerAK4PFPuppiSF_);
  }
  JME::JetParameters parameters;
  parameters.setJetPt(jet.pt());
  parameters.setJetEta(jet.eta());
  parameters.setRho(rhoJER);
  float relpterr = resolution.getResolution(parameters);
  if(genJetPt>0. && deltaR(jet.eta(),jet.phi(),jet.genJet()->eta(),jet.genJet()->phi())<0.2
     && (abs(jet.pt()-jet.genJet()->pt())<3*relpterr*jet.pt())) {
    JERScaleFactor     = (std::max(0., genJetPt + cFactorJER*diffPt))/recoJetPt;
    JERScaleFactorUP   = (std::max(0., genJetPt + cFactorJERup*diffPt))/recoJetPt;
    JERScaleFactorDOWN = (std::max(0., genJetPt + cFactorJERdown*diffPt))/recoJetPt;
  } else {
    JERScaleFactor     = 1.;
    JERScaleFactorUP   = 1.;
    JERScaleFactorDOWN = 1.;
  }
}


bool BJetnessSelector::is_softLep_jetelectron(const pat::Electron &lele){
  bool isele = false;
  const HitPattern &hitPattern = lele.gsfTrack().get()->hitPattern();
  uint32_t hit = hitPattern.getHitPattern(HitPattern::TRACK_HITS, 0);
  bool hitCondition = !(HitPattern::validHitFilter(hit) && ((HitPattern::pixelBarrelHitFilter(hit) && HitPattern::getLayer(hit) < 3) || HitPattern::pixelEndcapHitFilter(hit)));
  if(!hitCondition && lele.passConversionVeto()) isele = true;
  return isele;
}
bool BJetnessSelector::is_vetoPOGNoIPNoIso_jetelectron(const pat::Electron &lele){
  bool isele = false;
  double ooEmooP = 999;
  if(lele.ecalEnergy()==0)                   ooEmooP = 1e30;
  else if(!std::isfinite(lele.ecalEnergy())) ooEmooP = 1e30;
  else                                       ooEmooP = fabs(1.0/lele.ecalEnergy() - lele.eSuperClusterOverP()/lele.ecalEnergy());
  if(lele.full5x5_sigmaIetaIeta()<0.012 && fabs(lele.deltaEtaSuperClusterTrackAtVtx())<0.0126 && fabs(lele.deltaPhiSuperClusterTrackAtVtx())<0.107
     && lele.hcalOverEcal()<0.186 && ooEmooP<0.239
     //&& fabs((-1) * lele.gsfTrack()->dxy(vtx.position()))<0.0227 && fabs(lele.gsfTrack()->dz(vtx.position()))<0.379
     && lele.gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS)<=2 && lele.passConversionVeto()
    ) isele = true;
  return isele;
}
bool BJetnessSelector::is_loosePOGNoIPNoIso_jetelectron(const pat::Electron &lele){
  bool isele = false;
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
                                   const edm::Event& iEvent, edm::Handle<edm::View<pat::Electron> > electron_pat, edm::Handle<edm::View<pat::Muon> > muon_h,
                                   double& bjetness_num_pdgid_eles, double& bjetness_num_soft_eles, double& bjetness_num_vetonoipnoiso_eles, double& bjetness_num_loosenoipnoiso_eles, double& bjetness_num_veto_eles, double& bjetness_num_loose_eles, double& bjetness_num_medium_eles, double& bjetness_num_tight_eles, double& bjetness_num_mvatrig_eles, double& bjetness_num_mvanontrig_eles, double& bjetness_num_mvatrigwp90_eles, double& bjetness_num_mvanontrigwp90_eles, double& bjetness_num_heep_eles, double& bjetness_num_pdgid_mus, double& bjetness_num_loose_mus, double& bjetness_num_soft_mus, double& bjetness_num_medium_mus, double& bjetness_num_tight_mus, double& bjetness_num_highpt_mus, double& bjetness_num_POGisGood_mus
                                  ){
  //Access jet daughters
  vector<CandidatePtr> jdaus(jet.daughterPtrVector());
  sort(jdaus.begin(), jdaus.end(), [](const reco::CandidatePtr &p1, const reco::CandidatePtr &p2) {return p1->pt() > p2->pt();});
  for(uint jd=0; jd<jdaus.size(); jd++){
    const pat::PackedCandidate &jcand = dynamic_cast<const pat::PackedCandidate &>(*jdaus[jd]);
    if(deltaR(jcand.p4(),jet.p4())>0.4) continue;
    /////
    //   Trk information
    /////
    Track trk = Track(jcand.pseudoTrack());
    bool isgoodtrk = is_goodtrk(trk,vtx);
    //Minimal conditions for a track
    if((isgoodtrk || jcand.trackHighPurity()) && jcand.charge()!=0 && jcand.fromPV()>1 && jcand.vertexRef().key()==0){
      BJetness_jetschpvass.push_back(jcand.pvAssociationQuality());
      BJetness_jetschfrompv.push_back(jcand.fromPV());
      TransientTrack ttrk = ttrkbuilder.build(&trk);
      double valtemp = IPTools::absoluteImpactParameter3D(ttrk,vtx).second.value();
      if(valtemp==valtemp) BJetness_jetschip3dval.push_back(valtemp);
      else                 BJetness_jetschip3dval.push_back(0);
      valtemp = IPTools::absoluteImpactParameter3D(ttrk,vtx).second.significance();
      if(valtemp==valtemp) BJetness_jetschip3dsig.push_back(valtemp);
      else                 BJetness_jetschip3dsig.push_back(0);
      valtemp = IPTools::absoluteTransverseImpactParameter(ttrk,vtx).second.value();
      if(valtemp==valtemp) BJetness_jetschip2dval.push_back(valtemp);
      else                 BJetness_jetschip2dval.push_back(0);
      valtemp = IPTools::absoluteTransverseImpactParameter(ttrk,vtx).second.significance();
      if(valtemp==valtemp) BJetness_jetschip2dsig.push_back(valtemp);
      else                 BJetness_jetschip2dsig.push_back(0);
      BJetness_jetschisgoodtrk.push_back(isgoodtrk);
      BJetness_jetschtrkpur.push_back(jcand.trackHighPurity());
      BJetness_jetschpt.push_back(jcand.pt());
      BJetness_jetscheta.push_back(jcand.eta());
      BJetness_jetschen.push_back(jcand.energy());
      //Trk info
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
    /////
    //   Leptons information
    /////
    if(jcand.charge()!=0){
      //ele
      if(fabs(jcand.pdgId())==11) bjetness_num_pdgid_eles++;
      edm::Handle<edm::ValueMap<bool>  > veto_id_decisions;
      edm::Handle<edm::ValueMap<bool>  > loose_id_decisions;
      edm::Handle<edm::ValueMap<bool>  > medium_id_decisions;
      edm::Handle<edm::ValueMap<bool>  > tight_id_decisions;
      edm::Handle<edm::ValueMap<bool>  > heep_id_decisions;
      edm::Handle<edm::ValueMap<bool>  > mvatrig_id_decisions;
      edm::Handle<edm::ValueMap<bool>  > mvanontrig_id_decisions;
      edm::Handle<edm::ValueMap<bool>  > mvatrigwp90_id_decisions;
      edm::Handle<edm::ValueMap<bool>  > mvanontrigwp90_id_decisions;
      iEvent.getByToken(electronVetoIdMapToken_,   veto_id_decisions);
      iEvent.getByToken(electronLooseIdMapToken_,  loose_id_decisions);
      iEvent.getByToken(electronMediumIdMapToken_, medium_id_decisions);
      iEvent.getByToken(electronTightIdMapToken_,  tight_id_decisions);
      iEvent.getByToken(eleMVATrigIdMapToken_,     mvatrig_id_decisions);
      iEvent.getByToken(eleMVAnonTrigIdMap_,       mvanontrig_id_decisions);
      iEvent.getByToken(eleMVATrigwp90IdMap_,      mvatrigwp90_id_decisions);
      iEvent.getByToken(eleMVAnonTrigwp90IdMap_,   mvanontrigwp90_id_decisions);
      iEvent.getByToken(eleHEEPIdMapToken_,        heep_id_decisions);
      for(edm::View<pat::Electron>::const_iterator ele = electron_pat->begin(); ele != electron_pat->end(); ele++){
        const pat::Electron &lele = *ele;
        if(deltaR(jcand.p4(),lele.p4())<0.1 && fabs(jcand.pt()-lele.pt())/lele.pt()<0.05){
          const Ptr<pat::Electron> elPtr(electron_pat, ele - electron_pat->begin() );
          if(is_softLep_jetelectron(lele)) bjetness_num_soft_eles++;
          if(is_vetoPOGNoIPNoIso_jetelectron(lele)) bjetness_num_vetonoipnoiso_eles++;
          if(is_loosePOGNoIPNoIso_jetelectron(lele)) bjetness_num_loosenoipnoiso_eles++;
          if((*veto_id_decisions)  [ elPtr ]) bjetness_num_veto_eles++;
          if((*loose_id_decisions) [ elPtr ]) bjetness_num_loose_eles++;
          if((*medium_id_decisions)[ elPtr ]) bjetness_num_medium_eles++;
          if((*tight_id_decisions) [ elPtr ]) bjetness_num_tight_eles++;
          if((*mvatrig_id_decisions) [ elPtr ]) bjetness_num_mvatrig_eles++;
          if((*mvanontrig_id_decisions) [ elPtr ]) bjetness_num_mvanontrig_eles++;
          if((*mvatrigwp90_id_decisions) [ elPtr ]) bjetness_num_mvatrigwp90_eles++;
          if((*mvanontrigwp90_id_decisions) [ elPtr ]) bjetness_num_mvanontrigwp90_eles++;
          if((*heep_id_decisions)  [ elPtr ]) bjetness_num_heep_eles++;
        }
      }
      //mu
      if(fabs(jcand.pdgId())==13) bjetness_num_pdgid_mus++;
      for(const pat::Muon &mu : *muon_h){
        if(deltaR(jcand.p4(),mu.p4())<0.1 && fabs(jcand.pt()-mu.pt())/mu.pt()<0.05){
          if(mu.isLooseMuon()) bjetness_num_loose_mus++;
          if(mu.isSoftMuon(vtx)) bjetness_num_soft_mus++;
          if(mu.isMediumMuon()) bjetness_num_medium_mus++;
          if(mu.isTightMuon(vtx)) bjetness_num_tight_mus++;
          if(mu.isHighPtMuon(vtx)) bjetness_num_highpt_mus++;
          if(muon::isGoodMuon(mu, muon::TMOneStationTight)) bjetness_num_POGisGood_mus++;
        }
      }
    }
  }//Loop on jet daus
}
//Get trk kin properties
void BJetnessSelector::get_avreljet(vector<Track> trks, vector<tuple<double, double, double> >& jetsdir, double& jetchtrks_avprel, double& jetchtrks_avppar, double& jetchtrks_avetarel, double& jetchtrks_avetapar, double& jetchtrks_avdr){
  for(uint t=0; t<trks.size(); t++){
    TVector3 trkdir(trks[t].px(),trks[t].py(),trks[t].pz());
    TVector3 axis(get<0>(jetsdir[t]),get<1>(jetsdir[t]),get<2>(jetsdir[t]));
    jetchtrks_avprel += trkdir.Perp(axis);//Projection of track vector perpendicular to jet "axis"
    jetchtrks_avppar += trkdir.Dot(axis.Unit());//Projection of track vector parallel to jet "axis"
    double energy = std::sqrt(trkdir.Mag2() + pow(0.13957,2));
    if((energy - trkdir.Perp(axis))!=0 && (energy + trkdir.Perp(axis))!=0) jetchtrks_avetarel += 0.5 * log((energy + trkdir.Perp(axis)) / (energy - trkdir.Perp(axis)));
    if((energy - trkdir.Dot(axis.Unit()))!=0 && (energy + trkdir.Dot(axis.Unit()))!=0) jetchtrks_avetapar += 0.5 * log((energy + trkdir.Dot(axis.Unit())) / (energy - trkdir.Dot(axis.Unit())));
    jetchtrks_avdr += trkdir.DeltaR(axis);
  }
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
void BJetnessSelector::get_avip3d(vector<Track> trks, const TransientTrackBuilder& ttrkbuilder, reco::Vertex vtx, vector<tuple<double, double, double> >& jetsdir, double& jetchtrks_avip3d_val, double& jetchtrks_avip3d_sig, double& jetchtrks_avsip3d_val, double& jetchtrks_avsip3d_sig, double& jetchtrks_numip3dpos, double& jetchtrks_numip3dneg){
  double valtemp = 0;
  vector<TransientTrack> ttrks = get_ttrks(trks,ttrkbuilder);
  for(uint t=0; t<ttrks.size(); t++){
    valtemp = IPTools::absoluteImpactParameter3D(ttrks[t],vtx).second.value();
    if(valtemp==valtemp) jetchtrks_avip3d_val  += valtemp;
    valtemp = IPTools::absoluteImpactParameter3D(ttrks[t],vtx).second.significance();
    if(valtemp==valtemp) jetchtrks_avip3d_sig  += valtemp;
    GlobalVector jetsdirgv(get<0>(jetsdir[t]),get<1>(jetsdir[t]),get<2>(jetsdir[t]));
    valtemp = IPTools::signedImpactParameter3D(ttrks[t],jetsdirgv,vtx).second.value();
    if(valtemp==valtemp){
      jetchtrks_avsip3d_val += valtemp;
      if(valtemp!=0){
        if(valtemp/fabs(valtemp)>=0) jetchtrks_numip3dpos++;
        if(valtemp/fabs(valtemp)<0)  jetchtrks_numip3dneg++;
      }
    }
    valtemp = IPTools::signedImpactParameter3D(ttrks[t],jetsdirgv,vtx).second.significance();
    if(valtemp==valtemp) jetchtrks_avsip3d_sig += valtemp;
  }
}
void BJetnessSelector::get_avip2d(vector<Track> trks, const TransientTrackBuilder& ttrkbuilder, reco::Vertex vtx, vector<tuple<double, double, double> >& jetsdir, double& jetchtrks_avip2d_val, double& jetchtrks_avip2d_sig, double& jetchtrks_avsip2d_val, double& jetchtrks_avsip2d_sig, double& jetchtrks_numip2dpos, double& jetchtrks_numip2dneg){
  double valtemp = 0;
  vector<TransientTrack> ttrks = get_ttrks(trks,ttrkbuilder);
  for(uint t=0; t<ttrks.size(); t++){
    valtemp = IPTools::absoluteTransverseImpactParameter(ttrks[t],vtx).second.value();
    if(valtemp==valtemp) jetchtrks_avip2d_val  += valtemp;
    valtemp = IPTools::absoluteTransverseImpactParameter(ttrks[t],vtx).second.significance();
    if(valtemp==valtemp) jetchtrks_avip2d_sig  += valtemp;
    GlobalVector jetsdirgv(get<0>(jetsdir[t]),get<1>(jetsdir[t]),get<2>(jetsdir[t]));
    valtemp = IPTools::signedTransverseImpactParameter(ttrks[t],jetsdirgv,vtx).second.value();
    if(valtemp==valtemp){
      jetchtrks_avsip2d_val += valtemp;
      if(valtemp!=0){
        if(valtemp/fabs(valtemp)>=0) jetchtrks_numip2dpos++;
        if(valtemp/fabs(valtemp)<0)  jetchtrks_numip2dneg++;
      }
    }
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
