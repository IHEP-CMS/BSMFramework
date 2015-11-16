#include "BSMFramework/BSM3G_TNT_Maker/interface/TopSubJetSelector.h"
TopSubJetSelector::TopSubJetSelector(std::string name, TTree* tree, bool debug, const pset& iConfig):baseTree(name,tree,debug){
  topsubjetToken_ = iConfig.getParameter<edm::InputTag>("topsubjets");
  _vertexInputTag = iConfig.getParameter<edm::InputTag>("vertices");
  jecPayloadNamesAK4PFchsMC1_   = iConfig.getParameter<edm::FileInPath>("jecPayloadNamesAK4PFchsMC1");
  jecPayloadNamesAK4PFchsMC2_   = iConfig.getParameter<edm::FileInPath>("jecPayloadNamesAK4PFchsMC2");
  jecPayloadNamesAK4PFchsMC3_   = iConfig.getParameter<edm::FileInPath>("jecPayloadNamesAK4PFchsMC3");
  jecPayloadNamesAK4PFchsMCUnc_ = iConfig.getParameter<edm::FileInPath>("jecPayloadNamesAK4PFchsMCUnc");
  jecPayloadNamesAK4PFchsDATA1_   = iConfig.getParameter<edm::FileInPath>("jecPayloadNamesAK4PFchsDATA1");
  jecPayloadNamesAK4PFchsDATA2_   = iConfig.getParameter<edm::FileInPath>("jecPayloadNamesAK4PFchsDATA2");
  jecPayloadNamesAK4PFchsDATA3_   = iConfig.getParameter<edm::FileInPath>("jecPayloadNamesAK4PFchsDATA3");
  jecPayloadNamesAK4PFchsDATAUnc_ = iConfig.getParameter<edm::FileInPath>("jecPayloadNamesAK4PFchsDATAUnc");
  _Jet_pt_min     = iConfig.getParameter<double>("Jet_pt_min");
  _is_data = iConfig.getParameter<bool>("is_data");
  JECInitialization();
  SetBranches();
}
TopSubJetSelector::~TopSubJetSelector(){
  delete tree_;
}
void TopSubJetSelector::Fill(const edm::Event& iEvent){
  Clear();
  /////
  //   Require a good vertex 
  /////  
  edm::Handle<pat::JetCollection> topsubjets;                                       
  iEvent.getByLabel(topsubjetToken_, topsubjets);
  edm::Handle<double> rhoHandle;
  iEvent.getByLabel("fixedGridRhoFastjetAll",rhoHandle);
  double rho = *rhoHandle;
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByLabel(_vertexInputTag, vertices);                                         
  /////
  //   Get TopSubJet information
  /////  
  for(const pat::Jet &j : *topsubjets){ 
    //Acceptance
    if(j.pt()<_Jet_pt_min) continue;
    //Kinematics
    TopSubjet_pt.push_back(j.pt());         
    TopSubjet_eta.push_back(j.eta());       
    TopSubjet_phi.push_back(j.phi());       
    TopSubjet_energy.push_back(j.energy());
    TopSubjet_mass.push_back(j.mass()); 
    //ID
    TopSubjet_Btag0.push_back(j.bDiscriminator("combinedSecondaryVertexBJetTags"));
    TopSubjet_Btag1.push_back(j.bDiscriminator("pfCombinedSecondaryVertexV2BJetTags"));
    TopSubjet_Btag2.push_back(j.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
    //Jet Energy Corrections and Uncertainties
    double corrAK4PFchs     = 1;
    double corrUpAK4PFchs   = 1;
    double corrDownAK4PFchs = 1;
    reco::Candidate::LorentzVector uncorrJetAK4PFchs = j.correctedP4(0);
    if(!_is_data){
      jecAK4PFchsMC_->setJetEta( uncorrJetAK4PFchs.eta()    );
      jecAK4PFchsMC_->setJetPt ( uncorrJetAK4PFchs.pt()     );
      jecAK4PFchsMC_->setJetE  ( uncorrJetAK4PFchs.energy() );
      jecAK4PFchsMC_->setRho	( rho  );
      jecAK4PFchsMC_->setNPV	( vertices->size()  );
      jecAK4PFchsMC_->setJetA  ( j.jetArea()	     );
      corrAK4PFchs = jecAK4PFchsMC_->getCorrection();
      jecAK4PFchsMCUnc_->setJetEta( uncorrJetAK4PFchs.eta() );
      jecAK4PFchsMCUnc_->setJetPt( corrAK4PFchs * uncorrJetAK4PFchs.pt() );
      corrUpAK4PFchs = corrAK4PFchs * (1 + fabs(jecAK4PFchsMCUnc_->getUncertainty(1)));
      jecAK4PFchsMCUnc_->setJetEta( uncorrJetAK4PFchs.eta() );
      jecAK4PFchsMCUnc_->setJetPt( corrAK4PFchs * uncorrJetAK4PFchs.pt() );
      corrDownAK4PFchs = corrAK4PFchs * ( 1 - fabs(jecAK4PFchsMCUnc_->getUncertainty(-1)) );
    } else {
      jecAK4PFchsDATA_->setJetEta( uncorrJetAK4PFchs.eta()    );
      jecAK4PFchsDATA_->setJetPt ( uncorrJetAK4PFchs.pt()     );
      jecAK4PFchsDATA_->setJetE  ( uncorrJetAK4PFchs.energy() );
      jecAK4PFchsDATA_->setRho	( rho  );
      jecAK4PFchsDATA_->setNPV	( vertices->size()  );
      jecAK4PFchsDATA_->setJetA  ( j.jetArea()	     );
      corrAK4PFchs = jecAK4PFchsDATA_->getCorrection();
      jecAK4PFchsDATAUnc_->setJetEta( uncorrJetAK4PFchs.eta() );
      jecAK4PFchsDATAUnc_->setJetPt( corrAK4PFchs * uncorrJetAK4PFchs.pt() );
      corrUpAK4PFchs = corrAK4PFchs * (1 + fabs(jecAK4PFchsDATAUnc_->getUncertainty(1)));
      jecAK4PFchsDATAUnc_->setJetEta( uncorrJetAK4PFchs.eta() );
      jecAK4PFchsDATAUnc_->setJetPt( corrAK4PFchs * uncorrJetAK4PFchs.pt() );
      corrDownAK4PFchs = corrAK4PFchs * ( 1 - fabs(jecAK4PFchsDATAUnc_->getUncertainty(-1)) );
    }
    TopSubjet_JesSF.push_back(corrAK4PFchs);
    TopSubjet_JesSFup.push_back(corrUpAK4PFchs);
    TopSubjet_JesSFdown.push_back(corrDownAK4PFchs);
  } 
}
void TopSubJetSelector::JECInitialization(){
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
  std::vector<JetCorrectorParameters> vParAK4PFchsDATA;
  for ( std::vector<std::string>::const_iterator payloadBegin = jecPayloadNamesAK4PFchsDATA_.begin(),
	  payloadEnd = jecPayloadNamesAK4PFchsDATA_.end(), ipayload = payloadBegin; ipayload != payloadEnd; ++ipayload ) {
    JetCorrectorParameters pars(*ipayload);
    vParAK4PFchsDATA.push_back(pars);
  }
  jecAK4PFchsDATA_    = boost::shared_ptr<FactorizedJetCorrector>  ( new FactorizedJetCorrector(vParAK4PFchsDATA) );
  jecAK4PFchsDATAUnc_ = boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(jecPayloadNamesAK4PFchsDATAUnc_.fullPath()) );
}
void TopSubJetSelector::SetBranches(){
  if(debug_) std::cout<<"setting branches: calling AddBranch of baseTree"<<std::endl;
  //Kinematics
  AddBranch(&TopSubjet_pt,       "TopSubjet_pt");
  AddBranch(&TopSubjet_eta,      "TopSubjet_eta");
  AddBranch(&TopSubjet_phi,      "TopSubjet_phi");
  AddBranch(&TopSubjet_energy,   "TopSubjet_energy");
  AddBranch(&TopSubjet_mass,     "TopSubjet_mass");
  //ID
  AddBranch(&TopSubjet_Btag0,    "TopSubjet_Btag0");
  AddBranch(&TopSubjet_Btag1,    "TopSubjet_Btag1");
  AddBranch(&TopSubjet_Btag2,    "TopSubjet_Btag2");
  //Jet Energy Corrections and Uncertainties
  AddBranch(&TopSubjet_JesSF    ,"TopSubjet_JesSF");
  AddBranch(&TopSubjet_JesSFup  ,"TopSubjet_JesSFup");
  AddBranch(&TopSubjet_JesSFdown,"TopSubjet_JesSFdown");
  if(debug_) std::cout<<"set branches"<<std::endl;
}
void TopSubJetSelector::Clear(){
  //Kinematics
  TopSubjet_pt.clear();
  TopSubjet_eta.clear();
  TopSubjet_phi.clear();
  TopSubjet_energy.clear();
  TopSubjet_mass.clear();
  //ID
  TopSubjet_Btag0.clear();
  TopSubjet_Btag1.clear();
  TopSubjet_Btag2.clear();
  //Jet Energy Corrections and Uncertainties
  TopSubjet_JesSF.clear();
  TopSubjet_JesSFup.clear();
  TopSubjet_JesSFdown.clear();
}
