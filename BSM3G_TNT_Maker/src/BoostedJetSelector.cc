#include "BSMFramework/BSM3G_TNT_Maker/interface/BoostedJetSelector.h"
BoostedJetSelector::BoostedJetSelector(std::string name, TTree* tree, bool debug, const pset& iConfig):baseTree(name,tree,debug){
  fatjetToken_ = iConfig.getParameter<edm::InputTag>("fatjets");
  _vertexInputTag = iConfig.getParameter<edm::InputTag>("vertices");
  jecPayloadNamesAK8PFchsMC1_   = iConfig.getParameter<edm::FileInPath>("jecPayloadNamesAK8PFchsMC1");
  jecPayloadNamesAK8PFchsMC2_   = iConfig.getParameter<edm::FileInPath>("jecPayloadNamesAK8PFchsMC2");
  jecPayloadNamesAK8PFchsMC3_   = iConfig.getParameter<edm::FileInPath>("jecPayloadNamesAK8PFchsMC3");
  jecPayloadNamesAK8PFchsMCUnc_ = iConfig.getParameter<edm::FileInPath>("jecPayloadNamesAK8PFchsMCUnc");
  jecPayloadNamesAK8PFchsDATA1_   = iConfig.getParameter<edm::FileInPath>("jecPayloadNamesAK8PFchsDATA1");
  jecPayloadNamesAK8PFchsDATA2_   = iConfig.getParameter<edm::FileInPath>("jecPayloadNamesAK8PFchsDATA2");
  jecPayloadNamesAK8PFchsDATA3_   = iConfig.getParameter<edm::FileInPath>("jecPayloadNamesAK8PFchsDATA3");
  jecPayloadNamesAK8PFchsDATAUnc_ = iConfig.getParameter<edm::FileInPath>("jecPayloadNamesAK8PFchsDATAUnc");
  _is_data = iConfig.getParameter<bool>("is_data");
  JECInitialization();
  SetBranches();
}
BoostedJetSelector::~BoostedJetSelector(){
  delete tree_;
}
void BoostedJetSelector::Fill(const edm::Event& iEvent){
  Clear();
  /////
  //   Recall collections
  /////  
  edm::Handle<pat::JetCollection> fatjets;                                       
  iEvent.getByLabel(fatjetToken_, fatjets); 
  edm::Handle<double> rhoHandle;
  iEvent.getByLabel("fixedGridRhoFastjetAll",rhoHandle);
  double rho = *rhoHandle;
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByLabel(_vertexInputTag, vertices);                                        
  /////
  //   Get fatjet information
  /////  
  for(const pat::Jet &j : *fatjets){ 
    //if (j.pt() < _Jet_pt_min) continue;
    //Kinematic
    BoostedJet_pt.push_back(j.pt());         
    BoostedJet_eta.push_back(j.eta());       
    BoostedJet_phi.push_back(j.phi());       
    BoostedJet_energy.push_back(j.energy());
    BoostedJet_mass.push_back(j.mass()); 
    BoostedJet_Uncorr_pt.push_back(j.correctedJet("Uncorrected").pt());    
    //ID
    BoostedJet_combinedSecondaryVertexBJetTags.push_back(j.bDiscriminator("combinedSecondaryVertexBJetTags"));              
    BoostedJet_pfCombinedSecondaryVertexV2BJetTags.push_back(j.bDiscriminator("pfCombinedSecondaryVertexV2BJetTags"));              
    BoostedJet_pfCombinedInclusiveSecondaryVertexV2BJetTags.push_back(j.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));              
    //Energy related variables
    BoostedJet_neutralHadEnergyFraction.push_back(j.neutralHadronEnergyFraction());                               
    BoostedJet_neutralEmEmEnergyFraction.push_back(j.neutralEmEnergyFraction());                                   
    BoostedJet_chargedHadronEnergyFraction.push_back(j.chargedHadronEnergyFraction());                               
    BoostedJet_chargedEmEnergyFraction.push_back(j.chargedEmEnergyFraction());                              
    BoostedJet_muonEnergyFraction.push_back(j.muonEnergyFraction());                                  
    BoostedJet_numberOfConstituents.push_back(j.chargedMultiplicity() + j.neutralMultiplicity());                                  
    BoostedJet_chargedMultiplicity.push_back(j.chargedMultiplicity());
    BoostedJet_electronEnergy.push_back(j.electronEnergy());                               
    BoostedJet_photonEnergy.push_back(j.photonEnergy());
    //Boosted jet prop 
    BoostedJet_tau1.push_back(j.userFloat("NjettinessAK8:tau1"));    //
    BoostedJet_tau2.push_back(j.userFloat("NjettinessAK8:tau2"));    //  Access the n-subjettiness variables
    BoostedJet_tau3.push_back(j.userFloat("NjettinessAK8:tau3"));    // 
    BoostedJet_softdrop_mass.push_back(j.userFloat("ak8PFJetsCHSSoftDropMass")); // access to filtered mass
    BoostedJet_trimmed_mass.push_back(j.userFloat("ak8PFJetsCHSTrimmedMass"));   // access to trimmed mass
    BoostedJet_pruned_mass.push_back(j.userFloat("ak8PFJetsCHSPrunedMass"));     // access to pruned mass
    BoostedJet_filtered_mass.push_back(j.userFloat("ak8PFJetsCHSFilteredMass")); // access to filtered mass
    //Jet Energy Corrections and Uncertainties
    double corrAK8PFchs     = 1;
    double corrUpAK8PFchs   = 1;
    double corrDownAK8PFchs = 1;
    reco::Candidate::LorentzVector uncorrJetAK8PFchs = j.correctedP4(0);
    if(!_is_data){
      jecAK8PFchsMC_->setJetEta( uncorrJetAK8PFchs.eta()    );
      jecAK8PFchsMC_->setJetPt ( uncorrJetAK8PFchs.pt()     );
      jecAK8PFchsMC_->setJetE  ( uncorrJetAK8PFchs.energy() );
      jecAK8PFchsMC_->setRho	( rho  );
      jecAK8PFchsMC_->setNPV	( vertices->size()  );
      jecAK8PFchsMC_->setJetA  ( j.jetArea()	     );
      corrAK8PFchs = jecAK8PFchsMC_->getCorrection();
      jecAK8PFchsMCUnc_->setJetEta( uncorrJetAK8PFchs.eta() );
      jecAK8PFchsMCUnc_->setJetPt( corrAK8PFchs * uncorrJetAK8PFchs.pt() );
      corrUpAK8PFchs = corrAK8PFchs * (1 + fabs(jecAK8PFchsMCUnc_->getUncertainty(1)));
      jecAK8PFchsMCUnc_->setJetEta( uncorrJetAK8PFchs.eta() );
      jecAK8PFchsMCUnc_->setJetPt( corrAK8PFchs * uncorrJetAK8PFchs.pt() );
      corrDownAK8PFchs = corrAK8PFchs * ( 1 - fabs(jecAK8PFchsMCUnc_->getUncertainty(-1)) );
    } else {
      jecAK8PFchsDATA_->setJetEta( uncorrJetAK8PFchs.eta()    );
      jecAK8PFchsDATA_->setJetPt ( uncorrJetAK8PFchs.pt()     );
      jecAK8PFchsDATA_->setJetE  ( uncorrJetAK8PFchs.energy() );
      jecAK8PFchsDATA_->setRho	( rho  );
      jecAK8PFchsDATA_->setNPV	( vertices->size()  );
      jecAK8PFchsDATA_->setJetA  ( j.jetArea()	     );
      corrAK8PFchs = jecAK8PFchsDATA_->getCorrection();
      jecAK8PFchsDATAUnc_->setJetEta( uncorrJetAK8PFchs.eta() );
      jecAK8PFchsDATAUnc_->setJetPt( corrAK8PFchs * uncorrJetAK8PFchs.pt() );
      corrUpAK8PFchs = corrAK8PFchs * (1 + fabs(jecAK8PFchsDATAUnc_->getUncertainty(1)));
      jecAK8PFchsDATAUnc_->setJetEta( uncorrJetAK8PFchs.eta() );
      jecAK8PFchsDATAUnc_->setJetPt( corrAK8PFchs * uncorrJetAK8PFchs.pt() );
      corrDownAK8PFchs = corrAK8PFchs * ( 1 - fabs(jecAK8PFchsDATAUnc_->getUncertainty(-1)) );
    }
    BoostedJet_JesSF.push_back(corrAK8PFchs);
    BoostedJet_JesSFup.push_back(corrUpAK8PFchs);
    BoostedJet_JesSFdown.push_back(corrDownAK8PFchs);
    //Variables for top-tagging
    double TopMass = -10.;
    double MinMass = -10.;
    double WMass = -10.;
    int NSubJets = -10;
    reco::CATopJetTagInfo const * tagInfo =  dynamic_cast<reco::CATopJetTagInfo const *>( j.tagInfo("caTop"));
    if( tagInfo != 0 ){
      TopMass  = tagInfo->properties().topMass;
      MinMass  = tagInfo->properties().minMass;
      WMass    = tagInfo->properties().wMass;
      NSubJets = tagInfo->properties().nSubJets;
    }
    TopTagging_topMass.push_back(TopMass);
    TopTagging_minMass.push_back(MinMass);
    TopTagging_wMass.push_back(WMass);
    TopTagging_nSubJets.push_back(NSubJets);
  } 
}
void BoostedJetSelector::JECInitialization(){
  //AK8chs - MC: Get the factorized jet corrector parameters. 
  std::vector<std::string> jecPayloadNamesAK8PFchsMC_;
  jecPayloadNamesAK8PFchsMC_.push_back(jecPayloadNamesAK8PFchsMC1_.fullPath());
  jecPayloadNamesAK8PFchsMC_.push_back(jecPayloadNamesAK8PFchsMC2_.fullPath());
  jecPayloadNamesAK8PFchsMC_.push_back(jecPayloadNamesAK8PFchsMC3_.fullPath());
  std::vector<JetCorrectorParameters> vParAK8PFchsMC;
  for ( std::vector<std::string>::const_iterator payloadBegin = jecPayloadNamesAK8PFchsMC_.begin(),
	  payloadEnd = jecPayloadNamesAK8PFchsMC_.end(), ipayload = payloadBegin; ipayload != payloadEnd; ++ipayload ) {
    JetCorrectorParameters pars(*ipayload);
    vParAK8PFchsMC.push_back(pars);
  }
  jecAK8PFchsMC_    = boost::shared_ptr<FactorizedJetCorrector>  ( new FactorizedJetCorrector(vParAK8PFchsMC) );
  jecAK8PFchsMCUnc_ = boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(jecPayloadNamesAK8PFchsMCUnc_.fullPath()) );
  //AK8chs - DATA: Get the factorized jet corrector parameters. 
  std::vector<std::string> jecPayloadNamesAK8PFchsDATA_;
  jecPayloadNamesAK8PFchsDATA_.push_back(jecPayloadNamesAK8PFchsDATA1_.fullPath());
  jecPayloadNamesAK8PFchsDATA_.push_back(jecPayloadNamesAK8PFchsDATA2_.fullPath());
  jecPayloadNamesAK8PFchsDATA_.push_back(jecPayloadNamesAK8PFchsDATA3_.fullPath());
  std::vector<JetCorrectorParameters> vParAK8PFchsDATA;
  for ( std::vector<std::string>::const_iterator payloadBegin = jecPayloadNamesAK8PFchsDATA_.begin(),
	  payloadEnd = jecPayloadNamesAK8PFchsDATA_.end(), ipayload = payloadBegin; ipayload != payloadEnd; ++ipayload ) {
    JetCorrectorParameters pars(*ipayload);
    vParAK8PFchsDATA.push_back(pars);
  }
  jecAK8PFchsDATA_    = boost::shared_ptr<FactorizedJetCorrector>  ( new FactorizedJetCorrector(vParAK8PFchsDATA) );
  jecAK8PFchsDATAUnc_ = boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(jecPayloadNamesAK8PFchsDATAUnc_.fullPath()) );
}
void BoostedJetSelector::SetBranches(){
  if(debug_)    std::cout<<"setting branches: calling AddBranch of baseTree"<<std::endl;
  //Kinematic
  AddBranch(&BoostedJet_pt,                          "BoostedJet_pt");
  AddBranch(&BoostedJet_eta,                         "BoostedJet_eta");
  AddBranch(&BoostedJet_phi,                         "BoostedJet_phi");
  AddBranch(&BoostedJet_energy,                      "BoostedJet_energy");
  AddBranch(&BoostedJet_mass,                        "BoostedJet_mass");
  AddBranch(&BoostedJet_Uncorr_pt ,                  "BoostedJet_Uncorr_pt");
  //ID
  AddBranch(&BoostedJet_combinedSecondaryVertexBJetTags,              "BoostedJet_combinedSecondaryVertexBJetTags");
  AddBranch(&BoostedJet_pfCombinedInclusiveSecondaryVertexV2BJetTags, "BoostedJet_pfCombinedSecondaryVertexV2BJetTags");
  AddBranch(&BoostedJet_pfCombinedInclusiveSecondaryVertexV2BJetTags, "BoostedJet_pfCombinedInclusiveSecondaryVertexV2BJetTags");
  //Energy related variables
  AddBranch(&BoostedJet_neutralHadEnergyFraction,    "BoostedJet_neutralHadEnergyFraction");
  AddBranch(&BoostedJet_neutralEmEmEnergyFraction,   "BoostedJet_neutralEmEmEnergyFraction");
  AddBranch(&BoostedJet_chargedHadronEnergyFraction, "BoostedJet_chargedHadronEnergyFraction");
  AddBranch(&BoostedJet_chargedEmEnergyFraction,     "BoostedJet_chargedEmEnergyFraction");
  AddBranch(&BoostedJet_muonEnergyFraction,          "BoostedJet_muonEnergyFraction");
  AddBranch(&BoostedJet_numberOfConstituents,        "BoostedJet_numberOfConstituents");
  AddBranch(&BoostedJet_chargedMultiplicity,         "BoostedJet_chargedMultiplicity");
  AddBranch(&BoostedJet_electronEnergy,              "BoostedJet_electronEnergy");
  AddBranch(&BoostedJet_photonEnergy,                "BoostedJet_photonEnergy");
  //Boosted jet prop 
  AddBranch(&BoostedJet_tau1,           "BoostedJet_tau1");
  AddBranch(&BoostedJet_tau2,           "BoostedJet_tau2");
  AddBranch(&BoostedJet_tau3,           "BoostedJet_tau3");
  AddBranch(&BoostedJet_softdrop_mass,  "BoostedJet_softdrop_mass");
  AddBranch(&BoostedJet_trimmed_mass,   "BoostedJet_trimmed_mass");
  AddBranch(&BoostedJet_pruned_mass,    "BoostedJet_pruned_mass");
  AddBranch(&BoostedJet_filtered_mass,  "BoostedJet_filtered_mass");
  //Jet Energy Corrections and Uncertainties
  AddBranch(&BoostedJet_JesSF                ,"BoostedJet_JesSF");
  AddBranch(&BoostedJet_JesSFup              ,"BoostedJet_JesSFup");
  AddBranch(&BoostedJet_JesSFdown            ,"BoostedJet_JesSFdown");
  //Variables for top-tagging
  AddBranch(&TopTagging_topMass,  "TopTagging_topMass");
  AddBranch(&TopTagging_minMass,  "TopTagging_minMass");
  AddBranch(&TopTagging_wMass,    "TopTagging_wMass");
  AddBranch(&TopTagging_nSubJets, "TopTagging_nSubJets");
  if(debug_)    std::cout<<"set branches"<<std::endl;
}
void BoostedJetSelector::Clear(){
  //Kinematic
  BoostedJet_pt.clear();
  BoostedJet_eta.clear();
  BoostedJet_phi.clear();
  BoostedJet_energy.clear();
  BoostedJet_mass.clear();
  BoostedJet_Uncorr_pt.clear();
  //ID
  BoostedJet_combinedSecondaryVertexBJetTags.clear();
  BoostedJet_pfCombinedSecondaryVertexV2BJetTags.clear();
  BoostedJet_pfCombinedInclusiveSecondaryVertexV2BJetTags.clear();
  //Energy related variables
  BoostedJet_neutralHadEnergyFraction.clear();
  BoostedJet_neutralEmEmEnergyFraction.clear();
  BoostedJet_chargedHadronEnergyFraction.clear();
  BoostedJet_chargedEmEnergyFraction.clear();
  BoostedJet_muonEnergyFraction.clear();
  BoostedJet_numberOfConstituents.clear();
  BoostedJet_chargedMultiplicity.clear();
  BoostedJet_electronEnergy.clear();
  BoostedJet_photonEnergy.clear();
  //Boosted jet prop 
  BoostedJet_tau1.clear();
  BoostedJet_tau2.clear();
  BoostedJet_tau3.clear();
  BoostedJet_softdrop_mass.clear();
  BoostedJet_trimmed_mass.clear();
  BoostedJet_pruned_mass.clear();
  BoostedJet_filtered_mass.clear();
  //Jet Energy Corrections and Uncertainties
  BoostedJet_JesSF.clear();
  BoostedJet_JesSFup.clear();
  BoostedJet_JesSFdown.clear();
  //Variables for top-tagging
  TopTagging_topMass.clear();
  TopTagging_minMass.clear();
  TopTagging_wMass.clear();
  TopTagging_nSubJets.clear();
}
