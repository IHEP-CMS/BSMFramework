#include "BSMFramework/BSM3G_TNT_Maker/interface/BoostedJetSelector.h"
BoostedJetSelector::BoostedJetSelector(std::string name, TTree* tree, bool debug, const pset& iConfig, edm::ConsumesCollector && ic):
  baseTree(name,tree,debug)
{
  vtx_h_        = ic.consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"));
  fatjets_      = ic.consumes<pat::JetCollection >(iConfig.getParameter<edm::InputTag>("fatjets"));
  rhopogHandle_ = ic.consumes<double>(edm::InputTag("fixedGridRhoFastjetAll"));
  jecPayloadNamesAK8PFchsMC1_   = iConfig.getParameter<edm::FileInPath>("jecPayloadNamesAK8PFchsMC1");
  jecPayloadNamesAK8PFchsMC2_   = iConfig.getParameter<edm::FileInPath>("jecPayloadNamesAK8PFchsMC2");
  jecPayloadNamesAK8PFchsMC3_   = iConfig.getParameter<edm::FileInPath>("jecPayloadNamesAK8PFchsMC3");
  jecPayloadNamesAK8PFchsMCUnc_ = iConfig.getParameter<edm::FileInPath>("jecPayloadNamesAK8PFchsMCUnc");
  jecPayloadNamesAK8PFchsDATA1_   = iConfig.getParameter<edm::FileInPath>("jecPayloadNamesAK8PFchsDATA1");
  jecPayloadNamesAK8PFchsDATA2_   = iConfig.getParameter<edm::FileInPath>("jecPayloadNamesAK8PFchsDATA2");
  jecPayloadNamesAK8PFchsDATA3_   = iConfig.getParameter<edm::FileInPath>("jecPayloadNamesAK8PFchsDATA3");
  jecPayloadNamesAK8PFchsDATA4_   = iConfig.getParameter<edm::FileInPath>("jecPayloadNamesAK8PFchsDATA4");
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
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vtx_h_, vertices);                                        
  edm::Handle<pat::JetCollection> fatjets;                                       
  iEvent.getByToken(fatjets_, fatjets); 
  edm::Handle<double> rhoHandle;
  iEvent.getByToken(rhopogHandle_,rhoHandle);
  double rho = *rhoHandle;
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
    BoostedJet_pfJetProbabilityBJetTags.push_back(j.bDiscriminator("pfJetProbabilityBJetTags"));              
    BoostedJet_pfCombinedMVAV2BJetTags.push_back(j.bDiscriminator("pfCombinedMVAV2BJetTags"));              
    BoostedJet_pfCombinedInclusiveSecondaryVertexV2BJetTags.push_back(j.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));              
    BoostedJet_pfCombinedCvsLJetTags.push_back(j.bDiscriminator("pfCombinedCvsLJetTags"));
    BoostedJet_pfCombinedCvsBJetTags.push_back(j.bDiscriminator("pfCombinedCvsBJetTags"));
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
    BoostedJet_softdrop_mass.push_back(j.userFloat("ak8PFJetsCHSSoftDropMass")); // access to soft drop mass
    //BoostedJet_trimmed_mass.push_back(j.userFloat("ak8PFJetsCHSTrimmedMass"));   // access to trimmed mass
    BoostedJet_pruned_mass.push_back(j.userFloat("ak8PFJetsCHSPrunedMass"));     // access to pruned mass
    //BoostedJet_filtered_mass.push_back(j.userFloat("ak8PFJetsCHSFilteredMass")); // access to filtered mass
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
    //JER scale factor and uncertainties
    float JERScaleFactor     = 1; 
    float JERScaleFactorUP   = 1;
    float JERScaleFactorDOWN = 1;
    if(!_is_data) GetJER(j, corrAK8PFchs, JERScaleFactor, JERScaleFactorUP, JERScaleFactorDOWN);
    BoostedJet_JerSF.push_back(JERScaleFactor);
    BoostedJet_JerSFup.push_back(JERScaleFactorUP);
    BoostedJet_JerSFdown.push_back(JERScaleFactorDOWN);
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
  jecPayloadNamesAK8PFchsDATA_.push_back(jecPayloadNamesAK8PFchsDATA4_.fullPath());
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
  AddBranch(&BoostedJet_pfJetProbabilityBJetTags,                     "BoostedJet_pfJetProbabilityBJetTags");
  AddBranch(&BoostedJet_pfCombinedMVAV2BJetTags,                      "BoostedJet_pfCombinedMVAV2BJetTags");
  AddBranch(&BoostedJet_pfCombinedInclusiveSecondaryVertexV2BJetTags, "BoostedJet_pfCombinedInclusiveSecondaryVertexV2BJetTags");
  AddBranch(&BoostedJet_pfCombinedCvsLJetTags                        ,"BoostedJet_pfCombinedCvsLJetTags");
  AddBranch(&BoostedJet_pfCombinedCvsBJetTags                        ,"BoostedJet_pfCombinedCvsBJetTags");
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
  //AddBranch(&BoostedJet_trimmed_mass,   "BoostedJet_trimmed_mass");
  AddBranch(&BoostedJet_pruned_mass,    "BoostedJet_pruned_mass");
  //AddBranch(&BoostedJet_filtered_mass,  "BoostedJet_filtered_mass");
  //Jet Energy Corrections and Uncertainties
  AddBranch(&BoostedJet_JesSF                ,"BoostedJet_JesSF");
  AddBranch(&BoostedJet_JesSFup              ,"BoostedJet_JesSFup");
  AddBranch(&BoostedJet_JesSFdown            ,"BoostedJet_JesSFdown");
  AddBranch(&BoostedJet_JerSF                ,"BoostedJet_JerSF");
  AddBranch(&BoostedJet_JerSFup              ,"BoostedJet_JerSFup");
  AddBranch(&BoostedJet_JerSFdown            ,"BoostedJet_JerSFdown");
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
  BoostedJet_pfJetProbabilityBJetTags.clear();
  BoostedJet_pfCombinedMVAV2BJetTags.clear();
  BoostedJet_pfCombinedInclusiveSecondaryVertexV2BJetTags.clear();
  BoostedJet_pfCombinedCvsLJetTags.clear();
  BoostedJet_pfCombinedCvsBJetTags.clear();
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
  //BoostedJet_trimmed_mass.clear();
  BoostedJet_pruned_mass.clear();
  //BoostedJet_filtered_mass.clear();
  //Jet Energy Corrections and Uncertainties
  BoostedJet_JesSF.clear();
  BoostedJet_JesSFup.clear();
  BoostedJet_JesSFdown.clear();
  BoostedJet_JerSF.clear();
  BoostedJet_JerSFup.clear();
  BoostedJet_JerSFdown.clear();
  //Variables for top-tagging
  TopTagging_topMass.clear();
  TopTagging_minMass.clear();
  TopTagging_wMass.clear();
  TopTagging_nSubJets.clear();
}
void BoostedJetSelector::GetJER(pat::Jet jet, float JesSF, float &JERScaleFactor, float &JERScaleFactorUP, float &JERScaleFactorDOWN){
  if(!jet.genJet()) return;
  double jetEta=fabs(jet.eta());
  double cFactorJER = 1.0; 
  double cFactorJERdown = 1.0;
  double cFactorJERup = 1.0;
  //https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution#JER_Scaling_factors_and_Unce_AN1
  if( jetEta<0.5 ){ 
    cFactorJER = 1.095; 
    cFactorJERdown = 1.095-0.018;
    cFactorJERup   = 1.095+0.018; 
  } else if( jetEta<0.8 ){ 
    cFactorJER = 1.120; 
    cFactorJERdown = 1.120-0.028;
    cFactorJERup   = 1.120+0.028; 
  } else if( jetEta<1.1 ){ 
    cFactorJER = 1.097; 
    cFactorJERdown = 1.097-0.017;
    cFactorJERup   = 1.097+0.017; 
  } else if( jetEta<1.3 ){ 
    cFactorJER = 1.103; 
    cFactorJERdown = 1.103-0.033;
    cFactorJERup   = 1.103+0.033; 
  } else if( jetEta<1.7 ){ 
    cFactorJER = 1.118; 
    cFactorJERdown = 1.118-0.014;
    cFactorJERup   = 1.118+0.014; 
  } else if( jetEta<1.9 ){ 
    cFactorJER = 1.100; 
    cFactorJERdown = 1.100-0.033;
    cFactorJERup   = 1.100+0.033; 
  } else if( jetEta<2.1 ){ 
    cFactorJER = 1.162; 
    cFactorJERdown = 1.162-0.044;
    cFactorJERup   = 1.162+0.044; 
  } else if( jetEta<2.3 ){ 
    cFactorJER = 1.160; 
    cFactorJERdown = 1.160-0.048;
    cFactorJERup   = 1.160+0.048; 
  } else if( jetEta<2.5 ){ 
    cFactorJER = 1.161; 
    cFactorJERdown = 1.161-0.060;
    cFactorJERup   = 1.161+0.060; 
  } else if( jetEta<2.8 ){ 
    cFactorJER = 1.209; 
    cFactorJERdown = 1.209-0.059;
    cFactorJERup   = 1.209+0.059; 
  } else if( jetEta<3.0 ){ 
    cFactorJER = 1.564; 
    cFactorJERdown = 1.564-0.321;
    cFactorJERup   = 1.564+0.321; 
  } else if( jetEta<3.2 ){ 
    cFactorJER = 1.384; 
    cFactorJERdown = 1.384-0.033;
    cFactorJERup   = 1.384+0.033; 
  } else if( jetEta<5.0 ){ 
    cFactorJER = 1.216; 
    cFactorJERdown = 1.216-0.050;
    cFactorJERup   = 1.216+0.050;
  }
  double recoJetPt = (jet.correctedJet("Uncorrected").pt())*JesSF;
  double genJetPt  = jet.genJet()->pt();
  double diffPt    = recoJetPt - genJetPt;
  if(genJetPt>0.){
    JERScaleFactor     = (std::max(0., genJetPt + cFactorJER*diffPt))/recoJetPt;
    JERScaleFactorUP   = (std::max(0., genJetPt + cFactorJERup*diffPt))/recoJetPt;
    JERScaleFactorDOWN = (std::max(0., genJetPt + cFactorJERdown*diffPt))/recoJetPt;
  } else {
    JERScaleFactor     = 1.;
    JERScaleFactorUP   = 1.;
    JERScaleFactorDOWN = 1.;
  } 
}
