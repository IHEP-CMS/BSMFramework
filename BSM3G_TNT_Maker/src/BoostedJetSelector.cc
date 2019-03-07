#include "BSMFramework/BSM3G_TNT_Maker/interface/BoostedJetSelector.h"
BoostedJetSelector::BoostedJetSelector(std::string name, TTree* tree, bool debug, const pset& iConfig, edm::ConsumesCollector && ic):
  baseTree(name,tree,debug)
{
  vtx_h_        = ic.consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"));
  fatjets_      = ic.consumes<pat::JetCollection >(iConfig.getParameter<edm::InputTag>("fatjets"));
  rhopogHandle_ = ic.consumes<double>(edm::InputTag("fixedGridRhoFastjetAll"));
  rhoJERHandle_ = ic.consumes<double>(edm::InputTag("fixedGridRhoAll"));
  jecPayloadNamesAK8PFchsMC1_   = iConfig.getParameter<edm::FileInPath>("jecPayloadNamesAK8PFchsMC1");
  jecPayloadNamesAK8PFchsMC2_   = iConfig.getParameter<edm::FileInPath>("jecPayloadNamesAK8PFchsMC2");
  jecPayloadNamesAK8PFchsMC3_   = iConfig.getParameter<edm::FileInPath>("jecPayloadNamesAK8PFchsMC3");
  jecPayloadNamesAK8PFchsMCUnc_ = iConfig.getParameter<edm::FileInPath>("jecPayloadNamesAK8PFchsMCUnc");
  jecPayloadNamesAK8PFchsDATA1_   = iConfig.getParameter<edm::FileInPath>("jecPayloadNamesAK8PFchsDATA1");
  jecPayloadNamesAK8PFchsDATA2_   = iConfig.getParameter<edm::FileInPath>("jecPayloadNamesAK8PFchsDATA2");
  jecPayloadNamesAK8PFchsDATA3_   = iConfig.getParameter<edm::FileInPath>("jecPayloadNamesAK8PFchsDATA3");
  jecPayloadNamesAK8PFchsDATA4_   = iConfig.getParameter<edm::FileInPath>("jecPayloadNamesAK8PFchsDATA4");
  jecPayloadNamesAK8PFchsDATAUnc_ = iConfig.getParameter<edm::FileInPath>("jecPayloadNamesAK8PFchsDATAUnc");
  jerAK8PFchs_     = iConfig.getParameter<edm::FileInPath>("jerAK8PFchs").fullPath();
  jerAK8PFchsSF_   = iConfig.getParameter<edm::FileInPath>("jerAK8PFchsSF").fullPath();
  jerAK8PFPuppi_   = iConfig.getParameter<edm::FileInPath>("jerAK8PFPuppi").fullPath();
  jerAK8PFPuppiSF_ = iConfig.getParameter<edm::FileInPath>("jerAK8PFPuppiSF").fullPath();
  _is_data = iConfig.getParameter<bool>("is_data");
  _reduced = iConfig.getParameter<bool>("reduced");
  _MC2016  = iConfig.getParameter<bool>("MC2016");
  PuppiWeightFilePath_ = iConfig.getParameter<edm::FileInPath>("PuppiWeightFilePath");
  const char *filePath = PuppiWeightFilePath_.fullPath().c_str();
  PuppiWeightFile = new TFile (filePath,"READ");
  puppisd_corrGEN      = (TF1*)PuppiWeightFile->Get("puppiJECcorr_gen");
  puppisd_corrRECO_cen = (TF1*)PuppiWeightFile->Get("puppiJECcorr_reco_0eta1v3");
  puppisd_corrRECO_for = (TF1*)PuppiWeightFile->Get("puppiJECcorr_reco_1v3eta2v5");
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
  edm::Handle<double> rhoJERHandle;
  iEvent.getByToken(rhoJERHandle_,rhoJERHandle);
  double rhoJER = *rhoJERHandle;
  /////
  //   Get fatjet information
  /////  
  for(const pat::Jet &j : *fatjets){ 
    if (j.pt() <= 170.) continue;
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
    BoostedJet_neutralEmEnergyFraction.push_back(j.neutralEmEnergyFraction());                                   
    BoostedJet_chargedHadronEnergyFraction.push_back(j.chargedHadronEnergyFraction());                               
    BoostedJet_chargedEmEnergyFraction.push_back(j.chargedEmEnergyFraction());                              
    BoostedJet_muonEnergyFraction.push_back(j.muonEnergyFraction());                                  
    BoostedJet_numberOfConstituents.push_back(j.chargedMultiplicity() + j.neutralMultiplicity());                                  
    BoostedJet_chargedMultiplicity.push_back(j.chargedMultiplicity());
    BoostedJet_electronEnergy.push_back(j.electronEnergy());                               
    BoostedJet_photonEnergy.push_back(j.photonEnergy());
    //Boosted jet prop
    BoostedJet_tau1.push_back(j.userFloat("ak8PFJetsCHSValueMap:NjettinessAK8CHSTau1"));    //
    BoostedJet_tau2.push_back(j.userFloat("ak8PFJetsCHSValueMap:NjettinessAK8CHSTau2"));    //  Access the n-subjettiness variables
    BoostedJet_tau3.push_back(j.userFloat("ak8PFJetsCHSValueMap:NjettinessAK8CHSTau3"));    // 
    BoostedJet_softdrop_mass.push_back(j.userFloat("ak8PFJetsCHSValueMap:ak8PFJetsCHSSoftDropMass")); // access to soft drop mass
    BoostedJet_pruned_mass.push_back(j.userFloat("ak8PFJetsCHSValueMap:ak8PFJetsCHSPrunedMass"));     // access to pruned mass
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
    if(!_is_data) GetJER(j, corrAK8PFchs, rhoJER, true, JERScaleFactor, JERScaleFactorUP, JERScaleFactorDOWN);
    BoostedJet_JerSF.push_back(JERScaleFactor);
    BoostedJet_JerSFup.push_back(JERScaleFactorUP);
    BoostedJet_JerSFdown.push_back(JERScaleFactorDOWN);
    //PUPPI Softdrop
    BoostedJet_puppi_softdrop_mass.push_back(j.userFloat("ak8PFJetsPuppiSoftDropMass"));
    BoostedJet_puppi_tau1.push_back(j.userFloat("NjettinessAK8Puppi:tau1"));
    BoostedJet_puppi_tau2.push_back(j.userFloat("NjettinessAK8Puppi:tau2"));
    BoostedJet_puppi_tau3.push_back(j.userFloat("NjettinessAK8Puppi:tau3"));
    BoostedJet_puppi_pt.push_back(-99);//j.userFloat("ak8PFJetsPuppiValueMap:pt"));
    BoostedJet_puppi_eta.push_back(-99);//j.userFloat("ak8PFJetsPuppiValueMap:eta"));
    BoostedJet_puppi_phi.push_back(-99);//j.userFloat("ak8PFJetsPuppiValueMap:phi"));
    TLorentzVector puppi_softdrop, puppi_softdrop_subjet;
    auto const & sdSubjetsPuppi = j.subjets("SoftDropPuppi");
    for ( auto const & it : sdSubjetsPuppi ) {
      puppi_softdrop_subjet.SetPtEtaPhiM(it->correctedP4(0).pt(),it->correctedP4(0).eta(),it->correctedP4(0).phi(),it->correctedP4(0).mass());
      puppi_softdrop+=puppi_softdrop_subjet;
    }
    float puppiCorr = getPUPPIweight( j.correctedJet("Uncorrected").pt()*corrAK8PFchs*JERScaleFactor ,j.eta() );
    BoostedJet_puppi_softdrop_masscorr.push_back(puppi_softdrop.M()*puppiCorr);
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
  if(!_reduced){
    AddBranch(&BoostedJet_pfJetProbabilityBJetTags,                     "BoostedJet_pfJetProbabilityBJetTags");
    AddBranch(&BoostedJet_pfCombinedMVAV2BJetTags,                      "BoostedJet_pfCombinedMVAV2BJetTags");
    AddBranch(&BoostedJet_pfCombinedInclusiveSecondaryVertexV2BJetTags, "BoostedJet_pfCombinedInclusiveSecondaryVertexV2BJetTags");
    AddBranch(&BoostedJet_pfCombinedCvsLJetTags                        ,"BoostedJet_pfCombinedCvsLJetTags");
    AddBranch(&BoostedJet_pfCombinedCvsBJetTags                        ,"BoostedJet_pfCombinedCvsBJetTags");
  }
  //Energy related variables
  AddBranch(&BoostedJet_neutralHadEnergyFraction,    "BoostedJet_neutralHadEnergyFraction");
  AddBranch(&BoostedJet_neutralEmEnergyFraction,     "BoostedJet_neutralEmEnergyFraction");
  AddBranch(&BoostedJet_chargedHadronEnergyFraction, "BoostedJet_chargedHadronEnergyFraction");
  AddBranch(&BoostedJet_chargedEmEnergyFraction,     "BoostedJet_chargedEmEnergyFraction");
  AddBranch(&BoostedJet_muonEnergyFraction,          "BoostedJet_muonEnergyFraction");
  AddBranch(&BoostedJet_numberOfConstituents,        "BoostedJet_numberOfConstituents");
  AddBranch(&BoostedJet_chargedMultiplicity,         "BoostedJet_chargedMultiplicity");
  AddBranch(&BoostedJet_electronEnergy,              "BoostedJet_electronEnergy");
  AddBranch(&BoostedJet_photonEnergy,                "BoostedJet_photonEnergy");
  //Boosted jet prop 
  AddBranch(&BoostedJet_tau1,               "BoostedJet_tau1");
  AddBranch(&BoostedJet_tau2,               "BoostedJet_tau2");
  AddBranch(&BoostedJet_tau3,               "BoostedJet_tau3");
  AddBranch(&BoostedJet_softdrop_mass,      "BoostedJet_softdrop_mass");
  AddBranch(&BoostedJet_pruned_mass,        "BoostedJet_pruned_mass");
  AddBranch(&BoostedJet_puppi_softdrop_mass,"BoostedJet_puppi_softdrop_mass");
  AddBranch(&BoostedJet_puppi_tau1,         "BoostedJet_puppi_tau1");
  AddBranch(&BoostedJet_puppi_tau2,         "BoostedJet_puppi_tau2");
  AddBranch(&BoostedJet_puppi_tau3,         "BoostedJet_puppi_tau3");
  if(!_reduced){
    AddBranch(&BoostedJet_puppi_pt,  "BoostedJet_puppi_pt");
    AddBranch(&BoostedJet_puppi_eta, "BoostedJet_puppi_eta");
    AddBranch(&BoostedJet_puppi_phi, "BoostedJet_puppi_phi");
    AddBranch(&BoostedJet_puppi_softdrop_masscorr,"BoostedJet_puppi_softdrop_masscorr");
  }
  //Jet Energy Corrections and Uncertainties
  AddBranch(&BoostedJet_JesSF                ,"BoostedJet_JesSF");
  AddBranch(&BoostedJet_JesSFup              ,"BoostedJet_JesSFup");
  AddBranch(&BoostedJet_JesSFdown            ,"BoostedJet_JesSFdown");
  AddBranch(&BoostedJet_JerSF                ,"BoostedJet_JerSF");
  AddBranch(&BoostedJet_JerSFup              ,"BoostedJet_JerSFup");
  AddBranch(&BoostedJet_JerSFdown            ,"BoostedJet_JerSFdown");
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
  BoostedJet_neutralEmEnergyFraction.clear();
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
  BoostedJet_pruned_mass.clear();
  BoostedJet_puppi_pt.clear();
  BoostedJet_puppi_softdrop_mass.clear();
  BoostedJet_puppi_eta.clear();
  BoostedJet_puppi_phi.clear();
  BoostedJet_puppi_tau1.clear();
  BoostedJet_puppi_tau2.clear();
  BoostedJet_puppi_tau3.clear();
  BoostedJet_puppi_softdrop_masscorr.clear();
  //Jet Energy Corrections and Uncertainties
  BoostedJet_JesSF.clear();
  BoostedJet_JesSFup.clear();
  BoostedJet_JesSFdown.clear();
  BoostedJet_JerSF.clear();
  BoostedJet_JerSFup.clear();
  BoostedJet_JerSFdown.clear();
}
void BoostedJetSelector::GetJER(pat::Jet jet, float JesSF, float rhoJER, bool AK8PFchs, float &JERScaleFactor, float &JERScaleFactorUP, float &JERScaleFactorDOWN){
  if(!jet.genJet()) return;
  double jetEta=fabs(jet.eta());
  double cFactorJER = 1.0; 
  double cFactorJERdown = 1.0;
  double cFactorJERup = 1.0;
  //https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution#JER_Scaling_factors_and_Unce_AN1
  if(!_MC2016){
    if( jetEta<0.522 ){ 
      cFactorJER = 1.1432; 
      cFactorJERdown = 1.1432-0.0222;
      cFactorJERup   = 1.1432+0.0222; 
    } else if( jetEta<0.783 ){ 
      cFactorJER = 1.1815; 
      cFactorJERdown = 1.1815-0.0484;
      cFactorJERup   = 1.1815+0.0484; 
    } else if( jetEta<1.131 ){ 
      cFactorJER = 1.0989; 
      cFactorJERdown = 1.0989-0.0456;
      cFactorJERup   = 1.0989+0.0456; 
    } else if( jetEta<1.305 ){ 
      cFactorJER = 1.1137; 
      cFactorJERdown = 1.1137-0.1397;
      cFactorJERup   = 1.1137+0.1397; 
    } else if( jetEta<1.740 ){ 
      cFactorJER = 1.1307; 
      cFactorJERdown = 1.1307-0.1470;
      cFactorJERup   = 1.1307+0.1470; 
    } else if( jetEta<1.930 ){ 
      cFactorJER = 1.1600; 
      cFactorJERdown = 1.1600-0.0976;
      cFactorJERup   = 1.1600+0.0976; 
    } else if( jetEta<2.043 ){ 
      cFactorJER = 1.2393; 
      cFactorJERdown = 1.2393-0.1909;
      cFactorJERup   = 1.2393+0.1909; 
    } else if( jetEta<2.322 ){ 
      cFactorJER = 1.2604; 
      cFactorJERdown = 1.2604-0.1501;
      cFactorJERup   = 1.2604+0.1501; 
    } else if( jetEta<2.500 ){ 
      cFactorJER = 1.4085; 
      cFactorJERdown = 1.4085-0.2020;
      cFactorJERup   = 1.4085+0.2020; 
    } else if( jetEta<2.853 ){ 
      cFactorJER = 1.9909; 
      cFactorJERdown = 1.9909-0.5684;
      cFactorJERup   = 1.9909+0.5684; 
    } else if( jetEta<2.964 ){ 
      cFactorJER = 2.2923; 
      cFactorJERdown = 2.2923-0.3743;
      cFactorJERup   = 2.2923+0.3743; 
    } else if( jetEta<3.139 ){ 
      cFactorJER = 1.2696; 
      cFactorJERdown = 1.2696-0.1089;
      cFactorJERup   = 1.2696+0.1089; 
    } else if( jetEta<5.191 ){ 
      cFactorJER = 1.1542; 
      cFactorJERdown = 1.1542-0.1524;
      cFactorJERup   = 1.1542+0.1524;
    }
  } else {
    if( jetEta<0.522 ){ 
      cFactorJER = 1.1595; 
      cFactorJERdown = 1.1595-0.0645;
      cFactorJERup   = 1.1595+0.0645; 
    } else if( jetEta<0.783 ){ 
      cFactorJER = 1.1948; 
      cFactorJERdown = 1.1948-0.0652;
      cFactorJERup   = 1.1948+0.0652; 
    } else if( jetEta<1.131 ){ 
      cFactorJER = 1.1464; 
      cFactorJERdown = 1.1464-0.0632;
      cFactorJERup   = 1.1464+0.0632; 
    } else if( jetEta<1.305 ){ 
      cFactorJER = 1.1609; 
      cFactorJERdown = 1.1609-0.1025;
      cFactorJERup   = 1.1609+0.1025; 
    } else if( jetEta<1.740 ){ 
      cFactorJER = 1.1278; 
      cFactorJERdown = 1.1278-0.0986;
      cFactorJERup   = 1.1278+0.0986; 
    } else if( jetEta<1.930 ){ 
      cFactorJER = 1.1000; 
      cFactorJERdown = 1.1000-0.1079;
      cFactorJERup   = 1.1000+0.1079; 
    } else if( jetEta<2.043 ){ 
      cFactorJER = 1.1426; 
      cFactorJERdown = 1.1426-0.1214;
      cFactorJERup   = 1.1426+0.1214; 
    } else if( jetEta<2.322 ){ 
      cFactorJER = 1.1512; 
      cFactorJERdown = 1.1512-0.1140;
      cFactorJERup   = 1.1512+0.1140; 
    } else if( jetEta<2.500 ){ 
      cFactorJER = 1.2963; 
      cFactorJERdown = 1.2963-0.2371;
      cFactorJERup   = 1.2963+0.2371; 
    } else if( jetEta<2.853 ){ 
      cFactorJER = 1.3418; 
      cFactorJERdown = 1.3418-0.2091;
      cFactorJERup   = 1.3418+0.2091; 
    } else if( jetEta<2.964 ){ 
      cFactorJER = 1.7788; 
      cFactorJERdown = 1.7788-0.2008;
      cFactorJERup   = 1.7788+0.2008; 
    } else if( jetEta<3.139 ){ 
      cFactorJER = 1.1869; 
      cFactorJERdown = 1.1869-0.1243;
      cFactorJERup   = 1.1869+0.1243; 
    } else if( jetEta<5.191 ){ 
      cFactorJER = 1.1922; 
      cFactorJERdown = 1.1922-0.1488;
      cFactorJERup   = 1.1922+0.1488;
    }
  }
  //double recoJetPt = jet.pt();//(jet.correctedJet("Uncorrected").pt())*JesSF;
  double recoJetPt = (jet.correctedJet("Uncorrected").pt())*JesSF;
  double genJetPt  = jet.genJet()->pt();
  double diffPt    = recoJetPt - genJetPt;
  JME::JetResolution resolution;
  JME::JetResolutionScaleFactor res_sf;
  if(AK8PFchs){
    resolution = JME::JetResolution(jerAK8PFchs_);
    res_sf = JME::JetResolutionScaleFactor(jerAK8PFchsSF_);
  } else {
    resolution = JME::JetResolution(jerAK8PFPuppi_);
    res_sf = JME::JetResolutionScaleFactor(jerAK8PFPuppiSF_);
  }
  JME::JetParameters parameters;
  parameters.setJetPt(jet.pt());
  parameters.setJetEta(jet.eta());
  parameters.setRho(rhoJER);
  float relpterr = resolution.getResolution(parameters);
  if(genJetPt>0. && deltaR(jet.eta(),jet.phi(),jet.genJet()->eta(),jet.genJet()->phi())<0.2 && (abs(diffPt)<3*relpterr*recoJetPt)){
    JERScaleFactor     = (std::max(0., genJetPt + cFactorJER*diffPt))/recoJetPt;
    JERScaleFactorUP   = (std::max(0., genJetPt + cFactorJERup*diffPt))/recoJetPt;
    JERScaleFactorDOWN = (std::max(0., genJetPt + cFactorJERdown*diffPt))/recoJetPt;
  } else {
    JERScaleFactor     = 1.;
    JERScaleFactorUP   = 1.;
    JERScaleFactorDOWN = 1.;
  }
}
float BoostedJetSelector::getPUPPIweight(float puppipt, float puppieta ){
  float genCorr  = 1.;
  float recoCorr = 1.;
  float totalWeight = 1.;
  genCorr =  puppisd_corrGEN->Eval( puppipt );
  if( fabs(puppieta)  <= 1.3 ){
    recoCorr = puppisd_corrRECO_cen->Eval( puppipt );
  }
  else{
    recoCorr = puppisd_corrRECO_for->Eval( puppipt );
  }
  totalWeight = genCorr * recoCorr;
  return totalWeight;
}
