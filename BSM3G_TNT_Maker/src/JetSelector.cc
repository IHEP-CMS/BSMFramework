#include "BSMFramework/BSM3G_TNT_Maker/interface/JetSelector.h"
JetSelector::JetSelector(std::string name, TTree* tree, bool debug, const pset& iConfig):baseTree(name,tree,debug){
  jetToken_       = iConfig.getParameter<edm::InputTag>("jets");
  _vertexInputTag = iConfig.getParameter<edm::InputTag>("vertices");
  _Jet_pt_min     = iConfig.getParameter<double>("Jet_pt_min");
  _super_TNT      = iConfig.getParameter<bool>("super_TNT");
  SetBranches();
}
JetSelector::~JetSelector(){
  delete tree_;
}
void JetSelector::Fill(const edm::Event& iEvent){
  Clear();
  JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty("../files/Summer13_V5_MC_Uncertainty_AK5PFchs.txt");
  /////
  //   Recall collections
  /////  
  edm::Handle<pat::JetCollection> jets;                                       
  iEvent.getByLabel(jetToken_, jets);                                         
  /////
  //   Get jet information
  /////  
  for(const pat::Jet &j : *jets){ 
    //Acceptance
    if(j.pt()<_Jet_pt_min) continue;
    //Kinematics
    Jet_pt.push_back(j.pt());  
    Jet_px.push_back(j.px());   
    Jet_py.push_back(j.py());          
    Jet_eta.push_back(j.eta());       
    Jet_phi.push_back(j.phi());       
    Jet_energy.push_back(j.energy());
    Jet_mass.push_back(j.mass()); 
    UncorrJet_pt.push_back(j.correctedJet("Uncorrected").pt());                
    //ID
    Jet_bDiscriminator.push_back(j.bDiscriminator("combinedSecondaryVertexBJetTags"));
    Jet_bDiscriminator1.push_back(j.bDiscriminator("combinedInclusiveSecondaryVertexV2BJetTags"));
    Jet_bDiscriminator2.push_back(j.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
    Jet_pileupId.push_back(j.userFloat("pileupJetId:fullDiscriminant"));
    //Jet Uncertainties
    float JesUncertainties=0;
    GetJESUncertainties(j, jecUnc, JesUncertainties);
    Jet_JesUp.push_back((1+JesUncertainties));         
    Jet_JesDown.push_back((1-JesUncertainties));
    //JER scale factor and uncertainties
    float JERScaleFactor    =1; 
    float JERScaleFactorUP  =1;
    float JERScaleFactorDOWN=1;
    GetJER(j, JERScaleFactor, JERScaleFactorUP, JERScaleFactorDOWN);
    Jet_JerSF.push_back(JERScaleFactor);
    Jet_JerSFup.push_back(JERScaleFactorUP);
    Jet_JerSFdown.push_back(JERScaleFactorDOWN);
    //Energy related variables
    if(!_super_TNT){
      //Jet_neutralHadEnergy.push_back(j.neutralHadronEnergy());                               
      Jet_neutralHadEnergyFraction.push_back(j.neutralHadronEnergyFraction());                               
      //Jet_neutralEmEmEnergy.push_back(j.neutralEmEnergy());                                   
      Jet_neutralEmEmEnergyFraction.push_back(j.neutralEmEnergyFraction());                                   
      //Jet_chargedHadronEnergy.push_back(j.chargedHadronEnergy());                               
      Jet_chargedHadronEnergyFraction.push_back(j.chargedHadronEnergyFraction());                               
      //Jet_chargedEmEnergy.push_back(j.chargedEmEnergy());                              
      Jet_chargedEmEnergyFraction.push_back(j.chargedEmEnergyFraction());                              
      //Jet_muonEnergy.push_back(j.muonEnergy());                                  
      Jet_muonEnergyFraction.push_back(j.muonEnergyFraction());                                  
      Jet_electronEnergy.push_back(j.electronEnergy());                               
      Jet_photonEnergy.push_back(j.photonEnergy());                    
      if(j.isCaloJet()) Jet_emEnergyFraction.push_back(j.emEnergyFraction());
      else              Jet_emEnergyFraction.push_back(99);
      //Jet constituent multiplicity
      Jet_numberOfConstituents.push_back(j.chargedMultiplicity() + j.neutralMultiplicity());                                  
      Jet_chargedMultiplicity.push_back(j.chargedMultiplicity());
      Jet_isPFJet.push_back(j.isPFJet());
    }
  } 
}
void JetSelector::SetBranches(){
  if(debug_) std::cout<<"setting branches: calling AddBranch of baseTree"<<std::endl;
  //Kinematics
  AddBranch(&Jet_pt,                  "Jet_pt");
  AddBranch(&Jet_px,                  "Jet_px");
  AddBranch(&Jet_py,                  "Jet_py");
  AddBranch(&Jet_eta,                 "Jet_eta");
  AddBranch(&Jet_phi,                 "Jet_phi");
  AddBranch(&Jet_energy,              "Jet_energy");
  AddBranch(&Jet_mass,                "Jet_mass");
  AddBranch(&UncorrJet_pt,            "UncorrJet_pt");
  //ID
  AddBranch(&Jet_bDiscriminator,      "Jet_bDiscriminator");
  AddBranch(&Jet_bDiscriminator1,     "Jet_bDiscriminator1");
  AddBranch(&Jet_bDiscriminator2,     "Jet_bDiscriminator2");
  AddBranch(&Jet_pileupId,            "Jet_pileupId");
  //Jet Uncertainties
  AddBranch(&Jet_JesUp,               "Jet_JesUp");
  AddBranch(&Jet_JesDown,             "Jet_JesDown");
  //JER scale factor and uncertainties
  AddBranch(&Jet_JerSF,               "Jet_JerSF");
  AddBranch(&Jet_JerSFup,             "Jet_JerSFup");
  AddBranch(&Jet_JerSFdown,           "Jet_JerSFdown");
  //Energy related variables
  if(!_super_TNT){
    AddBranch(&Jet_neutralHadEnergyFraction,    "Jet_neutralHadEnergyFraction");
    AddBranch(&Jet_neutralEmEmEnergyFraction,   "Jet_neutralEmEmEnergyFraction");
    AddBranch(&Jet_chargedHadronEnergyFraction, "Jet_chargedHadronEnergyFraction");
    AddBranch(&Jet_chargedEmEnergyFraction,     "Jet_chargedEmEnergyFraction");
    AddBranch(&Jet_muonEnergyFraction,          "Jet_muonEnergyFraction");
    AddBranch(&Jet_electronEnergy,      "Jet_electronEnergy");
    AddBranch(&Jet_photonEnergy,        "Jet_photonEnergy");
    AddBranch(&Jet_emEnergyFraction,    "Jet_emEnergyFraction");
    //Jet constituent multiplicity
    AddBranch(&Jet_numberOfConstituents,"Jet_numberOfConstituents");
    AddBranch(&Jet_chargedMultiplicity, "Jet_chargedMultiplicity");
    AddBranch(&Jet_isPFJet,             "Jet_isPFJet");
  }
  if(debug_) std::cout<<"set branches"<<std::endl;
}
void JetSelector::Clear(){
  //Kinematics
  Jet_pt.clear();
  Jet_px.clear();
  Jet_py.clear();
  Jet_eta.clear();
  Jet_phi.clear();
  Jet_energy.clear();
  Jet_mass.clear();
  UncorrJet_pt.clear();
  //ID
  Jet_bDiscriminator.clear();
  Jet_bDiscriminator1.clear();
  Jet_bDiscriminator2.clear();
  Jet_pileupId.clear();
  //Jet Uncertainties
  Jet_JesUp.clear();
  Jet_JesDown.clear();
  //JER scale factor and uncertainties
  Jet_JerSF.clear();
  Jet_JerSFup.clear();
  Jet_JerSFdown.clear();
  //Energy related variables
  Jet_neutralHadEnergyFraction.clear();
  Jet_neutralEmEmEnergyFraction.clear();
  Jet_chargedHadronEnergyFraction.clear();
  Jet_chargedEmEnergyFraction.clear();
  Jet_muonEnergyFraction.clear();
  Jet_electronEnergy.clear();
  Jet_photonEnergy.clear();
  Jet_emEnergyFraction.clear();
  //Jet constituent multiplicity
  Jet_numberOfConstituents.clear();
  Jet_chargedMultiplicity.clear();
  Jet_isPFJet.clear();
}
void JetSelector::GetJESUncertainties(pat::Jet jet, JetCorrectionUncertainty *jecUnc, float &JesUncertainties){
  jecUnc->setJetEta(jet.eta());
  jecUnc->setJetPt(jet.pt()); // here you must use the CORRECTED jet pt
  JesUncertainties = jecUnc->getUncertainty(true);
}
void JetSelector::GetJER(pat::Jet jet, float &JERScaleFactor, float &JERScaleFactorUP, float &JERScaleFactorDOWN){
  if(!jet.genJet()) return;
  double jetEta=fabs(jet.eta());
  double cFactorJER = 1.0; 
  double cFactorJERdown = 1.0;
  double cFactorJERup = 1.0;
  //The following factors are derived from 8TeV but blessed for 13TeV, before new factors will be available
  //https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution#JER_Scaling_factors_and_Unce_AN1
  if( jetEta<0.5 ){ 
    cFactorJER = 1.079; 
    cFactorJERdown = 1.053;
    cFactorJERup = 1.105; 
  }
  else if( jetEta<1.1 ){ 
    cFactorJER = 1.099; 
    cFactorJERdown = 1.071;
    cFactorJERup = 1.127; 
  }
  else if( jetEta<1.7 ){ 
    cFactorJER = 1.121; 
    cFactorJERdown = 1.092;
    cFactorJERup = 1.150; 
  }
  else if( jetEta<2.3 ){ 
    cFactorJER = 1.208; 
    cFactorJERdown = 1.162;
    cFactorJERup = 1.254; 
  }
  else if( jetEta<2.8 ){ 
    cFactorJER = 1.254; 
    cFactorJERdown = 1.192;
    cFactorJERup = 1.316; 
  }
  else if( jetEta<3.2 ){ 
    cFactorJER = 1.395; 
    cFactorJERdown = 1.332;
    cFactorJERup = 1.458; 
  }
  else if( jetEta<5.0 ){ 
    cFactorJER = 1.056; 
    cFactorJERdown = 0.865;
    cFactorJERup = 1.247; 
  }
  double recoJetPt = jet.pt();
  double genJetPt  = jet.genJet()->pt();
  double diffPt    = recoJetPt - genJetPt;
  if(genJetPt>10.){
    JERScaleFactor     = (std::max(0., genJetPt + cFactorJER*diffPt))/recoJetPt;
    JERScaleFactorUP   = (std::max(0., genJetPt + cFactorJERup*diffPt))/recoJetPt;
    JERScaleFactorDOWN = (std::max(0., genJetPt + cFactorJERdown*diffPt))/recoJetPt;
  } else {
    JERScaleFactor     = 1.;
    JERScaleFactorUP   = 1.;
    JERScaleFactorDOWN = 1.;
  } 
}
