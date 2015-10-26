#include "BSMFramework/BSM3G_TNT_Maker/interface/JetSelector.h"
JetSelector::JetSelector(std::string name, TTree* tree, bool debug, const pset& iConfig):baseTree(name,tree,debug){
  jetToken_       = iConfig.getParameter<edm::InputTag>("jets");
  puppi_jetToken_ = iConfig.getParameter<edm::InputTag>("jetsPUPPI");
  _vertexInputTag = iConfig.getParameter<edm::InputTag>("vertices");
  jecfile_        = iConfig.getParameter<edm::FileInPath>("jecfile");   
  _Jet_pt_min     = iConfig.getParameter<double>("Jet_pt_min");
  _super_TNT      = iConfig.getParameter<bool>("super_TNT");
  SetBranches();
}
JetSelector::~JetSelector(){
  delete tree_;
}
void JetSelector::Fill(const edm::Event& iEvent){
  Clear();
  /////
  //   Recall collections
  /////  
  edm::Handle<pat::JetCollection> jets;                                       
  iEvent.getByLabel(jetToken_, jets);                                         
  edm::Handle<pat::JetCollection> puppijets;                                       
  iEvent.getByLabel(puppi_jetToken_, puppijets); 
  /////
  //   Get jet information
  /////  
  ////slimmedJets
  for(const pat::Jet &j : *jets){ 
    //Acceptance
    if(j.pt()<_Jet_pt_min) continue;
    //Kinematics
    Jet_pt.push_back(j.pt());  
    Jet_eta.push_back(j.eta());       
    Jet_phi.push_back(j.phi());       
    Jet_energy.push_back(j.energy());
    Jet_mass.push_back(j.mass()); 
    Jet_px.push_back(j.px());   
    Jet_py.push_back(j.py());          
    Jet_pz.push_back(j.pz());          
    Jet_Uncorr_pt.push_back(j.correctedJet("Uncorrected").pt());                
    //ID
    Jet_pfCombinedInclusiveSecondaryVertexV2BJetTags.push_back(j.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
    Jet_pfCombinedMVABJetTags.push_back(j.bDiscriminator("pfCombinedMVABJetTags"));
    Jet_pfJetProbabilityBJetTags.push_back(j.bDiscriminator("pfJetProbabilityBJetTags"));
    Jet_pileupId.push_back(j.userFloat("pileupJetId:fullDiscriminant"));
    Jet_isPFJet.push_back(j.isPFJet());
    Jet_isCaloJet.push_back(j.isCaloJet());
    //Energy
    Jet_neutralHadEnergyFraction.push_back(j.neutralHadronEnergyFraction());                               
    Jet_neutralEmEnergyFraction.push_back(j.neutralEmEnergyFraction());                                   
    Jet_chargedHadronEnergyFraction.push_back(j.chargedHadronEnergyFraction());                               
    Jet_chargedEmEnergyFraction.push_back(j.chargedEmEnergyFraction());                              
    Jet_muonEnergyFraction.push_back(j.muonEnergyFraction());                                  
    Jet_electronEnergy.push_back(j.electronEnergy());                               
    Jet_photonEnergy.push_back(j.photonEnergy());                                 
    if(j.isCaloJet()) Jet_emEnergyFraction.push_back(j.emEnergyFraction());
    else              Jet_emEnergyFraction.push_back(-999);
    //Other prop
    Jet_numberOfConstituents.push_back(j.chargedMultiplicity() + j.neutralMultiplicity());                                  
    Jet_chargedMultiplicity.push_back(j.chargedMultiplicity());
    Jet_vtxMass.push_back(j.userFloat("vtxMass"));
    Jet_vtxNtracks.push_back(j.userFloat("vtxNtracks"));
    Jet_vtx3DVal.push_back(j.userFloat("vtx3DVal"));
    Jet_vtx3DSig.push_back(j.userFloat("vtx3DSig"));
    //Corrections/Systematics
    JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty(jecfile_.fullPath());
    //Jet Uncertainties
    float JesUncertainties = 0;
    GetJESUncertainties(j, jecUnc, JesUncertainties);
    Jet_JesUp.push_back((1+JesUncertainties));         
    Jet_JesDown.push_back((1-JesUncertainties));
    //JER scale factor and uncertainties
    float JERScaleFactor     = 1; 
    float JERScaleFactorUP   = 1;
    float JERScaleFactorDOWN = 1;
    GetJER(j, JERScaleFactor, JERScaleFactorUP, JERScaleFactorDOWN);
    Jet_JerSF.push_back(JERScaleFactor);
    Jet_JerSFup.push_back(JERScaleFactorUP);
    Jet_JerSFdown.push_back(JERScaleFactorDOWN);
    delete jecUnc;
    //MC
    Jet_partonFlavour.push_back(j.partonFlavour());
    Jet_hadronFlavour.push_back(j.hadronFlavour());
  } 
  ////slimmedJetsPuppi
  for(const pat::Jet &j : *puppijets){ 
    //Acceptance
    if(j.pt() < _Jet_pt_min) continue;
    //Kinematics
    Jet_puppi_pt.push_back(j.pt());  
    Jet_puppi_eta.push_back(j.eta());       
    Jet_puppi_phi.push_back(j.phi());       
    Jet_puppi_energy.push_back(j.energy());
    Jet_puppi_mass.push_back(j.mass()); 
    Jet_puppi_px.push_back(j.px());   
    Jet_puppi_py.push_back(j.py());          
    Jet_puppi_pz.push_back(j.pz());          
    Jet_puppi_Uncorr_pt.push_back(j.correctedJet("Uncorrected").pt());                
    //ID
    Jet_puppi_pfCombinedInclusiveSecondaryVertexV2BJetTags.push_back(j.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
    Jet_puppi_pfCombinedMVABJetTags.push_back(j.bDiscriminator("pfCombinedMVABJetTags"));
    Jet_puppi_pfJetProbabilityBJetTags.push_back(j.bDiscriminator("pfJetProbabilityBJetTags"));
    Jet_puppi_pileupId.push_back(j.userFloat("pileupJetId:fullDiscriminant"));
    Jet_puppi_isPFJet.push_back(j.isPFJet());
    Jet_puppi_isCaloJet.push_back(j.isCaloJet());
    //Energy
    Jet_puppi_neutralHadEnergyFraction.push_back(j.neutralHadronEnergyFraction());                               
    Jet_puppi_neutralEmEnergyFraction.push_back(j.neutralEmEnergyFraction());                                   
    Jet_puppi_chargedHadronEnergyFraction.push_back(j.chargedHadronEnergyFraction());                               
    Jet_puppi_chargedEmEnergyFraction.push_back(j.chargedEmEnergyFraction());                              
    Jet_puppi_muonEnergyFraction.push_back(j.muonEnergyFraction());                                  
    Jet_puppi_electronEnergy.push_back(j.electronEnergy());                               
    Jet_puppi_photonEnergy.push_back(j.photonEnergy());                                 
    if(j.isCaloJet()) Jet_puppi_emEnergyFraction.push_back(j.emEnergyFraction());
    else              Jet_puppi_emEnergyFraction.push_back(-999);
    //Other prop
    Jet_puppi_numberOfConstituents.push_back(j.chargedMultiplicity() + j.neutralMultiplicity());                                  
    Jet_puppi_chargedMultiplicity.push_back(j.chargedMultiplicity());
    Jet_puppi_vtxMass.push_back(j.userFloat("vtxMass"));
    Jet_puppi_vtxNtracks.push_back(j.userFloat("vtxNtracks"));
    Jet_puppi_vtx3DVal.push_back(j.userFloat("vtx3DVal"));
    Jet_puppi_vtx3DSig.push_back(j.userFloat("vtx3DSig"));
    //Corrections/Systematics
    JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty(jecfile_.fullPath());
    //Jet Uncertainties
    float JesUncertainties = 0;
    GetJESUncertainties(j, jecUnc, JesUncertainties);
    Jet_puppi_JesUp.push_back((1+JesUncertainties));         
    Jet_puppi_JesDown.push_back((1-JesUncertainties));
    //JER scale factor and uncertainties
    float JERScaleFactor     = 1; 
    float JERScaleFactorUP   = 1;
    float JERScaleFactorDOWN = 1;
    GetJER(j, JERScaleFactor, JERScaleFactorUP, JERScaleFactorDOWN);
    Jet_puppi_JerSF.push_back(JERScaleFactor);
    Jet_puppi_JerSFup.push_back(JERScaleFactorUP);
    Jet_puppi_JerSFdown.push_back(JERScaleFactorDOWN);
    delete jecUnc;
    //MC
    Jet_puppi_partonFlavour.push_back(j.partonFlavour());
    Jet_puppi_hadronFlavour.push_back(j.hadronFlavour());
  } 
}
void JetSelector::SetBranches(){
  if(debug_) std::cout<<"setting branches: calling AddBranch of baseTree"<<std::endl;
  ////slimmedJets
  //Kinematics
  AddBranch(&Jet_pt        ,"Jet_pt");
  AddBranch(&Jet_eta       ,"Jet_eta");
  AddBranch(&Jet_phi       ,"Jet_phi");
  AddBranch(&Jet_energy    ,"Jet_energy");
  AddBranch(&Jet_mass      ,"Jet_mass");
  AddBranch(&Jet_px        ,"Jet_px");
  AddBranch(&Jet_py        ,"Jet_py");
  AddBranch(&Jet_pz        ,"Jet_pz");
  AddBranch(&Jet_Uncorr_pt ,"Jet_Uncorr_pt");
  //ID
  AddBranch(&Jet_pfCombinedInclusiveSecondaryVertexV2BJetTags ,"Jet_pfCombinedInclusiveSecondaryVertexV2BJetTags");
  AddBranch(&Jet_pfCombinedMVABJetTags                        ,"Jet_pfCombinedMVABJetTags");
  AddBranch(&Jet_pfJetProbabilityBJetTags                     ,"Jet_pfJetProbabilityBJetTags");
  AddBranch(&Jet_pileupId                                     ,"Jet_pileupId");
  AddBranch(&Jet_isPFJet                                      ,"Jet_isPFJet");
  AddBranch(&Jet_isCaloJet                                    ,"Jet_isCaloJet");
  //Energy
  AddBranch(&Jet_neutralHadEnergyFraction    ,"Jet_neutralHadEnergyFraction");
  AddBranch(&Jet_neutralEmEnergyFraction     ,"Jet_neutralEmEnergyFraction");
  AddBranch(&Jet_chargedHadronEnergyFraction ,"Jet_chargedHadronEnergyFraction");
  AddBranch(&Jet_chargedEmEnergyFraction     ,"Jet_chargedEmEnergyFraction");
  AddBranch(&Jet_muonEnergyFraction          ,"Jet_muonEnergyFraction");
  AddBranch(&Jet_electronEnergy              ,"Jet_electronEnergy");
  AddBranch(&Jet_photonEnergy                ,"Jet_photonEnergy");
  AddBranch(&Jet_emEnergyFraction            ,"Jet_emEnergyFraction");
  //Other prop
  AddBranch(&Jet_numberOfConstituents ,"Jet_numberOfConstituents");
  AddBranch(&Jet_chargedMultiplicity  ,"Jet_chargedMultiplicity");
  AddBranch(&Jet_vtxMass              ,"Jet_vtxMass");
  AddBranch(&Jet_vtxNtracks           ,"Jet_vtxNtracks");
  AddBranch(&Jet_vtx3DVal             ,"Jet_vtx3DVal");
  AddBranch(&Jet_vtx3DSig             ,"Jet_vtx3DSig");
  //Corrections/Systematics
  AddBranch(&Jet_JesUp                ,"Jet_JesUp");
  AddBranch(&Jet_JesDown              ,"Jet_JesDown");
  AddBranch(&Jet_JerSF                ,"Jet_JerSF");
  AddBranch(&Jet_JerSFup              ,"Jet_JerSFup");
  AddBranch(&Jet_JerSFdown            ,"Jet_JerSFdown");
  //MC
  AddBranch(&Jet_partonFlavour        ,"Jet_partonFlavour");
  AddBranch(&Jet_hadronFlavour        ,"Jet_hadronFlavour");
  ////slimmedJetsPuppi
  //Kinematics
  AddBranch(&Jet_puppi_pt        ,"Jet_puppi_pt");
  AddBranch(&Jet_puppi_eta       ,"Jet_puppi_eta");
  AddBranch(&Jet_puppi_phi       ,"Jet_puppi_phi");
  AddBranch(&Jet_puppi_energy    ,"Jet_puppi_energy");
  AddBranch(&Jet_puppi_mass      ,"Jet_puppi_mass");
  AddBranch(&Jet_puppi_px        ,"Jet_puppi_px");
  AddBranch(&Jet_puppi_py        ,"Jet_puppi_py");
  AddBranch(&Jet_puppi_pz        ,"Jet_puppi_pz");
  AddBranch(&Jet_puppi_Uncorr_pt ,"Jet_puppi_Uncorr_pt");
  //ID
  AddBranch(&Jet_puppi_pfCombinedInclusiveSecondaryVertexV2BJetTags ,"Jet_puppi_pfCombinedInclusiveSecondaryVertexV2BJetTags");
  AddBranch(&Jet_puppi_pfCombinedMVABJetTags                        ,"Jet_puppi_pfCombinedMVABJetTags");
  AddBranch(&Jet_puppi_pfJetProbabilityBJetTags                     ,"Jet_puppi_pfJetProbabilityBJetTags");
  AddBranch(&Jet_puppi_pileupId                                     ,"Jet_puppi_pileupId");
  AddBranch(&Jet_puppi_isPFJet                                      ,"Jet_puppi_isPFJet");
  AddBranch(&Jet_puppi_isCaloJet                                    ,"Jet_puppi_isCaloJet");
  //Energy
  AddBranch(&Jet_puppi_neutralHadEnergyFraction    ,"Jet_puppi_neutralHadEnergyFraction");
  AddBranch(&Jet_puppi_neutralEmEnergyFraction     ,"Jet_puppi_neutralEmEnergyFraction");
  AddBranch(&Jet_puppi_chargedHadronEnergyFraction ,"Jet_puppi_chargedHadronEnergyFraction");
  AddBranch(&Jet_puppi_chargedEmEnergyFraction     ,"Jet_puppi_chargedEmEnergyFraction");
  AddBranch(&Jet_puppi_muonEnergyFraction          ,"Jet_puppi_muonEnergyFraction");
  AddBranch(&Jet_puppi_electronEnergy              ,"Jet_puppi_electronEnergy");
  AddBranch(&Jet_puppi_photonEnergy                ,"Jet_puppi_photonEnergy");
  AddBranch(&Jet_puppi_emEnergyFraction            ,"Jet_puppi_emEnergyFraction");
  //Other prop
  AddBranch(&Jet_puppi_numberOfConstituents ,"Jet_puppi_numberOfConstituents");
  AddBranch(&Jet_puppi_chargedMultiplicity  ,"Jet_puppi_chargedMultiplicity");
  AddBranch(&Jet_puppi_vtxMass              ,"Jet_puppi_vtxMass");
  AddBranch(&Jet_puppi_vtxNtracks           ,"Jet_puppi_vtxNtracks");
  AddBranch(&Jet_puppi_vtx3DVal             ,"Jet_puppi_vtx3DVal");
  AddBranch(&Jet_puppi_vtx3DSig             ,"Jet_puppi_vtx3DSig");
  //Corrections/Systematics
  AddBranch(&Jet_puppi_JesUp                ,"Jet_puppi_JesUp");
  AddBranch(&Jet_puppi_JesDown              ,"Jet_puppi_JesDown");
  AddBranch(&Jet_puppi_JerSF                ,"Jet_puppi_JerSF");
  AddBranch(&Jet_puppi_JerSFup              ,"Jet_puppi_JerSFup");
  AddBranch(&Jet_puppi_JerSFdown            ,"Jet_puppi_JerSFdown");
  //MC
  AddBranch(&Jet_puppi_partonFlavour        ,"Jet_puppi_partonFlavour");
  AddBranch(&Jet_puppi_hadronFlavour        ,"Jet_puppi_hadronFlavour");
  if(debug_) std::cout<<"set branches"<<std::endl;
}
void JetSelector::Clear(){
  ////slimmedJets
  //Kinematics
  Jet_pt.clear();
  Jet_eta.clear();
  Jet_phi.clear();
  Jet_energy.clear();
  Jet_mass.clear();
  Jet_px.clear();
  Jet_py.clear();
  Jet_pz.clear();
  Jet_Uncorr_pt.clear();
  //ID
  Jet_pfCombinedInclusiveSecondaryVertexV2BJetTags.clear();
  Jet_pfCombinedMVABJetTags.clear();
  Jet_pfJetProbabilityBJetTags.clear();
  Jet_pileupId.clear();
  Jet_isPFJet.clear();
  Jet_isCaloJet.clear();
  //Energy
  Jet_neutralHadEnergyFraction.clear();
  Jet_neutralEmEnergyFraction.clear();
  Jet_chargedHadronEnergyFraction.clear();
  Jet_chargedEmEnergyFraction.clear();
  Jet_muonEnergyFraction.clear();
  Jet_electronEnergy.clear();
  Jet_photonEnergy.clear();
  Jet_emEnergyFraction.clear();
  //Other prop
  Jet_numberOfConstituents.clear();
  Jet_chargedMultiplicity.clear();
  Jet_vtxMass.clear();
  Jet_vtxNtracks.clear();
  Jet_vtx3DVal.clear();
  Jet_vtx3DSig.clear();
  //Corrections/Systematics
  Jet_JesUp.clear();
  Jet_JesDown.clear();
  Jet_JerSF.clear();
  Jet_JerSFup.clear();
  Jet_JerSFdown.clear(); 
  //MC
  Jet_partonFlavour.clear();
  Jet_hadronFlavour.clear();
  ////slimmedJetsPuppi
  //Kinematics
  Jet_puppi_pt.clear();
  Jet_puppi_eta.clear();
  Jet_puppi_phi.clear();
  Jet_puppi_energy.clear();
  Jet_puppi_mass.clear();
  Jet_puppi_px.clear();
  Jet_puppi_py.clear();
  Jet_puppi_pz.clear();
  Jet_puppi_Uncorr_pt.clear();
  //ID
  Jet_puppi_pfCombinedInclusiveSecondaryVertexV2BJetTags.clear();
  Jet_puppi_pfCombinedMVABJetTags.clear();
  Jet_puppi_pfJetProbabilityBJetTags.clear();
  Jet_puppi_pileupId.clear();
  Jet_puppi_isPFJet.clear();
  Jet_puppi_isCaloJet.clear();
  //Energy
  Jet_puppi_neutralHadEnergyFraction.clear();
  Jet_puppi_neutralEmEnergyFraction.clear();
  Jet_puppi_chargedHadronEnergyFraction.clear();
  Jet_puppi_chargedEmEnergyFraction.clear();
  Jet_puppi_muonEnergyFraction.clear();
  Jet_puppi_electronEnergy.clear();
  Jet_puppi_photonEnergy.clear();
  Jet_puppi_emEnergyFraction.clear();
  //Other prop
  Jet_puppi_numberOfConstituents.clear();
  Jet_puppi_chargedMultiplicity.clear();
  Jet_puppi_vtxMass.clear();
  Jet_puppi_vtxNtracks.clear();
  Jet_puppi_vtx3DVal.clear();
  Jet_puppi_vtx3DSig.clear();
  //Corrections/Systematics
  Jet_puppi_JesUp.clear();
  Jet_puppi_JesDown.clear();
  Jet_puppi_JerSF.clear();
  Jet_puppi_JerSFup.clear();
  Jet_puppi_JerSFdown.clear(); 
  //MC
  Jet_puppi_partonFlavour.clear();
  Jet_puppi_hadronFlavour.clear();
}
void JetSelector::GetJESUncertainties(pat::Jet jet, JetCorrectionUncertainty *jecUnc, float &JesUncertainties){
  jecUnc->setJetPt(jet.pt()); // here you must use the CORRECTED jet pt
  jecUnc->setJetEta(jet.eta());
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
