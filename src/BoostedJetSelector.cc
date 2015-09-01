#include "NtupleMaker/BSM3G_TNT_Maker/interface/BoostedJetSelector.h"
BoostedJetSelector::BoostedJetSelector(std::string name, TTree* tree, bool debug, const pset& iConfig):baseTree(name,tree,debug){
 fatjetToken_ = iConfig.getParameter<edm::InputTag>("fatjets");
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
  //ID
  BoostedJet_Btag0.push_back(j.bDiscriminator("combinedSecondaryVertexBJetTags"));              
  BoostedJet_Btag1.push_back(j.bDiscriminator("pfCombinedSecondaryVertexV2BJetTags"));              
  BoostedJet_Btag2.push_back(j.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));              
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

void BoostedJetSelector::SetBranches(){
 if(debug_)    std::cout<<"setting branches: calling AddBranch of baseTree"<<std::endl;
 //Kinematic
 AddBranch(&BoostedJet_pt,                          "BoostedJet_pt");
 AddBranch(&BoostedJet_eta,                         "BoostedJet_eta");
 AddBranch(&BoostedJet_phi,                         "BoostedJet_phi");
 AddBranch(&BoostedJet_energy,                      "BoostedJet_energy");
 AddBranch(&BoostedJet_mass,                        "BoostedJet_mass");
 //ID
 AddBranch(&BoostedJet_Btag0,                       "BoostedJet_Btag0");
 AddBranch(&BoostedJet_Btag2,                       "BoostedJet_Btag1");
 AddBranch(&BoostedJet_Btag2,                       "BoostedJet_Btag2");
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
 //ID
 BoostedJet_Btag0.clear();
 BoostedJet_Btag1.clear();
 BoostedJet_Btag2.clear();
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
 //Variables for top-tagging
 TopTagging_topMass.clear();
 TopTagging_minMass.clear();
 TopTagging_wMass.clear();
 TopTagging_nSubJets.clear();
}
