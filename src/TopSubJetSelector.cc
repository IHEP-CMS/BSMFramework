#include "NtupleMaker/BSM3G_TNT_Maker/interface/TopSubJetSelector.h"
TopSubJetSelector::TopSubJetSelector(std::string name, TTree* tree, bool debug, const pset& iConfig):baseTree(name,tree,debug){
 topsubjetToken_ = iConfig.getParameter<edm::InputTag>("topsubjets");
 _Jet_pt_min     = iConfig.getParameter<double>("Jet_pt_min");
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
 } 
}
void TopSubJetSelector::SetBranches(){
 if(debug_) std::cout<<"setting branches: calling AddBranch of baseTree"<<std::endl;
 //Kinematics
 AddBranch(&TopSubjet_pt,     "TopSubjet_pt");
 AddBranch(&TopSubjet_eta,    "TopSubjet_eta");
 AddBranch(&TopSubjet_phi,    "TopSubjet_phi");
 AddBranch(&TopSubjet_energy, "TopSubjet_energy");
 AddBranch(&TopSubjet_mass,   "TopSubjet_mass");
 //ID
 AddBranch(&TopSubjet_Btag0,  "TopSubjet_Btag0");
 AddBranch(&TopSubjet_Btag1,  "TopSubjet_Btag1");
 AddBranch(&TopSubjet_Btag2,  "TopSubjet_Btag2");
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
}
