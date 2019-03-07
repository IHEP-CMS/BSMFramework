#include "BSMFramework/BSM3G_TNT_Maker/interface/GenParticleSelector.h"
GenParticleSelector::GenParticleSelector(std::string name, TTree* tree, bool debug, const pset& iConfig, edm::ConsumesCollector && ic):
  baseTree(name,tree,debug)
{
  prunedGenToken_ = ic.consumes<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("pruned"));
  _reduced = iConfig.getParameter<bool>("reduced");
  if(debug) std::cout<<"in GenParticleSelector constructor"<<std::endl;
  if(debug) std::cout<<"in pileup constructor: calling SetBrances()"<<std::endl;
  SetBranches();
}
GenParticleSelector::~GenParticleSelector(){
  delete tree_;
}
void GenParticleSelector::Fill(const edm::Event& iEvent){
  Clear(); 
  if(debug_) std::cout<<"getting gen particle info"<<std::endl;
  /////
  //   Recall collections
  /////
  Handle<edm::View<reco::GenParticle> > pruned;
  iEvent.getByToken(prunedGenToken_, pruned);
  /////
  //   Get gen information
  /////
  for(size_t i=0; i<pruned->size(); i++){
    const Candidate * genparticles = &(*pruned)[i];
    //Kinematics
    Gen_pt.push_back(genparticles->pt());
    Gen_eta.push_back(genparticles->eta()); 
    Gen_phi.push_back(genparticles->phi());
    Gen_p.push_back(genparticles->p());
    Gen_energy.push_back(genparticles->energy());
    //Charge
    Gen_charge.push_back(genparticles->charge());
    //Vertex
    Gen_vx.push_back(genparticles->vx());
    Gen_vy.push_back(genparticles->vy());
    Gen_vz.push_back(genparticles->vz());
    //Origin
    Gen_status.push_back(genparticles->status());
    Gen_pdg_id.push_back(genparticles->pdgId());
    Gen_motherpdg_id.push_back(genparticles->numberOfMothers() > 0 ? genparticles->mother(0)->pdgId() : -999999);
    Gen_numDaught.push_back(genparticles->numberOfDaughters());
    Gen_numMother.push_back(genparticles->numberOfMothers());
    int idx = -1;
    for(size_t k = 0; k < pruned->size(); k++){
      const Candidate * mit = &(*pruned)[k];
      if(genparticles->mother() == mit){
        idx = k;
	break;
      }
    }
    Gen_BmotherIndex.push_back(idx);
    for(size_t j = 0; j < genparticles->numberOfMothers(); ++j){
      const reco::Candidate* m = genparticles->mother(j);
      for(size_t k = 0; k < pruned->size(); k++){
        const Candidate * mit = &(*pruned)[k];
	if(m == mit){ 
	  int idx = k;
	  Gen_BmotherIndices.push_back(idx);
	  break;
	}
      }
    }
    for(size_t j = 0; j < genparticles->numberOfDaughters(); ++j){
      const reco::Candidate* d = genparticles->daughter(j);
      for(size_t k = 0; k < pruned->size(); k++){
        const Candidate * mit = &(*pruned)[k];
	if(d == mit){ 
	  int idx = k;
	  Gen_BdaughtIndices.push_back(idx);
	  break;
	}
      }
    }
  }
  if(debug_) std::cout<<"got gen particle  info"<<std::endl;
}
void GenParticleSelector::SetBranches(){
  if(debug_) std::cout<<"setting branches: calling AddBranch of baseTree"<<std::endl;
  //Kinematics
  AddBranch(&Gen_pt               ,"Gen_pt");
  AddBranch(&Gen_eta              ,"Gen_eta");
  AddBranch(&Gen_phi              ,"Gen_phi");
  if(!_reduced){
    AddBranch(&Gen_p                ,"Gen_p");
    AddBranch(&Gen_energy           ,"Gen_energy");
    //Charge
    AddBranch(&Gen_charge           ,"Gen_charge");
    //Vertex
    AddBranch(&Gen_vx               ,"Gen_vx");
    AddBranch(&Gen_vy               ,"Gen_vy");
    AddBranch(&Gen_vz               ,"Gen_vz");
  }
  //Origin
  AddBranch(&Gen_status           ,"Gen_status");
  AddBranch(&Gen_pdg_id           ,"Gen_pdg_id");
  AddBranch(&Gen_motherpdg_id     ,"Gen_motherpdg_id");
  if(!_reduced){
    AddBranch(&Gen_numDaught        ,"Gen_numDaught");
    AddBranch(&Gen_numMother        ,"Gen_numMother");
    AddBranch(&Gen_BmotherIndex     ,"Gen_BmotherIndex");
    AddBranch(&Gen_BmotherIndices   ,"Gen_BmotherIndices");
    AddBranch(&Gen_BdaughtIndices   ,"Gen_BdaughtIndices");
  }
  if(debug_) std::cout<<"set branches genparticle"<<std::endl;
}
void GenParticleSelector::Clear(){
  //Kinematics
  Gen_pt.clear();
  Gen_eta.clear();
  Gen_phi.clear();
  Gen_p.clear();
  Gen_energy.clear();
  //Charge
  Gen_charge.clear();
  //Vertex
  Gen_vx.clear();
  Gen_vy.clear();
  Gen_vz.clear();
  //Origin
  Gen_status.clear();
  Gen_pdg_id.clear();
  Gen_motherpdg_id.clear();
  Gen_numDaught.clear();
  Gen_numMother.clear();
  Gen_BmotherIndex.clear();
  Gen_BmotherIndices.clear();
  Gen_BdaughtIndices.clear();
}
