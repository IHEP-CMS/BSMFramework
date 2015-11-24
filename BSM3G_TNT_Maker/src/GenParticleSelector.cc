#include "BSMFramework/BSM3G_TNT_Maker/interface/GenParticleSelector.h"
GenParticleSelector::GenParticleSelector(std::string name, TTree* tree, bool debug, const pset& iConfig):baseTree(name,tree,debug){
  _tthlepVar = iConfig.getParameter<bool>("tthlepVar");
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
  Handle< reco::GenParticleCollection > _Gen_collection;
  iEvent.getByLabel("prunedGenParticles", _Gen_collection);
  /////
  //   Get gen information
  /////
  for(reco::GenParticleCollection::const_iterator genparticles = _Gen_collection->begin(); genparticles !=  _Gen_collection->end(); ++genparticles){
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
    for(reco::GenParticleCollection::const_iterator mit = _Gen_collection->begin(); mit != _Gen_collection->end(); ++mit){
      if(genparticles->mother() == &(*mit)){
	idx = std::distance(_Gen_collection->begin(),mit);
	break;
      }
    }
    Gen_BmotherIndex.push_back(idx);
    for(size_t j = 0; j < genparticles->numberOfMothers(); ++j){
      const reco::Candidate* m = genparticles->mother(j);
      for(reco::GenParticleCollection::const_iterator mit = _Gen_collection->begin();  mit != _Gen_collection->end(); ++mit){
	if(m == &(*mit)){ 
	  int idx = std::distance(_Gen_collection->begin(), mit);
	  Gen_BmotherIndices.push_back(idx);
	  break;
	}
      }
    }
    for(size_t j = 0; j < genparticles->numberOfDaughters(); ++j){
      const reco::Candidate* d = genparticles->daughter(j);
      for(reco::GenParticleCollection::const_iterator mit = _Gen_collection->begin();  mit != _Gen_collection->end(); ++mit){
	if(d == &(*mit)){ 
	  int idx = std::distance(_Gen_collection->begin(), mit);
	  Gen_BdaughtIndices.push_back(idx);
	  break;
	}
      }
    }
  }
  if(_tthlepVar){
    //Higgs decays
    int Hdecay = -1;
    if(_Gen_collection.isValid() ){
      Hdecay=0;
      for(size_t k = 0; k < _Gen_collection->size(); k++){
        //get current candidate
        const reco::Candidate & mcParticle = (*_Gen_collection)[k];
        int status = mcParticle.status();
        int pdgId  = mcParticle.pdgId();
        int absId  = abs(pdgId);
        int numdgt = mcParticle.numberOfDaughters();
        //it must be a Higgs and status 62(pythia 8), with at least 2 daughters
        if(absId!=25 || status!=62) continue;
        if(!(numdgt>1)) continue;
        //Get the daughters
        int ind0 = 0; int ind1 = 1;
        if(numdgt>2){
          if(mcParticle.daughter(0)->pdgId()==pdgId){
            ind0 = 1;
            ind1 = 2;
          }
          if(mcParticle.daughter(1)->pdgId()==pdgId){
            ind0 = 0;
            ind1 = 2;
          }
        }
        int d0 = -99; int d1 = -99;
        d0 = mcParticle.daughter(ind0)->pdgId();
        d1 = mcParticle.daughter(ind1)->pdgId();
        d0 = abs(d0);
        d1 = abs(d1);

        if(d0==5  && d1==5 ) Hdecay = 1; //bb
        if(d0==24 && d1==24) Hdecay = 2; //WW
        if(d0==15 && d1==15) Hdecay = 3; //TauTau
        if(d0==21 && d1==21) Hdecay = 4; //glueglue
        if(d0==4  && d1==4 ) Hdecay = 5; //cc
        if(d0==23 && d1==23) Hdecay = 6; //ZZ

        if(d0==22 && d1==23) Hdecay = 7; //Zy
        if(d0==23 && d1==22) Hdecay = 7; //Zy
 
        if(d0==22 && d1==22) Hdecay = 8;  //yy
        if(d0==21 && d1==22) Hdecay = 9;  //gy
        if(d0==22 && d1==21) Hdecay = 9;  //gy
        if(d0==3  && d1==3)  Hdecay = 10; //ss

        if(d0==13 && d1==13) Hdecay = 11; //mumu

        if((Hdecay==0) && (d0==22 || d1==22)) Hdecay = 12; //?y
        if((Hdecay==0) && (d0==21 || d1==21)) Hdecay = 13; //?g
        if((Hdecay==0) && (d0>100 || d1>100)) Hdecay = 14; //?Hadron

        if(d0==1  && d1==1)  Hdecay = 15; //uu
        if(d0==2  && d1==2)  Hdecay = 16; //dd
        if(d0==6  && d1==6)  Hdecay = 17; //tt
        if(d0==11 && d1==11) Hdecay = 18; //ee
      }
    }
    HiggsDecay = Hdecay;
  }
  if(debug_) std::cout<<"got gen particle  info"<<std::endl;
}
void GenParticleSelector::SetBranches(){
  if(debug_) std::cout<<"setting branches: calling AddBranch of baseTree"<<std::endl;
  //Kinematics
  AddBranch(&Gen_pt               ,"Gen_pt");
  AddBranch(&Gen_eta              ,"Gen_eta");
  AddBranch(&Gen_phi              ,"Gen_phi");
  AddBranch(&Gen_p                ,"Gen_p");
  AddBranch(&Gen_energy           ,"Gen_energy");
  //Charge
  AddBranch(&Gen_charge           ,"Gen_charge");
  //Vertex
  AddBranch(&Gen_vx               ,"Gen_vx");
  AddBranch(&Gen_vy               ,"Gen_vy");
  AddBranch(&Gen_vz               ,"Gen_vz");
  //Origin
  AddBranch(&Gen_status           ,"Gen_status");
  AddBranch(&Gen_pdg_id           ,"Gen_pdg_id");
  AddBranch(&Gen_motherpdg_id     ,"Gen_motherpdg_id");
  AddBranch(&Gen_numDaught        ,"Gen_numDaught");
  AddBranch(&Gen_numMother        ,"Gen_numMother");
  AddBranch(&Gen_BmotherIndex     ,"Gen_BmotherIndex");
  AddBranch(&Gen_BmotherIndices   ,"Gen_BmotherIndices");
  AddBranch(&Gen_BdaughtIndices   ,"Gen_BdaughtIndices");
  //TTHLep
  if(_tthlepVar){
    AddBranch(&HiggsDecay           ,"HiggsDecay");
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
  //TTHLep
  if(_tthlepVar){
    HiggsDecay = -9999;
  }
}
