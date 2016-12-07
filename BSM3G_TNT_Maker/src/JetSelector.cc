#include "BSMFramework/BSM3G_TNT_Maker/interface/JetSelector.h"
JetSelector::JetSelector(std::string name, TTree* tree, bool debug, const pset& iConfig, edm::ConsumesCollector && ic):
  baseTree(name,tree,debug)
{
  vtx_h_        = ic.consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"));
  jets_         = ic.consumes<pat::JetCollection >(iConfig.getParameter<edm::InputTag>("jets"));
  puppijets_    = ic.consumes<pat::JetCollection >(iConfig.getParameter<edm::InputTag>("jetsPUPPI"));
  qgToken_      = ic.consumes<edm::ValueMap<float>>(edm::InputTag("QGTagger", "qgLikelihood"));
  axis2Token_   = ic.consumes<edm::ValueMap<float>>(edm::InputTag("QGTagger", "axis2"));
  ptDToken_     = ic.consumes<edm::ValueMap<float>>(edm::InputTag("QGTagger", "ptD"));
  multToken_    = ic.consumes<edm::ValueMap<int>>(edm::InputTag("QGTagger", "mult"));
  rhopogHandle_ = ic.consumes<double>(edm::InputTag("fixedGridRhoFastjetAll"));
  rhoJERHandle_ = ic.consumes<double>(edm::InputTag("fixedGridRhoAll"));
  jecPayloadNamesAK4PFchsMC1_   = iConfig.getParameter<edm::FileInPath>("jecPayloadNamesAK4PFchsMC1");
  jecPayloadNamesAK4PFchsMC2_   = iConfig.getParameter<edm::FileInPath>("jecPayloadNamesAK4PFchsMC2");
  jecPayloadNamesAK4PFchsMC3_   = iConfig.getParameter<edm::FileInPath>("jecPayloadNamesAK4PFchsMC3");
  jecPayloadNamesAK4PFchsMCUnc_ = iConfig.getParameter<edm::FileInPath>("jecPayloadNamesAK4PFchsMCUnc");
  jecPayloadNamesAK4PFchsDATA1_   = iConfig.getParameter<edm::FileInPath>("jecPayloadNamesAK4PFchsDATA1");
  jecPayloadNamesAK4PFchsDATA2_   = iConfig.getParameter<edm::FileInPath>("jecPayloadNamesAK4PFchsDATA2");
  jecPayloadNamesAK4PFchsDATA3_   = iConfig.getParameter<edm::FileInPath>("jecPayloadNamesAK4PFchsDATA3");
  jecPayloadNamesAK4PFchsDATA4_   = iConfig.getParameter<edm::FileInPath>("jecPayloadNamesAK4PFchsDATA4");
  jecPayloadNamesAK4PFchsDATAUnc_ = iConfig.getParameter<edm::FileInPath>("jecPayloadNamesAK4PFchsDATAUnc");
  jecPayloadNamesAK4PFPuppiMC1_   = iConfig.getParameter<edm::FileInPath>("jecPayloadNamesAK4PFPuppiMC1");
  jecPayloadNamesAK4PFPuppiMC2_   = iConfig.getParameter<edm::FileInPath>("jecPayloadNamesAK4PFPuppiMC2");
  jecPayloadNamesAK4PFPuppiMC3_   = iConfig.getParameter<edm::FileInPath>("jecPayloadNamesAK4PFPuppiMC3");
  jecPayloadNamesAK4PFPuppiMCUnc_ = iConfig.getParameter<edm::FileInPath>("jecPayloadNamesAK4PFPuppiMCUnc");
  jecPayloadNamesAK4PFPuppiDATA1_   = iConfig.getParameter<edm::FileInPath>("jecPayloadNamesAK4PFPuppiDATA1");
  jecPayloadNamesAK4PFPuppiDATA2_   = iConfig.getParameter<edm::FileInPath>("jecPayloadNamesAK4PFPuppiDATA2");
  jecPayloadNamesAK4PFPuppiDATA3_   = iConfig.getParameter<edm::FileInPath>("jecPayloadNamesAK4PFPuppiDATA3");
  jecPayloadNamesAK4PFPuppiDATA4_   = iConfig.getParameter<edm::FileInPath>("jecPayloadNamesAK4PFPuppiDATA4");
  jecPayloadNamesAK4PFPuppiDATAUnc_ = iConfig.getParameter<edm::FileInPath>("jecPayloadNamesAK4PFPuppiDATAUnc");
  jerAK4PFchs_     = iConfig.getParameter<edm::FileInPath>("jerAK4PFchs").fullPath();
  jerAK4PFchsSF_   = iConfig.getParameter<edm::FileInPath>("jerAK4PFchsSF").fullPath();
  jerAK4PFPuppi_   = iConfig.getParameter<edm::FileInPath>("jerAK4PFPuppi").fullPath();
  jerAK4PFPuppiSF_ = iConfig.getParameter<edm::FileInPath>("jerAK4PFPuppiSF").fullPath();
  _Jet_pt_min     = iConfig.getParameter<double>("Jet_pt_min");
  _super_TNT      = iConfig.getParameter<bool>("super_TNT");
  _is_data = iConfig.getParameter<bool>("is_data");
  _PuppiVar = iConfig.getParameter<bool>("PuppiVar");
  JECInitialization();
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
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vtx_h_, vertices);
  edm::Handle<pat::JetCollection> jets;                                       
  iEvent.getByToken(jets_, jets);                                         
  edm::Handle<edm::ValueMap<float>> qgHandle;
  iEvent.getByToken(qgToken_, qgHandle);
  edm::Handle<edm::ValueMap<float>> axis2Handle;
  iEvent.getByToken(axis2Token_, axis2Handle);
  edm::Handle<edm::ValueMap<float>> ptDHandle;
  iEvent.getByToken(ptDToken_, ptDHandle);
  edm::Handle<edm::ValueMap<int>> multHandle;
  iEvent.getByToken(multToken_, multHandle);
  edm::Handle<pat::JetCollection> puppijets;                                       
  iEvent.getByToken(puppijets_, puppijets); 
  edm::Handle<double> rhoHandle;
  iEvent.getByToken(rhopogHandle_,rhoHandle);
  double rho = *rhoHandle;
  edm::Handle<double> rhoJERHandle;
  iEvent.getByToken(rhoJERHandle_,rhoJERHandle);
  double rhoJER = *rhoJERHandle;
  /////
  //   Get jet information
  /////  
  //bool ajet = false;
  ////slimmedJets
  int ij = 0;
  for(const pat::Jet &j : *jets){ 
    //Acceptance
    if(j.pt()<_Jet_pt_min){ij++; continue;}
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
    if(!_is_data){
      if(j.genJet()){ 
        Jet_genpt.push_back(j.genJet()->pt());
        Jet_geneta.push_back(j.genJet()->eta());
        Jet_genphi.push_back(j.genJet()->phi());
        Jet_genenergy.push_back(j.genJet()->energy());
        Jet_genmass.push_back(j.genJet()->mass());
        Jet_genpx.push_back(j.genJet()->px());
        Jet_genpy.push_back(j.genJet()->py());
        Jet_genpz.push_back(j.genJet()->pz());
      }else{
        Jet_genpt.push_back(-999);
        Jet_geneta.push_back(-999);
        Jet_genphi.push_back(-999);
        Jet_genenergy.push_back(-999);
        Jet_genmass.push_back(-999);
        Jet_genpx.push_back(-999);
        Jet_genpy.push_back(-999);
        Jet_genpz.push_back(-999);
      }    
    }
    //ID
    Jet_newpfCombinedInclusiveSecondaryVertexV2BJetTags.push_back(j.bDiscriminator("newpfCombinedInclusiveSecondaryVertexV2BJetTags"));
    Jet_newpfCombinedMVAV2BJetTags.push_back(j.bDiscriminator("newpfCombinedMVAV2BJetTags"));
    Jet_newpfJetProbabilityBJetTags.push_back(j.bDiscriminator("newpfJetProbabilityBJetTags"));
    Jet_newpfCombinedCvsLJetTags.push_back(j.bDiscriminator("newpfCombinedCvsLJetTags"));
    Jet_newpfCombinedCvsBJetTags.push_back(j.bDiscriminator("newpfCombinedCvsBJetTags"));
    Jet_pfCombinedInclusiveSecondaryVertexV2BJetTags.push_back(j.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
    Jet_pfCombinedMVAV2BJetTags.push_back(j.bDiscriminator("pfCombinedMVAV2BJetTags"));
    Jet_pfJetProbabilityBJetTags.push_back(j.bDiscriminator("pfJetProbabilityBJetTags"));
    Jet_pfCombinedCvsLJetTags.push_back(j.bDiscriminator("pfCombinedCvsLJetTags"));
    Jet_pfCombinedCvsBJetTags.push_back(j.bDiscriminator("pfCombinedCvsBJetTags"));
    Jet_pileupId.push_back(j.userFloat("pileupJetId:fullDiscriminant"));
    Jet_isPFJet.push_back(j.isPFJet());
    Jet_isCaloJet.push_back(j.isCaloJet());
    edm::Ref<pat::JetCollection> jetRef(jets, ij);
    Jet_qg.push_back((*qgHandle)[jetRef]);
    Jet_axis2.push_back((*axis2Handle)[jetRef]);
    Jet_ptD.push_back((*ptDHandle)[jetRef]);
    Jet_mult.push_back((*multHandle)[jetRef]);
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
    Jet_JesSF.push_back(corrAK4PFchs);
    Jet_JesSFup.push_back(corrUpAK4PFchs);
    Jet_JesSFdown.push_back(corrDownAK4PFchs);
    //JER scale factor and uncertainties
    float JERScaleFactor     = 1; 
    float JERScaleFactorUP   = 1;
    float JERScaleFactorDOWN = 1;
    if(!_is_data) GetJER(j, corrAK4PFchs, rhoJER, true, JERScaleFactor, JERScaleFactorUP, JERScaleFactorDOWN);
    Jet_JerSF.push_back(JERScaleFactor);
    Jet_JerSFup.push_back(JERScaleFactorUP);
    Jet_JerSFdown.push_back(JERScaleFactorDOWN);
    //MC
    if(!_is_data) {
      Jet_partonFlavour.push_back(j.partonFlavour());
      Jet_hadronFlavour.push_back(j.hadronFlavour());
    }
    /////
    //   TTH variables
    /////
    //cout<<setiosflags(ios::fixed)<<setprecision(5);
    if(ij==0){
      //cout<<"jet val are "<<setw(20)<<iEvent.id().event()<<setw(20)<<j.pt()<<setw(20)<<j.eta()<<setw(20)<<j.phi()<<setw(20)<<j.energy()<<setw(20)<<j.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags")<<setw(20);
    //  ajet = true;
    }
    ij++;
  } 
  ////slimmedJetsPuppi
  if(_PuppiVar){
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
      Jet_puppi_pfCombinedMVAV2BJetTags.push_back(j.bDiscriminator("pfCombinedMVAV2BJetTags"));
      Jet_puppi_pfJetProbabilityBJetTags.push_back(j.bDiscriminator("pfJetProbabilityBJetTags"));
      Jet_puppi_pfCombinedCvsLJetTags.push_back(j.bDiscriminator("pfCombinedCvsLJetTags"));
      Jet_puppi_pfCombinedCvsBJetTags.push_back(j.bDiscriminator("pfCombinedCvsBJetTags"));
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
      //Jet Energy Corrections and Uncertainties
      double corrAK4PFPuppi     = 1;
      double corrUpAK4PFPuppi   = 1;
      double corrDownAK4PFPuppi = 1;
      reco::Candidate::LorentzVector uncorrJetAK4PFPuppi = j.correctedP4(0);
      if(!_is_data){
	jecAK4PFPuppiMC_->setJetEta( uncorrJetAK4PFPuppi.eta()    );
	jecAK4PFPuppiMC_->setJetPt ( uncorrJetAK4PFPuppi.pt()     );
	jecAK4PFPuppiMC_->setJetE  ( uncorrJetAK4PFPuppi.energy() );
	jecAK4PFPuppiMC_->setRho	( rho  );
	jecAK4PFPuppiMC_->setNPV	( vertices->size()  );
	jecAK4PFPuppiMC_->setJetA  ( j.jetArea()	     );
	corrAK4PFPuppi = jecAK4PFPuppiMC_->getCorrection();
	jecAK4PFPuppiMCUnc_->setJetEta( uncorrJetAK4PFPuppi.eta() );
	jecAK4PFPuppiMCUnc_->setJetPt( corrAK4PFPuppi * uncorrJetAK4PFPuppi.pt() );
	corrUpAK4PFPuppi = corrAK4PFPuppi * (1 + fabs(jecAK4PFPuppiMCUnc_->getUncertainty(1)));
	jecAK4PFPuppiMCUnc_->setJetEta( uncorrJetAK4PFPuppi.eta() );
	jecAK4PFPuppiMCUnc_->setJetPt( corrAK4PFPuppi * uncorrJetAK4PFPuppi.pt() );
	corrDownAK4PFPuppi = corrAK4PFPuppi * ( 1 - fabs(jecAK4PFPuppiMCUnc_->getUncertainty(-1)) );
      } else {
	jecAK4PFPuppiDATA_->setJetEta( uncorrJetAK4PFPuppi.eta()    );
	jecAK4PFPuppiDATA_->setJetPt ( uncorrJetAK4PFPuppi.pt()     );
	jecAK4PFPuppiDATA_->setJetE  ( uncorrJetAK4PFPuppi.energy() );
	jecAK4PFPuppiDATA_->setRho	( rho  );
	jecAK4PFPuppiDATA_->setNPV	( vertices->size()  );
	jecAK4PFPuppiDATA_->setJetA  ( j.jetArea()	     );
	corrAK4PFPuppi = jecAK4PFPuppiDATA_->getCorrection();
	jecAK4PFPuppiDATAUnc_->setJetEta( uncorrJetAK4PFPuppi.eta() );
	jecAK4PFPuppiDATAUnc_->setJetPt( corrAK4PFPuppi * uncorrJetAK4PFPuppi.pt() );
	corrUpAK4PFPuppi = corrAK4PFPuppi * (1 + fabs(jecAK4PFPuppiDATAUnc_->getUncertainty(1)));
	jecAK4PFPuppiDATAUnc_->setJetEta( uncorrJetAK4PFPuppi.eta() );
	jecAK4PFPuppiDATAUnc_->setJetPt( corrAK4PFPuppi * uncorrJetAK4PFPuppi.pt() );
	corrDownAK4PFPuppi = corrAK4PFPuppi * ( 1 - fabs(jecAK4PFPuppiDATAUnc_->getUncertainty(-1)) );
      }
      Jet_puppi_JesSF.push_back(corrAK4PFPuppi);
      Jet_puppi_JesSFup.push_back(corrUpAK4PFPuppi);
      Jet_puppi_JesSFdown.push_back(corrDownAK4PFPuppi);
      //JER scale factor and uncertainties
      float JERScaleFactor     = 1; 
      float JERScaleFactorUP   = 1;
      float JERScaleFactorDOWN = 1;
      if(!_is_data) GetJER(j, corrAK4PFPuppi, rhoJER, false, JERScaleFactor, JERScaleFactorUP, JERScaleFactorDOWN);
      Jet_puppi_JerSF.push_back(JERScaleFactor);
      Jet_puppi_JerSFup.push_back(JERScaleFactorUP);
      Jet_puppi_JerSFdown.push_back(JERScaleFactorDOWN);
      //delete jecUnc;
      //MC
      if(!_is_data) {
	Jet_puppi_partonFlavour.push_back(j.partonFlavour());
	Jet_puppi_hadronFlavour.push_back(j.hadronFlavour());
      } 
    }
  }
}
void JetSelector::JECInitialization(){
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
  jecPayloadNamesAK4PFchsDATA_.push_back(jecPayloadNamesAK4PFchsDATA4_.fullPath());
  std::vector<JetCorrectorParameters> vParAK4PFchsDATA;
  for ( std::vector<std::string>::const_iterator payloadBegin = jecPayloadNamesAK4PFchsDATA_.begin(),
	  payloadEnd = jecPayloadNamesAK4PFchsDATA_.end(), ipayload = payloadBegin; ipayload != payloadEnd; ++ipayload ) {
    JetCorrectorParameters pars(*ipayload);
    vParAK4PFchsDATA.push_back(pars);
  }
  jecAK4PFchsDATA_    = boost::shared_ptr<FactorizedJetCorrector>  ( new FactorizedJetCorrector(vParAK4PFchsDATA) );
  jecAK4PFchsDATAUnc_ = boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(jecPayloadNamesAK4PFchsDATAUnc_.fullPath()) );
  //AK4Puppi - MC: Get the factorized jet corrector parameters. 
  std::vector<std::string> jecPayloadNamesAK4PFPuppiMC_;
  jecPayloadNamesAK4PFPuppiMC_.push_back(jecPayloadNamesAK4PFPuppiMC1_.fullPath());
  jecPayloadNamesAK4PFPuppiMC_.push_back(jecPayloadNamesAK4PFPuppiMC2_.fullPath());
  jecPayloadNamesAK4PFPuppiMC_.push_back(jecPayloadNamesAK4PFPuppiMC3_.fullPath());
  std::vector<JetCorrectorParameters> vParAK4PFPuppiMC;
  for ( std::vector<std::string>::const_iterator payloadBegin = jecPayloadNamesAK4PFPuppiMC_.begin(),
	  payloadEnd = jecPayloadNamesAK4PFPuppiMC_.end(), ipayload = payloadBegin; ipayload != payloadEnd; ++ipayload ) {
    JetCorrectorParameters pars(*ipayload);
    vParAK4PFPuppiMC.push_back(pars);
  }
  jecAK4PFPuppiMC_    = boost::shared_ptr<FactorizedJetCorrector>  ( new FactorizedJetCorrector(vParAK4PFPuppiMC) );
  jecAK4PFPuppiMCUnc_ = boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(jecPayloadNamesAK4PFPuppiMCUnc_.fullPath()) );
  //AK4Puppi - DATA: Get the factorized jet corrector parameters. 
  std::vector<std::string> jecPayloadNamesAK4PFPuppiDATA_;
  jecPayloadNamesAK4PFPuppiDATA_.push_back(jecPayloadNamesAK4PFPuppiDATA1_.fullPath());
  jecPayloadNamesAK4PFPuppiDATA_.push_back(jecPayloadNamesAK4PFPuppiDATA2_.fullPath());
  jecPayloadNamesAK4PFPuppiDATA_.push_back(jecPayloadNamesAK4PFPuppiDATA3_.fullPath());
  jecPayloadNamesAK4PFPuppiDATA_.push_back(jecPayloadNamesAK4PFPuppiDATA4_.fullPath());
  std::vector<JetCorrectorParameters> vParAK4PFPuppiDATA;
  for ( std::vector<std::string>::const_iterator payloadBegin = jecPayloadNamesAK4PFPuppiDATA_.begin(),
	  payloadEnd = jecPayloadNamesAK4PFPuppiDATA_.end(), ipayload = payloadBegin; ipayload != payloadEnd; ++ipayload ) {
    JetCorrectorParameters pars(*ipayload);
    vParAK4PFPuppiDATA.push_back(pars);
  }
  jecAK4PFPuppiDATA_    = boost::shared_ptr<FactorizedJetCorrector>  ( new FactorizedJetCorrector(vParAK4PFPuppiDATA) );
  jecAK4PFPuppiDATAUnc_ = boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(jecPayloadNamesAK4PFPuppiDATAUnc_.fullPath()) );
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
  AddBranch(&Jet_genpt     ,"Jet_genpt");
  AddBranch(&Jet_geneta    ,"Jet_geneta");
  AddBranch(&Jet_genphi    ,"Jet_genphi");
  AddBranch(&Jet_genenergy ,"Jet_genenergy");
  AddBranch(&Jet_genmass   ,"Jet_genmass");
  AddBranch(&Jet_genpx     ,"Jet_genpx");
  AddBranch(&Jet_genpy     ,"Jet_genpy");
  AddBranch(&Jet_genpz     ,"Jet_genpz");
  //ID
  AddBranch(&Jet_newpfCombinedInclusiveSecondaryVertexV2BJetTags ,"Jet_newpfCombinedInclusiveSecondaryVertexV2BJetTags");
  AddBranch(&Jet_newpfCombinedMVAV2BJetTags                      ,"Jet_newpfCombinedMVAV2BJetTags");
  AddBranch(&Jet_newpfJetProbabilityBJetTags                     ,"Jet_newpfJetProbabilityBJetTags");
  AddBranch(&Jet_newpfCombinedCvsLJetTags                        ,"Jet_newpfCombinedCvsLJetTags");
  AddBranch(&Jet_newpfCombinedCvsBJetTags                        ,"Jet_newpfCombinedCvsBJetTags");
  AddBranch(&Jet_pfCombinedInclusiveSecondaryVertexV2BJetTags ,"Jet_pfCombinedInclusiveSecondaryVertexV2BJetTags");
  AddBranch(&Jet_pfCombinedMVAV2BJetTags                      ,"Jet_pfCombinedMVAV2BJetTags");
  AddBranch(&Jet_pfJetProbabilityBJetTags                     ,"Jet_pfJetProbabilityBJetTags");
  AddBranch(&Jet_pfCombinedCvsLJetTags                        ,"Jet_pfCombinedCvsLJetTags");
  AddBranch(&Jet_pfCombinedCvsBJetTags                        ,"Jet_pfCombinedCvsBJetTags");
  AddBranch(&Jet_pileupId                                     ,"Jet_pileupId");
  AddBranch(&Jet_isPFJet                                      ,"Jet_isPFJet");
  AddBranch(&Jet_isCaloJet                                    ,"Jet_isCaloJet");
  AddBranch(&Jet_qg               ,"Jet_qg");
  AddBranch(&Jet_axis2            ,"Jet_axis2");
  AddBranch(&Jet_ptD              ,"Jet_ptD");
  AddBranch(&Jet_mult             ,"Jet_mult");
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
  //Jet Energy Corrections and Uncertainties
  AddBranch(&Jet_JesSF                ,"Jet_JesSF");
  AddBranch(&Jet_JesSFup              ,"Jet_JesSFup");
  AddBranch(&Jet_JesSFdown            ,"Jet_JesSFdown");
  AddBranch(&Jet_JerSF                ,"Jet_JerSF");
  AddBranch(&Jet_JerSFup              ,"Jet_JerSFup");
  AddBranch(&Jet_JerSFdown            ,"Jet_JerSFdown");
  //MC
  if(!_is_data) {
    AddBranch(&Jet_partonFlavour        ,"Jet_partonFlavour");
    AddBranch(&Jet_hadronFlavour        ,"Jet_hadronFlavour");
  }
  ////slimmedJetsPuppi
  if(_PuppiVar){
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
    AddBranch(&Jet_puppi_pfCombinedMVAV2BJetTags                      ,"Jet_puppi_pfCombinedMVAV2BJetTags");
    AddBranch(&Jet_puppi_pfJetProbabilityBJetTags                     ,"Jet_puppi_pfJetProbabilityBJetTags");
    AddBranch(&Jet_puppi_pfCombinedCvsLJetTags                        ,"Jet_puppi_pfCombinedCvsLJetTags");
    AddBranch(&Jet_puppi_pfCombinedCvsBJetTags                        ,"Jet_puppi_pfCombinedCvsBJetTags");
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
    //Jet Energy Corrections and Uncertainties
    AddBranch(&Jet_puppi_JesSF                ,"Jet_puppi_JesSF");
    AddBranch(&Jet_puppi_JesSFup              ,"Jet_puppi_JesSFup");
    AddBranch(&Jet_puppi_JesSFdown            ,"Jet_puppi_JesSFdown");
    AddBranch(&Jet_puppi_JerSF                ,"Jet_puppi_JerSF");
    AddBranch(&Jet_puppi_JerSFup              ,"Jet_puppi_JerSFup");
    AddBranch(&Jet_puppi_JerSFdown            ,"Jet_puppi_JerSFdown");
    //MC
    if(!_is_data) {
      AddBranch(&Jet_puppi_partonFlavour        ,"Jet_puppi_partonFlavour");
      AddBranch(&Jet_puppi_hadronFlavour        ,"Jet_puppi_hadronFlavour");
    }
  }
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
  Jet_genpt.clear();
  Jet_geneta.clear();
  Jet_genphi.clear();
  Jet_genenergy.clear();
  Jet_genmass.clear();
  Jet_genpx.clear();
  Jet_genpy.clear();
  Jet_genpz.clear();
  //ID
  Jet_newpfCombinedInclusiveSecondaryVertexV2BJetTags.clear();
  Jet_newpfCombinedMVAV2BJetTags.clear();
  Jet_newpfJetProbabilityBJetTags.clear();
  Jet_newpfCombinedCvsLJetTags.clear();
  Jet_newpfCombinedCvsBJetTags.clear();
  Jet_pfCombinedInclusiveSecondaryVertexV2BJetTags.clear();
  Jet_pfCombinedMVAV2BJetTags.clear();
  Jet_pfJetProbabilityBJetTags.clear();
  Jet_pfCombinedCvsLJetTags.clear();
  Jet_pfCombinedCvsBJetTags.clear();
  Jet_pileupId.clear();
  Jet_isPFJet.clear();
  Jet_isCaloJet.clear();
  Jet_qg.clear();
  Jet_axis2.clear();
  Jet_ptD.clear();
  Jet_mult.clear();
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
  //Jet Energy Corrections and Uncertainties
  Jet_JesSF.clear();
  Jet_JesSFup.clear();
  Jet_JesSFdown.clear();
  Jet_JerSF.clear();
  Jet_JerSFup.clear();
  Jet_JerSFdown.clear(); 
  //MC
  if(!_is_data) {
    Jet_partonFlavour.clear();
    Jet_hadronFlavour.clear();
  }
  ////slimmedJetsPuppi
  if(_PuppiVar){
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
    Jet_puppi_pfCombinedMVAV2BJetTags.clear();
    Jet_puppi_pfJetProbabilityBJetTags.clear();
    Jet_puppi_pfCombinedCvsLJetTags.clear();
    Jet_puppi_pfCombinedCvsBJetTags.clear();
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
    Jet_puppi_JesSF.clear();
    Jet_puppi_JesSFup.clear();
    Jet_puppi_JesSFdown.clear();
    Jet_puppi_JerSF.clear();
    Jet_puppi_JerSFup.clear();
    Jet_puppi_JerSFdown.clear(); 
    //MC
    if(!_is_data) {
      Jet_puppi_partonFlavour.clear();
      Jet_puppi_hadronFlavour.clear();
    }
  }
}
void JetSelector::GetJER(pat::Jet jet, float JesSF, float rhoJER, bool AK4PFchs, float &JERScaleFactor, float &JERScaleFactorUP, float &JERScaleFactorDOWN){
  if(!jet.genJet()) return;
  double jetEta=fabs(jet.eta());
  double cFactorJER = 1.0; 
  double cFactorJERdown = 1.0;
  double cFactorJERup = 1.0;
  //https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution#JER_Scaling_factors_and_Unce_AN1
  if( jetEta<0.5 ){ 
    cFactorJER = 1.122; 
    cFactorJERdown = 1.122-0.026;
    cFactorJERup   = 1.122+0.026; 
  } else if( jetEta<0.8 ){ 
    cFactorJER = 1.167; 
    cFactorJERdown = 1.167-0.048;
    cFactorJERup   = 1.167+0.048; 
  } else if( jetEta<1.1 ){ 
    cFactorJER = 1.168; 
    cFactorJERdown = 1.168-0.046;
    cFactorJERup   = 1.168+0.046; 
  } else if( jetEta<1.3 ){ 
    cFactorJER = 1.029; 
    cFactorJERdown = 1.029-0.066;
    cFactorJERup   = 1.029+0.066; 
  } else if( jetEta<1.7 ){ 
    cFactorJER = 1.115; 
    cFactorJERdown = 1.115-0.030;
    cFactorJERup   = 1.115+0.030; 
  } else if( jetEta<1.9 ){ 
    cFactorJER = 1.041; 
    cFactorJERdown = 1.041-0.062;
    cFactorJERup   = 1.041+0.062; 
  } else if( jetEta<2.1 ){ 
    cFactorJER = 1.167; 
    cFactorJERdown = 1.167-0.086;
    cFactorJERup   = 1.167+0.086; 
  } else if( jetEta<2.3 ){ 
    cFactorJER = 1.094; 
    cFactorJERdown = 1.094-0.093;
    cFactorJERup   = 1.094+0.093; 
  } else if( jetEta<2.5 ){ 
    cFactorJER = 1.168; 
    cFactorJERdown = 1.168-0.120;
    cFactorJERup   = 1.168+0.120; 
  } else if( jetEta<2.8 ){ 
    cFactorJER = 1.266; 
    cFactorJERdown = 1.266-0.132;
    cFactorJERup   = 1.266+0.132; 
  } else if( jetEta<3.0 ){ 
    cFactorJER = 1.595; 
    cFactorJERdown = 1.595-0.175;
    cFactorJERup   = 1.595+0.175; 
  } else if( jetEta<3.2 ){ 
    cFactorJER = 0.998; 
    cFactorJERdown = 0.998-0.066;
    cFactorJERup   = 0.998+0.066; 
  } else if( jetEta<5.0 ){ 
    cFactorJER = 1.226; 
    cFactorJERdown = 1.226-0.145;
    cFactorJERup   = 1.226+0.145;
  }
  //double recoJetPt = jet.pt();//(jet.correctedJet("Uncorrected").pt())*JesSF;
  double recoJetPt = (jet.correctedJet("Uncorrected").pt())*JesSF;
  double genJetPt  = jet.genJet()->pt();
  double diffPt    = recoJetPt - genJetPt;
  JME::JetResolution resolution;
  JME::JetResolutionScaleFactor res_sf;
  if(AK4PFchs){
    resolution = JME::JetResolution(jerAK4PFchs_);
    res_sf = JME::JetResolutionScaleFactor(jerAK4PFchsSF_);
  } else {
    resolution = JME::JetResolution(jerAK4PFPuppi_);
    res_sf = JME::JetResolutionScaleFactor(jerAK4PFPuppiSF_);
  }
  JME::JetParameters parameters;
  parameters.setJetPt(jet.pt());
  parameters.setJetEta(jet.eta());
  parameters.setRho(rhoJER);
  float relpterr = resolution.getResolution(parameters);
  if(genJetPt>0. && deltaR(jet.eta(),jet.phi(),jet.genJet()->eta(),jet.genJet()->phi())<0.2
     && (abs(jet.pt()-jet.genJet()->pt())<3*relpterr*jet.pt())) {
    JERScaleFactor     = (std::max(0., genJetPt + cFactorJER*diffPt))/recoJetPt;
    JERScaleFactorUP   = (std::max(0., genJetPt + cFactorJERup*diffPt))/recoJetPt;
    JERScaleFactorDOWN = (std::max(0., genJetPt + cFactorJERdown*diffPt))/recoJetPt;
  } else {
    JERScaleFactor     = 1.;
    JERScaleFactorUP   = 1.;
    JERScaleFactorDOWN = 1.;
  } 
}
