//Please note that 
// else if(jetPt >=160 && jetPt<100000) iPt = 5; was originally else if(jetPt >=160 && jetPt<10000) iPt = 5;
// but it could crash if jetPt>10000, so I temporarily changed the value to 100000
// I do not know what it changes in the analysis, but since we are not really use this class right now, I think it is temporarily ok
#include "BSMFramework/BSM3G_TNT_Maker/interface/BTagReweight.h"
BTagReweight::BTagReweight(std::string name, TTree* tree, bool debug, const pset& iConfig, edm::ConsumesCollector && ic):
  baseTree(name,tree,debug),
  eleMVATrigIdMapToken_(ic.consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleMVATrigIdMap")))
{
  if(debug) std::cout<<"in BTagReweight constructor"<<std::endl;
  if(debug) std::cout<<"in BTagReweight constructor: calling SetBrances()"<<std::endl;
  vtx_h_               = ic.consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"));
  electron_pat_        = ic.consumes<edm::View<pat::Electron> >(iConfig.getParameter<edm::InputTag>("patElectrons"));
  muon_h_              = ic.consumes<edm::View<pat::Muon> >(iConfig.getParameter<edm::InputTag>("muons"));
  jets_                = ic.consumes<pat::JetCollection >(iConfig.getParameter<edm::InputTag>("jets"));
  rhopogHandle_        = ic.consumes<double>(edm::InputTag("fixedGridRhoFastjetAll"));
  BTAGReweightfile1_ = iConfig.getParameter<edm::FileInPath>("BTAGReweightfile1");
  BTAGReweightfile2_ = iConfig.getParameter<edm::FileInPath>("BTAGReweightfile2");
  _vtx_ndof_min       = iConfig.getParameter<int>("vtx_ndof_min");
  _vtx_rho_max        = iConfig.getParameter<int>("vtx_rho_max");
  _vtx_position_z_max = iConfig.getParameter<double>("vtx_position_z_max");
  _is_data = iConfig.getParameter<bool>("is_data");
  const char *filePathHF = BTAGReweightfile1_.fullPath().c_str();
  const char *filePathLF = BTAGReweightfile2_.fullPath().c_str();
  TFile* f_CSVwgt_HF = new TFile (filePathHF);
  TFile* f_CSVwgt_LF = new TFile (filePathLF);
  fillCSVhistos(f_CSVwgt_HF, f_CSVwgt_LF);	
  SetBranches();
}
BTagReweight::~BTagReweight(){
  delete tree_;
}

void BTagReweight::Fill(const edm::Event& iEvent){
  if(debug_) std::cout<<"getting BTagReweight info"<<std::endl;
  if(!_is_data){
    edm::Handle<pat::JetCollection> jets;
    iEvent.getByToken(jets_, jets);
    std::vector<double> jetPts;
    std::vector<double> jetEtas;
    std::vector<double> jetCSVs;
    std::vector<int> jetFlavors;
    vector<TLorentzVector> LeptonsForDeltaRWithJets;
    GetLeptonsForDeltaRWithJets(LeptonsForDeltaRWithJets,iEvent);
    for(const pat::Jet &j : *jets){ 
      float JERScaleFactor     = 1; 
      float JERScaleFactorUP   = 1;
      float JERScaleFactorDOWN = 1;
      GetJER(j, 1, JERScaleFactor, JERScaleFactorUP, JERScaleFactorDOWN);
      if((j.pt()*JERScaleFactor) < 20) continue;
      if(LeptonsForDeltaRWithJets.size()==1){if((j.pt()*JERScaleFactor) < 30) continue;}
      if(fabs(j.eta()) > 2.4) continue;
      if(!(j.neutralHadronEnergyFraction()<0.99)) continue;
      if(!(j.chargedEmEnergyFraction()<0.99)) continue;
      if(!(j.neutralEmEnergyFraction()<0.99)) continue;
      if(!((j.chargedMultiplicity() + j.neutralMultiplicity())>1)) continue;
      if(!(j.chargedHadronEnergyFraction()>0.0)) continue;
      if(!(j.chargedMultiplicity()>0.0)) continue;
      bool deltaRJetLepBoolean = true;
      for (size_t k = 0; k < LeptonsForDeltaRWithJets.size(); ++k){
        float deltaEta = LeptonsForDeltaRWithJets[k].Eta()-j.eta();
        float deltaPhi = fabs(LeptonsForDeltaRWithJets[k].Phi()-j.phi());
        if(deltaPhi>TMath::Pi()) deltaPhi = 2*TMath::Pi() - deltaPhi;
        if(sqrt(deltaEta*deltaEta + deltaPhi*deltaPhi)<0.4) deltaRJetLepBoolean=false;
      }
      if(!(deltaRJetLepBoolean==true)) continue;
      jetPts.push_back(j.pt()*JERScaleFactor);
      jetEtas.push_back(j.eta());
      if(j.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags")<1.0) jetCSVs.push_back(j.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
      else jetCSVs.push_back(1);
      //jetFlavors.push_back(j.partonFlavour());
      jetFlavors.push_back(j.hadronFlavour());
    }
    int iSys = 0;
    double wgt_csv_hf, wgt_csv_lf, wgt_csv_cf;
    double wgt_csv = get_csv_wgt(jetPts, jetEtas, jetCSVs, jetFlavors, iSys, wgt_csv_hf, wgt_csv_lf, wgt_csv_cf);	
    bWeight = wgt_csv;
  } else {	
    bWeight = 1.;
  }
  if(debug_) std::cout<<"got BTagReweight info"<<std::endl;
}

void BTagReweight::SetBranches(){
  if(debug_) std::cout<<"setting branches: calling AddBranch of BTagReweight"<<std::endl;
  AddBranch(&bWeight,"bWeight");
}

//Fill the histograms (done once)
void BTagReweight::fillCSVhistos(TFile* fileHF, TFile* fileLF){
  for(int iSys=0; iSys<9; iSys++){
    for(int iPt=0; iPt<5; iPt++) h_csv_wgt_hf[iSys][iPt] = NULL;
    for(int iPt=0; iPt<3; iPt++){
      for(int iEta=0; iEta<3; iEta++)h_csv_wgt_lf[iSys][iPt][iEta] = NULL;
    }
  }
  for(int iSys=0; iSys<5; iSys++){
    for(int iPt=0; iPt<5; iPt++) c_csv_wgt_hf[iSys][iPt] = NULL;
  }
  //CSV reweighting /// only care about the nominal ones
  for(int iSys=0; iSys<9; iSys++){
    TString syst_csv_suffix_hf = "final";
    TString syst_csv_suffix_c = "final";
    TString syst_csv_suffix_lf = "final";
    switch(iSys){
    case 0:
      //this is the nominal case
      break;
    case 1:
      //JESUp
      syst_csv_suffix_hf = "final_JESUp"; syst_csv_suffix_lf = "final_JESUp";
      syst_csv_suffix_c  = "final_cErr1Up";
      break;
    case 2:
      //JESDown
      syst_csv_suffix_hf = "final_JESDown"; syst_csv_suffix_lf = "final_JESDown";
      syst_csv_suffix_c  = "final_cErr1Down";
      break;
    case 3:
      //purity up
      syst_csv_suffix_hf = "final_LFUp"; syst_csv_suffix_lf = "final_HFUp";
      syst_csv_suffix_c  = "final_cErr2Up";
      break;
    case 4:
      //purity down
      syst_csv_suffix_hf = "final_LFDown"; syst_csv_suffix_lf = "final_HFDown";
      syst_csv_suffix_c  = "final_cErr2Down";
      break;
    case 5:
      //stats1 up
      syst_csv_suffix_hf = "final_Stats1Up"; syst_csv_suffix_lf = "final_Stats1Up";
      break;
    case 6:
      //stats1 down
      syst_csv_suffix_hf = "final_Stats1Down"; syst_csv_suffix_lf = "final_Stats1Down";
      break;
    case 7:
      //stats2 up
      syst_csv_suffix_hf = "final_Stats2Up"; syst_csv_suffix_lf = "final_Stats2Up";
      break;
    case 8:
      //stats2 down
      syst_csv_suffix_hf = "final_Stats2Down"; syst_csv_suffix_lf = "final_Stats2Down";
      break;
    }
    for(int iPt=0; iPt<5; iPt++) h_csv_wgt_hf[iSys][iPt] = (TH1D*)fileHF->Get( Form("csv_ratio_Pt%i_Eta0_%s",iPt,syst_csv_suffix_hf.Data()) );
    if(iSys<5){
      for(int iPt=0; iPt<5; iPt++) c_csv_wgt_hf[iSys][iPt] = (TH1D*)fileHF->Get( Form("c_csv_ratio_Pt%i_Eta0_%s",iPt,syst_csv_suffix_c.Data()) );
    }
    for(int iPt=0; iPt<4; iPt++){
      for(int iEta=0; iEta<3; iEta++ )h_csv_wgt_lf[iSys][iPt][iEta] = (TH1D*)fileLF->Get( Form("csv_ratio_Pt%i_Eta%i_%s",iPt,iEta,syst_csv_suffix_lf.Data()) );
    }
  }
  return;
}

double BTagReweight::get_csv_wgt(std::vector<double> jetPts, std::vector<double> jetEtas, std::vector<double> jetCSVs, std::vector<int> jetFlavors, int iSys, double &csvWgtHF, double &csvWgtLF, double &csvWgtCF ){
  int iSysHF = 0;
  switch(iSys){
  case 7:  iSysHF=1; break; //JESUp
  case 8:  iSysHF=2; break; //JESDown
  case 9:  iSysHF=3; break; //LFUp
  case 10: iSysHF=4; break; //LFDown
  case 13: iSysHF=5; break; //Stats1Up
  case 14: iSysHF=6; break; //Stats1Down
  case 15: iSysHF=7; break; //Stats2Up
  case 16: iSysHF=8; break; //Stats2Down
  default : iSysHF = 0; break; //NoSys
  }
  int iSysC = 0;
  switch(iSys){
  case 21: iSysC=1; break;
  case 22: iSysC=2; break;
  case 23: iSysC=3; break;
  case 24: iSysC=4; break;
  default : iSysC = 0; break;
  }
  int iSysLF = 0;
  switch(iSys){
  case 7:  iSysLF=1; break; //JESUp
  case 8:  iSysLF=2; break; //JESDown
  case 11: iSysLF=3; break; //HFUp
  case 12: iSysLF=4; break; //HFDown
  case 17: iSysLF=5; break; //Stats1Up
  case 18: iSysLF=6; break; //Stats1Down
  case 19: iSysLF=7; break; //Stats2Up
  case 20: iSysLF=8; break; //Stats2Down
  default : iSysLF = 0; break; //NoSys
  }
  double csvWgthf = 1.;
  double csvWgtC  = 1.;
  double csvWgtlf = 1.;
  for(int iJet=0; iJet<int(jetPts.size()); iJet++){
    double csv = jetCSVs[iJet];
    double jetPt = jetPts[iJet];
    double jetAbsEta = fabs(jetEtas[iJet]);
    int flavor = jetFlavors[iJet];
    int iPt = -1; int iEta = -1;
    if (jetPt >=19.99 && jetPt<30) iPt = 0;
    else if (jetPt >=30 && jetPt<40) iPt = 1;
    else if (jetPt >=40 && jetPt<60) iPt = 2;
    else if (jetPt >=60 && jetPt<100) iPt = 3;
    else if (jetPt >=100) iPt = 4;
    if (jetAbsEta >=0 &&  jetAbsEta<0.8) iEta = 0;
    else if (jetAbsEta>=0.8 && jetAbsEta<1.6)  iEta = 1;
    else if (jetAbsEta>=1.6 && jetAbsEta<2.41) iEta = 2;
    if (iPt < 0 || iEta < 0) std::cout << "Error, couldn't find Pt, Eta bins for this b-flavor jet, jetPt = " << jetPt << ", jetAbsEta = " << jetAbsEta << std::endl;
    if (abs(flavor) == 5){
      int useCSVBin = (csv>=0.) ? h_csv_wgt_hf[iSysHF][iPt]->FindBin(csv) : 1;
      double iCSVWgtHF = h_csv_wgt_hf[iSysHF][iPt]->GetBinContent(useCSVBin);
      if(iCSVWgtHF!=0) csvWgthf *= iCSVWgtHF;
    }
    else if(abs(flavor) == 4){
      int useCSVBin = (csv>=0.) ? c_csv_wgt_hf[iSysC][iPt]->FindBin(csv) : 1;
      double iCSVWgtC = c_csv_wgt_hf[iSysC][iPt]->GetBinContent(useCSVBin);
      if(iCSVWgtC!=0) csvWgtC *= iCSVWgtC;
    }
    else {
      if (iPt >=3) iPt=3;       /// [30-40], [40-60] and [60-10000] only 3 Pt bins for lf
      int useCSVBin = (csv>=0.) ? h_csv_wgt_lf[iSysLF][iPt][iEta]->FindBin(csv) : 1;
      double iCSVWgtLF = h_csv_wgt_lf[iSysLF][iPt][iEta]->GetBinContent(useCSVBin);
      if(iCSVWgtLF!=0) csvWgtlf *= iCSVWgtLF;
    }
  }
  double csvWgtTotal = csvWgthf * csvWgtC * csvWgtlf;
  csvWgtHF = csvWgthf;
  csvWgtLF = csvWgtlf;
  csvWgtCF = csvWgtC;
  return csvWgtTotal;
}
void BTagReweight::GetJER(pat::Jet jet, float JesSF, float &JERScaleFactor, float &JERScaleFactorUP, float &JERScaleFactorDOWN){
  if(!jet.genJet()) return;
  double jetEta=fabs(jet.eta());
  double cFactorJER = 1.0; 
  double cFactorJERdown = 1.0;
  double cFactorJERup = 1.0;
  //https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution#JER_Scaling_factors_and_Unce_AN1
  string ERA="13TeV";
  if(ERA=="8TeV"){
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
  } else if(ERA=="13TeV"){
    if( jetEta<0.8 ){ 
      cFactorJER = 1.061; 
      cFactorJERdown = 1.061-0.023;
      cFactorJERup   = 1.061+0.023; 
    }
    else if( jetEta<1.3 ){ 
      cFactorJER = 1.088; 
      cFactorJERdown = 1.088-0.029;
      cFactorJERup   = 1.088+0.029; 
    }
    else if( jetEta<1.9 ){ 
      cFactorJER = 1.106; 
      cFactorJERdown = 1.106-0.030;
      cFactorJERup   = 1.106+0.030; 
    }
    else if( jetEta<2.5 ){ 
      cFactorJER = 1.126; 
      cFactorJERdown = 1.126-0.094;
      cFactorJERup   = 1.126+0.094; 
    }
    else if( jetEta<3.0 ){ 
      cFactorJER = 1.343; 
      cFactorJERdown = 1.343-0.123;
      cFactorJERup   = 1.343+0.123; 
    }
    else if( jetEta<3.2 ){ 
      cFactorJER = 1.303; 
      cFactorJERdown = 1.303-0.111;
      cFactorJERup   = 1.303+0.111; 
    }
    else if( jetEta<5.0 ){ 
      cFactorJER = 1.320; 
      cFactorJERdown = 1.320-0.286;
      cFactorJERup   = 1.320+0.286; 
    }
  }
  double recoJetPt = jet.pt();//(jet.correctedJet("Uncorrected").pt())*JesSF;
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

void BTagReweight::GetLeptonsForDeltaRWithJets(vector<TLorentzVector> &LeptonsForDeltaRWithJets, const edm::Event& iEvent){
  edm::Handle<reco::VertexCollection> vtx_h;
  iEvent.getByToken(vtx_h_, vtx_h);
  edm::Handle<edm::View<pat::Electron> > electron_pat;
  iEvent.getByToken(electron_pat_, electron_pat);
  edm::Handle<edm::View<pat::Muon> > muon_h;
  iEvent.getByToken(muon_h_, muon_h);
  edm::Handle<pat::JetCollection> jets;
  iEvent.getByToken(jets_, jets);
  edm::Handle<double> rhopogHandle;
  iEvent.getByToken(rhopogHandle_,rhopogHandle);
  double rhopog = *rhopogHandle;
  edm::Handle<edm::ValueMap<bool>  > mvatrig_id_decisions;
  iEvent.getByToken(eleMVATrigIdMapToken_,mvatrig_id_decisions);
  if(vtx_h->empty()) return; // skip the event if no PV found
  const reco::Vertex &firstGoodVertex = vtx_h->front();  
  bool isgoodvtx = isGoodVertex(firstGoodVertex);
  if(!isgoodvtx) return;
  //LEPTON SELECTION - MUON
  for(edm::View<pat::Muon>::const_iterator mu = muon_h->begin(); mu != muon_h->end(); mu++){
    double SumChHadPt  = mu->pfIsolationR04().sumChargedHadronPt;
    double SumNeuHadEt = mu->pfIsolationR04().sumNeutralHadronEt;
    double SumPhotonEt = mu->pfIsolationR04().sumPhotonEt;
    double SumPU       = mu->pfIsolationR04().sumPUPt;
    double SumNeutralCorrEt = std::max( 0.0, SumNeuHadEt+SumPhotonEt - 0.5*SumPU );
    double relIsoDeltaBeta = (SumChHadPt + SumNeutralCorrEt)/mu->pt();
    if(!(mu->pt()>15))                         continue;
    if(!(fabs(mu->eta())<2.4))                 continue;  
    if(!(mu->isTightMuon(firstGoodVertex)==1)) continue;
    if(!(relIsoDeltaBeta<0.25))                continue;
    TLorentzVector muon = TLorentzVector(mu->px(),mu->py(),mu->pz(),mu->p4().E());
    LeptonsForDeltaRWithJets.push_back(muon);
  }
  //LEPTON SELECTION - ELECTRON (tth)
  for(edm::View<pat::Electron>::const_iterator el = electron_pat->begin(); el != electron_pat->end(); el++){
    const Ptr<pat::Electron> elPtr(electron_pat, el - electron_pat->begin() );
    bool isPassMvatrig = (*mvatrig_id_decisions) [ elPtr ];
    double EleSCeta    = el->superCluster()->position().eta();
    double SumChHadPt  = el->pfIsolationVariables().sumChargedHadronPt;
    double SumNeuHadEt = el->pfIsolationVariables().sumNeutralHadronEt;
    double SumPhotonEt = el->pfIsolationVariables().sumPhotonEt; 
    double EffArea = get_effarea(el->superCluster()->position().eta());
    double SumNeutralCorrEt = std::max( 0.0, SumNeuHadEt+SumPhotonEt - rhopog*EffArea );
    double relIsoRhoEA = (SumChHadPt + SumNeutralCorrEt)/el->pt();
    if(!(el->pt()>15))                                   continue;
    if(!(fabs(el->eta())<2.4))                           continue; 
    if((fabs(EleSCeta)>1.4442 && fabs(EleSCeta)<1.5660)) continue; 
    if(!(isPassMvatrig==1))                              continue;
    if(!(relIsoRhoEA<0.15))                              continue;
    if(fabs(EleSCeta)<1.4442){
      if(!(el->full5x5_sigmaIetaIeta()<0.012))                 continue;
      if(!(el->hcalOverEcal()<0.09))                           continue;
      if(!((el->ecalPFClusterIso()/el->pt())<0.37))            continue;
      if(!((el->hcalPFClusterIso()/el->pt())<0.25))            continue;
      if(!((el->dr03TkSumPt()/el->pt())<0.18))                 continue;
      if(!(fabs(el->deltaEtaSuperClusterTrackAtVtx())<0.0095)) continue;
      if(!(fabs(el->deltaPhiSuperClusterTrackAtVtx())<0.065))  continue;
    }
    if(fabs(EleSCeta)>1.5660){
      if(!(el->full5x5_sigmaIetaIeta()<0.033))      continue;
      if(!(el->hcalOverEcal()<0.09))                continue;
      if(!((el->ecalPFClusterIso()/el->pt())<0.45)) continue;
      if(!((el->hcalPFClusterIso()/el->pt())<0.28)) continue;
      if(!((el->dr03TkSumPt()/el->pt())<0.18))      continue;
    }
    TLorentzVector electron = TLorentzVector(el->px(),el->py(),el->pz(),el->p4().E());
    LeptonsForDeltaRWithJets.push_back(electron);
  }
}

double BTagReweight::get_effarea(double eta){
  double effarea = -1;
  if(abs(eta) < 1.0)        effarea = 0.1752;
  else if(abs(eta) < 1.479) effarea = 0.1862;
  else if(abs(eta) < 2.0)   effarea = 0.1411;
  else if(abs(eta) < 2.2)   effarea = 0.1534;
  else if(abs(eta) < 2.3)   effarea = 0.1903;
  else if(abs(eta) < 2.4)   effarea = 0.2243;
  else                      effarea = 0.2687;
  return effarea;
}

bool BTagReweight::isGoodVertex(const reco::Vertex& vtx){
  if(vtx.isFake())                                   return false;
  if(vtx.ndof()<_vtx_ndof_min)                       return false;
  if(vtx.position().Rho()>_vtx_rho_max)              return false;
  if(fabs(vtx.position().Z()) > _vtx_position_z_max) return false;
  return true;
}
