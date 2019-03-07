#include "BSMFramework/BSM3G_TNT_Maker/interface/EventInfoSelector.h"
EventInfoSelector::EventInfoSelector(std::string name, TTree* tree, bool debug, const pset& iConfig, edm::ConsumesCollector && ic):
  baseTree(name,tree,debug)
{
  genEvtInfo_    = ic.consumes<GenEventInfoProduct>(edm::InputTag("generator"));
  lheEventProduct_ = ic.consumes<LHEEventProduct>(edm::InputTag("externalLHEProducer"));
  rhopogHandle_  = ic.consumes<double>(edm::InputTag("fixedGridRhoFastjetAll"));
  rhotthHandle_  = ic.consumes<double>(edm::InputTag("fixedGridRhoFastjetCentralNeutral"));
  fixedGridRhoFastjetCentralHandle_  = ic.consumes<double>(edm::InputTag("fixedGridRhoFastjetCentral"));
  fixedGridRhoFastjetCentralChargedPileUpHandle_  = ic.consumes<double>(edm::InputTag("fixedGridRhoFastjetCentralChargedPileUp"));
  fixedGridRhoFastjetCentralNeutralHandle_  = ic.consumes<double>(edm::InputTag("fixedGridRhoFastjetCentralNeutral"));
  metFilterBits_ = ic.consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("metFilterBits"));
  ecalBadCalibFilterUpdate_token = ic.consumes< bool >(edm::InputTag("ecalBadCalibReducedMINIAODFilter"));
  prefweight_token = ic.consumes< double >(edm::InputTag("prefiringweight:NonPrefiringProb"));
  prefweightup_token = ic.consumes< double >(edm::InputTag("prefiringweight:NonPrefiringProbUp"));
  prefweightdown_token = ic.consumes< double >(edm::InputTag("prefiringweight:NonPrefiringProbDown"));
  _is_data = iConfig.getParameter<bool>("is_data");
  _MC2016  = iConfig.getParameter<bool>("MC2016");
  if(debug) std::cout<<"in EventInfoSelector constructor"<<std::endl;
  SetBranches();
  read_PDFSet = new (LHAPDF::PDFSet)("NNPDF30_nlo_as_0118");
  _systPDFs = read_PDFSet->mkPDFs();
}
EventInfoSelector::~EventInfoSelector(){
  delete tree_;
}
void EventInfoSelector::Fill(const edm::Event& iEvent){
  Initialise();
  //Event quantities
  EVENT_event_     = iEvent.id().event();
  EVENT_run_       = iEvent.id().run();
  EVENT_lumiBlock_ = iEvent.id().luminosityBlock();
  EVENT_genWeight_ = 1;
  edm::Handle<GenEventInfoProduct> genEvtInfo;
  iEvent.getByToken(genEvtInfo_,genEvtInfo);
  edm::Handle<LHEEventProduct> lheEventProduct;
  iEvent.getByToken(lheEventProduct_, lheEventProduct);
  if(!_is_data){
    //PreFiring
    edm::Handle< double > theprefweight;
    iEvent.getByToken(prefweight_token, theprefweight ) ;
    EVENT_prefiringweight_ =(*theprefweight);
    edm::Handle< double > theprefweightup;
    iEvent.getByToken(prefweightup_token, theprefweightup ) ;
    EVENT_prefiringweightup_ =(*theprefweightup);
    edm::Handle< double > theprefweightdown;
    iEvent.getByToken(prefweightdown_token, theprefweightdown ) ;
    EVENT_prefiringweightdown_ =(*theprefweightdown);
    //Gen Weights
    EVENT_genWeight_ = genEvtInfo->weight();
    const GenEventInfoProduct& genEventInfoW = *(genEvtInfo.product());
    const gen::PdfInfo* pdf = genEventInfoW.pdf();
    EVENT_scalePDF_ = pdf->scalePDF;
    if(lheEventProduct.isValid()){
      EVENT_originalXWGTUP_ = lheEventProduct->originalXWGTUP();
      for (unsigned int i=0; i<lheEventProduct->weights().size(); i++){
	EVENT_genWeights_.push_back(lheEventProduct->weights()[i].wgt);
	//Q2 for ttHbb synchronization
	if (lheEventProduct->weights()[i].id == "1005") EVENT_Q2tthbbWeightUp_   = lheEventProduct->weights()[i].wgt/lheEventProduct->originalXWGTUP(); 
	if (lheEventProduct->weights()[i].id == "1009") EVENT_Q2tthbbWeightDown_ = lheEventProduct->weights()[i].wgt/lheEventProduct->originalXWGTUP(); 
      }
    } else {
      EVENT_genWeights_.push_back(-99);
      EVENT_originalXWGTUP_    = -99;
      EVENT_Q2tthbbWeightUp_   = -99;
      EVENT_Q2tthbbWeightDown_ = -99;
    }
    //PDF for ttHbb synchronization
    auto pdfInfos = genEvtInfo->pdf();
    double pdfNominal = pdfInfos->xPDF.first * pdfInfos->xPDF.second;
    std::vector<double> pdfs;
    for (size_t j = 0; j < _systPDFs.size(); ++j) {
      double xpdf1 = _systPDFs[j]->xfxQ(pdfInfos->id.first, pdfInfos->x.first, pdfInfos->scalePDF);
      double xpdf2 = _systPDFs[j]->xfxQ(pdfInfos->id.second, pdfInfos->x.second, pdfInfos->scalePDF);
      pdfs.push_back(xpdf1 * xpdf2);
    }
    const LHAPDF::PDFUncertainty pdfUnc = read_PDFSet->uncertainty(pdfs, 68.);
    double weight_up = 1.0;
    double weight_down = 1.0;
    if (std::isfinite(1./pdfNominal)) {
      weight_up = (pdfUnc.central + pdfUnc.errplus) / pdfNominal;
      weight_down = (pdfUnc.central - pdfUnc.errminus) / pdfNominal;
    }
    EVENT_PDFtthbbWeightUp_   = weight_up;
    EVENT_PDFtthbbWeightDown_ = weight_down;
    //gen Parton HT
    //Definition taken from https://github.com/cmkuo/ggAnalysis/blob/a24edc65be23b402d761c75545192ce79cddf316/ggNtuplizer/plugins/ggNtuplizer_genParticles.cc#L201 
    //Zaixing has a somehow different, but likely equivalent implementation
    //https://github.com/zaixingmao/FSA/blob/miniAOD_dev_7_4_14/DataFormats/src/PATFinalStateEvent.cc#L153
    double lheHt = 0.;
    double lheV_pt = 0.;
    bool lCheck = false;
    bool lbarCheck = false;
    bool vlCheck = false;
    bool vlbarCheck = false;
    TLorentzVector l, lbar, vl, vlbar, V_tlv;
    if(lheEventProduct.isValid()){
      const lhef::HEPEUP& lheEvent = lheEventProduct->hepeup();
      std::vector<lhef::HEPEUP::FiveVector> lheParticles = lheEvent.PUP;
      size_t numParticles = lheParticles.size();
      for(size_t idxParticle = 0; idxParticle < numParticles; ++idxParticle){
        int absPdgId = TMath::Abs(lheEvent.IDUP[idxParticle]);
        int id       = lheEvent.IDUP[idxParticle];
        int status = lheEvent.ISTUP[idxParticle];
        if(status == 1 && ((absPdgId >= 1 && absPdgId <= 6) || absPdgId == 21)){ // quarks and gluons
          lheHt += TMath::Sqrt(TMath::Power(lheParticles[idxParticle][0], 2.) + TMath::Power(lheParticles[idxParticle][1], 2.)); // first entry is px, second py
        }
	if(id==11){ l.SetPxPyPzE(lheParticles[idxParticle][0],lheParticles[idxParticle][1],lheParticles[idxParticle][2],lheParticles[idxParticle][3]); lCheck=true;}
	if(id==-11){ lbar.SetPxPyPzE(lheParticles[idxParticle][0],lheParticles[idxParticle][1],lheParticles[idxParticle][2],lheParticles[idxParticle][3]); lbarCheck=true;}
	if(id==12){ vl.SetPxPyPzE(lheParticles[idxParticle][0],lheParticles[idxParticle][1],lheParticles[idxParticle][2],lheParticles[idxParticle][3]); vlCheck=true;}
	if(id==-12){ vlbar.SetPxPyPzE(lheParticles[idxParticle][0],lheParticles[idxParticle][1],lheParticles[idxParticle][2],lheParticles[idxParticle][3]); vlbarCheck=true;}
	if(id==13){ l.SetPxPyPzE(lheParticles[idxParticle][0],lheParticles[idxParticle][1],lheParticles[idxParticle][2],lheParticles[idxParticle][3]); lCheck=true;}
	if(id==-13){ lbar.SetPxPyPzE(lheParticles[idxParticle][0],lheParticles[idxParticle][1],lheParticles[idxParticle][2],lheParticles[idxParticle][3]); lbarCheck=true;}
	if(id==14){ vl.SetPxPyPzE(lheParticles[idxParticle][0],lheParticles[idxParticle][1],lheParticles[idxParticle][2],lheParticles[idxParticle][3]); vlCheck=true;}
	if(id==-14){ vlbar.SetPxPyPzE(lheParticles[idxParticle][0],lheParticles[idxParticle][1],lheParticles[idxParticle][2],lheParticles[idxParticle][3]); vlbarCheck=true;}
	if(id==15){ l.SetPxPyPzE(lheParticles[idxParticle][0],lheParticles[idxParticle][1],lheParticles[idxParticle][2],lheParticles[idxParticle][3]); lCheck=true;}
	if(id==-15){ lbar.SetPxPyPzE(lheParticles[idxParticle][0],lheParticles[idxParticle][1],lheParticles[idxParticle][2],lheParticles[idxParticle][3]); lbarCheck=true;}
	if(id==16){ vl.SetPxPyPzE(lheParticles[idxParticle][0],lheParticles[idxParticle][1],lheParticles[idxParticle][2],lheParticles[idxParticle][3]); vlCheck=true;}
	if(id==-16){ vlbar.SetPxPyPzE(lheParticles[idxParticle][0],lheParticles[idxParticle][1],lheParticles[idxParticle][2],lheParticles[idxParticle][3]); vlbarCheck=true;}
	/*
	  if(id==11 || id==13 || id==15){
	  l.SetPxPyPzE(lheParticles[idxParticle][0], lheParticles[idxParticle][1], lheParticles[idxParticle][2], lheParticles[idxParticle][3]);
	  lCheck = true;
	  }
	  if(id==-11 || id==-13 || id==-15){
	  lbar.SetPxPyPzE(lheParticles[idxParticle][0], lheParticles[idxParticle][1], lheParticles[idxParticle][2], lheParticles[idxParticle][3]);
	  lbarCheck = true;
	  }
	  if(id==12 || id==14 || id==16){
	  vl.SetPxPyPzE(lheParticles[idxParticle][0], lheParticles[idxParticle][1], lheParticles[idxParticle][2], lheParticles[idxParticle][3]);
	  vlCheck = true;
	  }
	  if(id==-12 || id==-14 || id==-16){
	  vlbar.SetPxPyPzE(lheParticles[idxParticle][0], lheParticles[idxParticle][1], lheParticles[idxParticle][2], lheParticles[idxParticle][3]);
	  vlbarCheck = true;
	  }
	*/
      }
      if( lCheck && lbarCheck ) V_tlv = l + lbar; // ZtoLL
      if( vlCheck && vlbarCheck ) V_tlv = vl + vlbar; // ZtoNuNu
      if( lCheck && vlbarCheck ) V_tlv = l + vlbar; // WtoLNu
      if( lbarCheck && vlCheck ) V_tlv = lbar + vl; // WtoLNu
      lheV_pt = V_tlv.Pt(); // why not directly take from id==23 or 24?
    }
    EVENT_genHT = lheHt;
    EVENT_genPt = lheV_pt;
  } else {
    EVENT_originalXWGTUP_ = 1;
    EVENT_scalePDF_ = 1;
    EVENT_genWeights_.push_back(1);
    EVENT_Q2tthbbWeightUp_    = 1;
    EVENT_Q2tthbbWeightDown_  = 1;
    EVENT_PDFtthbbWeightUp_   = 1;
    EVENT_PDFtthbbWeightDown_ = 1;
  }
  edm::Handle<double> rhopogHandle;
  iEvent.getByToken(rhopogHandle_,rhopogHandle);
  double rhopog = *rhopogHandle;
  EVENT_rhopog_ = rhopog;

  edm::Handle<double> rhotthHandle;
  iEvent.getByToken(rhotthHandle_,rhotthHandle);
  double rhotth = *rhotthHandle;
  EVENT_rhotth_ = rhotth;

  edm::Handle<double> fixedGridRhoFastjetCentralHandle;
  iEvent.getByToken(fixedGridRhoFastjetCentralHandle_,fixedGridRhoFastjetCentralHandle);
  double fixedGridRhoFastjetCentral = *fixedGridRhoFastjetCentralHandle;
  EVENT_fixedGridRhoFastjetCentral = fixedGridRhoFastjetCentral;

  edm::Handle<double> fixedGridRhoFastjetCentralChargedPileUpHandle;
  iEvent.getByToken(fixedGridRhoFastjetCentralChargedPileUpHandle_,fixedGridRhoFastjetCentralChargedPileUpHandle);
  double fixedGridRhoFastjetCentralChargedPileUp = *fixedGridRhoFastjetCentralChargedPileUpHandle;
  EVENT_fixedGridRhoFastjetCentralChargedPileUp = fixedGridRhoFastjetCentralChargedPileUp;

  edm::Handle<double> fixedGridRhoFastjetCentralNeutralHandle;
  iEvent.getByToken(fixedGridRhoFastjetCentralNeutralHandle_,fixedGridRhoFastjetCentralNeutralHandle);
  double fixedGridRhoFastjetCentralNeutral = *fixedGridRhoFastjetCentralNeutralHandle;
  EVENT_fixedGridRhoFastjetCentralNeutral = fixedGridRhoFastjetCentralNeutral;

  //Event filters
  edm::Handle<edm::TriggerResults> metFilterBits;
  iEvent.getByToken(metFilterBits_, metFilterBits);
  const edm::TriggerNames &metNames = iEvent.triggerNames(*metFilterBits);
  for(unsigned int i = 0, n = metFilterBits->size(); i < n; ++i){ 
    if(metNames.triggerName(i)=="Flag_goodVertices")                       Flag_goodVertices                       = metFilterBits->accept(i);
    if(metNames.triggerName(i)=="Flag_globalSuperTightHalo2016Filter")     Flag_globalSuperTightHalo2016Filter     = metFilterBits->accept(i);
    if(metNames.triggerName(i)=="Flag_HBHENoiseFilter")                    Flag_HBHENoiseFilter                    = metFilterBits->accept(i);
    if(metNames.triggerName(i)=="Flag_HBHENoiseIsoFilter")                 Flag_HBHENoiseIsoFilter                 = metFilterBits->accept(i);
    if(metNames.triggerName(i)=="Flag_EcalDeadCellTriggerPrimitiveFilter") Flag_EcalDeadCellTriggerPrimitiveFilter = metFilterBits->accept(i);
    if(metNames.triggerName(i)=="Flag_BadPFMuonFilter")                    Flag_BadPFMuonFilter                    = metFilterBits->accept(i);
    if(metNames.triggerName(i)=="Flag_BadChargedCandidateFilter")          Flag_BadChargedCandidateFilter          = metFilterBits->accept(i);
    if(metNames.triggerName(i)=="Flag_eeBadScFilter")                      Flag_eeBadScFilter                      = metFilterBits->accept(i);
    if(metNames.triggerName(i)=="Flag_ecalBadCalibFilter")                 Flag_ecalBadCalibFilter                 = metFilterBits->accept(i);
  } //loop over met filters
  if(!_MC2016){
    edm::Handle< bool > passecalBadCalibFilterUpdate ;
    iEvent.getByToken(ecalBadCalibFilterUpdate_token,passecalBadCalibFilterUpdate);
    Flag_ecalBadCalibReducedMINIAODFilter =  (*passecalBadCalibFilterUpdate );
  }
}
void EventInfoSelector::SetBranches(){
  if(debug_) std::cout<<"setting branches: calling AddBranch of baseTree"<<std::endl;
  //Event quantities
  AddBranch(&EVENT_event_     ,"EVENT_event");
  AddBranch(&EVENT_run_       ,"EVENT_run");
  AddBranch(&EVENT_lumiBlock_ ,"EVENT_lumiBlock");
  AddBranch(&EVENT_genWeight_ ,"EVENT_genWeight");
  AddBranch(&EVENT_genWeights_,"EVENT_genWeights");
  AddBranch(&EVENT_genHT      ,"EVENT_genHT");
  AddBranch(&EVENT_genPt      ,"EVENT_genPt");
  AddBranch(&EVENT_prefiringweight_    ,"EVENT_prefiringweight");
  AddBranch(&EVENT_prefiringweightup_  ,"EVENT_prefiringweightup");
  AddBranch(&EVENT_prefiringweightdown_,"EVENT_prefiringweightdown");
  //Event filters
  AddBranch(&Flag_goodVertices                       ,"Flag_goodVertices");
  AddBranch(&Flag_globalSuperTightHalo2016Filter     ,"Flag_globalSuperTightHalo2016Filter");
  AddBranch(&Flag_HBHENoiseFilter                    ,"Flag_HBHENoiseFilter");
  AddBranch(&Flag_HBHENoiseIsoFilter                 ,"Flag_HBHENoiseIsoFilter");
  AddBranch(&Flag_EcalDeadCellTriggerPrimitiveFilter ,"Flag_EcalDeadCellTriggerPrimitiveFilter");
  AddBranch(&Flag_BadPFMuonFilter                    ,"Flag_BadPFMuonFilter");
  AddBranch(&Flag_BadChargedCandidateFilter          ,"Flag_BadChargedCandidateFilter");
  AddBranch(&Flag_eeBadScFilter                      ,"Flag_eeBadScFilter");
  AddBranch(&Flag_ecalBadCalibFilter                 ,"Flag_ecalBadCalibFilter");
  AddBranch(&Flag_ecalBadCalibReducedMINIAODFilter   ,"Flag_ecalBadCalibReducedMINIAODFilter");
  /*
    AddBranch(&EVENT_rhopog_    ,"EVENT_rhopog");
    AddBranch(&EVENT_rhotth_    ,"EVENT_rhotth");
    AddBranch(&EVENT_Q2tthbbWeightUp_    ,"EVENT_Q2tthbbWeightUp");
    AddBranch(&EVENT_Q2tthbbWeightDown_  ,"EVENT_Q2tthbbWeightDown");
    AddBranch(&EVENT_PDFtthbbWeightUp_   ,"EVENT_PDFtthbbWeightUp");
    AddBranch(&EVENT_PDFtthbbWeightDown_ ,"EVENT_PDFtthbbWeightDown");
    AddBranch(&EVENT_fixedGridRhoFastjetCentral               ,"EVENT_fixedGridRhoFastjetCentral");
    AddBranch(&EVENT_fixedGridRhoFastjetCentralChargedPileUp  ,"EVENT_fixedGridRhoFastjetCentralChargedPileUp");
    AddBranch(&EVENT_fixedGridRhoFastjetCentralNeutral        ,"EVENT_fixedGridRhoFastjetCentralNeutral");*/
}
void EventInfoSelector::Initialise(){
  //Event quantities
  EVENT_event_      = -9999;
  EVENT_run_        = -9999; 
  EVENT_lumiBlock_  = -9999;
  EVENT_genWeight_  = -9999;
  EVENT_genWeights_.clear();
  EVENT_genHT       = -9999;
  EVENT_genPt       = -9999;
  EVENT_prefiringweight_     = 1;
  EVENT_prefiringweightup_   = 1;
  EVENT_prefiringweightdown_ = 1;
  EVENT_rhopog_     = -9999;
  EVENT_rhotth_     = -9999; 
  EVENT_Q2tthbbWeightUp_    = -9999; 
  EVENT_Q2tthbbWeightDown_  = -9999; 
  EVENT_PDFtthbbWeightUp_   = -9999; 
  EVENT_PDFtthbbWeightDown_ = -9999; 
  EVENT_fixedGridRhoFastjetCentral              = -9999;
  EVENT_fixedGridRhoFastjetCentralChargedPileUp = -9999; 
  EVENT_fixedGridRhoFastjetCentralNeutral       = -9999;
  //Event filters
  Flag_goodVertices                       = -9999;
  Flag_globalSuperTightHalo2016Filter     = -9999;
  Flag_HBHENoiseFilter                    = -9999;
  Flag_HBHENoiseIsoFilter                 = -9999;
  Flag_EcalDeadCellTriggerPrimitiveFilter = -9999;
  Flag_BadPFMuonFilter                    = -9999;
  Flag_BadChargedCandidateFilter          = -9999;
  Flag_eeBadScFilter                      = -9999;
  Flag_ecalBadCalibFilter                 = -9999;
  Flag_ecalBadCalibReducedMINIAODFilter   = -9999;
}
