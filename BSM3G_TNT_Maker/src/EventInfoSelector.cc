#include "BSMFramework/BSM3G_TNT_Maker/interface/EventInfoSelector.h"
EventInfoSelector::EventInfoSelector(std::string name, TTree* tree, bool debug, const pset& iConfig, edm::ConsumesCollector && ic):
  baseTree(name,tree,debug)
{
  genEvtInfo_    = ic.consumes<GenEventInfoProduct>(edm::InputTag("generator"));
  lheEventProduct_ = ic.consumes<LHEEventProduct>(edm::InputTag("externalLHEProducer"));
  lheEventProductSource_ = ic.consumes<LHEEventProduct>(edm::InputTag("source"));// -> per genWeights di HCMN
  rhopogHandle_  = ic.consumes<double>(edm::InputTag("fixedGridRhoFastjetAll"));
  rhotthHandle_  = ic.consumes<double>(edm::InputTag("fixedGridRhoFastjetCentralNeutral"));
  fixedGridRhoFastjetCentralHandle_  = ic.consumes<double>(edm::InputTag("fixedGridRhoFastjetCentral"));
  fixedGridRhoFastjetCentralChargedPileUpHandle_  = ic.consumes<double>(edm::InputTag("fixedGridRhoFastjetCentralChargedPileUp"));
  fixedGridRhoFastjetCentralNeutralHandle_  = ic.consumes<double>(edm::InputTag("fixedGridRhoFastjetCentralNeutral"));
  metFilterBits_ = ic.consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("metFilterBits"));
  _is_data = iConfig.getParameter<bool>("is_data");
  if(debug) std::cout<<"in EventInfoSelector constructor"<<std::endl;
  SetBranches();
  read_PDFSet = new (LHAPDF::PDFSet)("NNPDF30_nlo_as_0118");
  _systPDFs = read_PDFSet->mkPDFs();
  bjetnessproducer_ = iConfig.getParameter<bool>("bjetnessproducer");
  bjetness_     = ic.consumes<vector<float> >(edm::InputTag("BJetness", "BJetnessValue"));
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
  edm::Handle<LHEEventProduct> lheEventProductSource;
  iEvent.getByToken(lheEventProductSource_, lheEventProductSource);
  if(!_is_data){
    EVENT_genWeight_ = genEvtInfo->weight();
    const GenEventInfoProduct& genEventInfoW = *(genEvtInfo.product());
    const gen::PdfInfo* pdf = genEventInfoW.pdf();
    EVENT_scalePDF_ = pdf->scalePDF;
    if(lheEventProduct.isValid()){
      EVENT_originalXWGTUP_ = lheEventProduct->originalXWGTUP();
      for (unsigned int i=0; i<lheEventProduct->weights().size(); i++){
        EVENT_genWeights_.push_back(lheEventProduct->weights()[i].wgt);
        //Q2 for ttHbb synchronization
        if(lheEventProduct->weights()[i].id == "1005") EVENT_Q2tthbbWeightUp_   = lheEventProduct->weights()[i].wgt/lheEventProduct->originalXWGTUP(); 
        if(lheEventProduct->weights()[i].id == "1009") EVENT_Q2tthbbWeightDown_ = lheEventProduct->weights()[i].wgt/lheEventProduct->originalXWGTUP(); 
      }
    }
    if(lheEventProductSource.isValid()){
      for(unsigned int i=0; i<lheEventProductSource->weights().size(); i++) EVENT_genWeightsSource_.push_back(lheEventProductSource->weights()[i].wgt);
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
    if(lheEventProduct.isValid()){
      const lhef::HEPEUP& lheEvent = lheEventProduct->hepeup();
      std::vector<lhef::HEPEUP::FiveVector> lheParticles = lheEvent.PUP;
      size_t numParticles = lheParticles.size();
      for(size_t idxParticle = 0; idxParticle < numParticles; ++idxParticle){
        int absPdgId = TMath::Abs(lheEvent.IDUP[idxParticle]);
        int status = lheEvent.ISTUP[idxParticle];
        if(status == 1 && ((absPdgId >= 1 && absPdgId <= 6) || absPdgId == 21)){ // quarks and gluons
          lheHt += TMath::Sqrt(TMath::Power(lheParticles[idxParticle][0], 2.) + TMath::Power(lheParticles[idxParticle][1], 2.)); // first entry is px, second py
        } 
      }
    }
    EVENT_genHT = lheHt;
  } else {
    EVENT_originalXWGTUP_ = 1;
    EVENT_scalePDF_ = 1;
    EVENT_genWeights_.push_back(1);
    EVENT_genWeightsSource_.push_back(1);
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
    if(metNames.triggerName(i)=="Flag_HBHENoiseFilter")                    Flag_HBHENoiseFilter                    = metFilterBits->accept(i);
    if(metNames.triggerName(i)=="Flag_HBHENoiseIsoFilter")                 Flag_HBHENoiseIsoFilter                 = metFilterBits->accept(i);
    if(metNames.triggerName(i)=="Flag_CSCTightHaloFilter")                 Flag_CSCTightHaloFilter                 = metFilterBits->accept(i);
    if(metNames.triggerName(i)=="Flag_CSCTightHaloTrkMuUnvetoFilter")      Flag_CSCTightHaloTrkMuUnvetoFilter      = metFilterBits->accept(i);
    if(metNames.triggerName(i)=="Flag_CSCTightHalo2015Filter")             Flag_CSCTightHalo2015Filter             = metFilterBits->accept(i);
    if(metNames.triggerName(i)=="Flag_HcalStripHaloFilter")                Flag_HcalStripHaloFilter                = metFilterBits->accept(i);
    if(metNames.triggerName(i)=="Flag_hcalLaserEventFilter")               Flag_hcalLaserEventFilter               = metFilterBits->accept(i);
    if(metNames.triggerName(i)=="Flag_EcalDeadCellTriggerPrimitiveFilter") Flag_EcalDeadCellTriggerPrimitiveFilter = metFilterBits->accept(i);
    if(metNames.triggerName(i)=="Flag_EcalDeadCellBoundaryEnergyFilter")   Flag_EcalDeadCellBoundaryEnergyFilter   = metFilterBits->accept(i);
    if(metNames.triggerName(i)=="Flag_goodVertices")                       Flag_goodVertices                       = metFilterBits->accept(i);
    if(metNames.triggerName(i)=="Flag_eeBadScFilter")                      Flag_eeBadScFilter                      = metFilterBits->accept(i);
    if(metNames.triggerName(i)=="Flag_ecalLaserCorrFilter")                Flag_ecalLaserCorrFilter                = metFilterBits->accept(i);
    if(metNames.triggerName(i)=="Flag_trkPOGFilters")                      Flag_trkPOGFilters                      = metFilterBits->accept(i);
    if(metNames.triggerName(i)=="Flag_chargedHadronTrackResolutionFilter") Flag_chargedHadronTrackResolutionFilter = metFilterBits->accept(i);
    if(metNames.triggerName(i)=="Flag_muonBadTrackFilter")                 Flag_muonBadTrackFilter                 = metFilterBits->accept(i);
    if(metNames.triggerName(i)=="Flag_trkPOG_manystripclus53X")            Flag_trkPOG_manystripclus53X            = metFilterBits->accept(i);
    if(metNames.triggerName(i)=="Flag_trkPOG_toomanystripclus53X")         Flag_trkPOG_toomanystripclus53X         = metFilterBits->accept(i);
    if(metNames.triggerName(i)=="Flag_trkPOG_logErrorTooManyClusters")     Flag_trkPOG_logErrorTooManyClusters     = metFilterBits->accept(i);
    if(metNames.triggerName(i)=="Flag_METFilters")                         Flag_METFilters                         = metFilterBits->accept(i);
  } //loop over met filters
  //BJetness variables
  if(bjetnessproducer_){
    Handle<vector<float> > bjetness;
    iEvent.getByToken(bjetness_, bjetness);
    BJetness_jetpt0          = (*bjetness)[0]; 
    BJetness_numleps         = (*bjetness)[1]; 
    BJetness_npvTrkOVcollTrk = (*bjetness)[2];
    BJetness_avip3dval       = (*bjetness)[3];
    BJetness_avip3dsig       = (*bjetness)[4];
    BJetness_avsip3dsig      = (*bjetness)[5];
    BJetness_avip1dsig       = (*bjetness)[6];
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
  AddBranch(&EVENT_genWeightsSource_,"EVENT_genWeightsSource");
  AddBranch(&EVENT_genHT      ,"EVENT_genHT");
  AddBranch(&EVENT_rhopog_    ,"EVENT_rhopog");
  AddBranch(&EVENT_rhotth_    ,"EVENT_rhotth");
  AddBranch(&EVENT_Q2tthbbWeightUp_    ,"EVENT_Q2tthbbWeightUp");
  AddBranch(&EVENT_Q2tthbbWeightDown_  ,"EVENT_Q2tthbbWeightDown");
  AddBranch(&EVENT_PDFtthbbWeightUp_   ,"EVENT_PDFtthbbWeightUp");
  AddBranch(&EVENT_PDFtthbbWeightDown_ ,"EVENT_PDFtthbbWeightDown");
  AddBranch(&EVENT_fixedGridRhoFastjetCentral               ,"EVENT_fixedGridRhoFastjetCentral");
  AddBranch(&EVENT_fixedGridRhoFastjetCentralChargedPileUp  ,"EVENT_fixedGridRhoFastjetCentralChargedPileUp");
  AddBranch(&EVENT_fixedGridRhoFastjetCentralNeutral        ,"EVENT_fixedGridRhoFastjetCentralNeutral");
  //Event filters
  AddBranch(&Flag_HBHENoiseFilter                    ,"Flag_HBHENoiseFilter");
  AddBranch(&Flag_HBHENoiseIsoFilter                 ,"Flag_HBHENoiseIsoFilter");
  AddBranch(&Flag_CSCTightHaloFilter                 ,"Flag_CSCTightHaloFilter");
  AddBranch(&Flag_CSCTightHaloTrkMuUnvetoFilter      ,"Flag_CSCTightHaloTrkMuUnvetoFilter");
  AddBranch(&Flag_CSCTightHalo2015Filter             ,"Flag_CSCTightHalo2015Filter");
  AddBranch(&Flag_HcalStripHaloFilter                ,"Flag_HcalStripHaloFilter");
  AddBranch(&Flag_hcalLaserEventFilter               ,"Flag_hcalLaserEventFilter");
  AddBranch(&Flag_EcalDeadCellTriggerPrimitiveFilter ,"Flag_EcalDeadCellTriggerPrimitiveFilter");
  AddBranch(&Flag_EcalDeadCellBoundaryEnergyFilter   ,"Flag_EcalDeadCellBoundaryEnergyFilter");
  AddBranch(&Flag_goodVertices                       ,"Flag_goodVertices");
  AddBranch(&Flag_eeBadScFilter                      ,"Flag_eeBadScFilter");
  AddBranch(&Flag_ecalLaserCorrFilter                ,"Flag_ecalLaserCorrFilter");
  AddBranch(&Flag_trkPOGFilters                      ,"Flag_trkPOGFilters");
  AddBranch(&Flag_chargedHadronTrackResolutionFilter ,"Flag_chargedHadronTrackResolutionFilter");
  AddBranch(&Flag_muonBadTrackFilter                 ,"Flag_muonBadTrackFilter");
  AddBranch(&Flag_trkPOG_manystripclus53X            ,"Flag_trkPOG_manystripclus53X");
  AddBranch(&Flag_trkPOG_toomanystripclus53X         ,"Flag_trkPOG_toomanystripclus53X");
  AddBranch(&Flag_trkPOG_logErrorTooManyClusters     ,"Flag_trkPOG_logErrorTooManyClusters");
  AddBranch(&Flag_METFilters                         ,"Flag_METFilters");
  //BJetness variables
  AddBranch(&BJetness_jetpt0               ,"BJetness_jetpt0");
  AddBranch(&BJetness_numleps              ,"BJetness_numleps");
  AddBranch(&BJetness_npvTrkOVcollTrk      ,"BJetness_npvTrkOVcollTrk");
  AddBranch(&BJetness_avip3dval            ,"BJetness_avip3dval");
  AddBranch(&BJetness_avip3dsig            ,"BJetness_avip3dsig");
  AddBranch(&BJetness_avsip3dsig           ,"BJetness_avsip3dsig");
  AddBranch(&BJetness_avip1dsig            ,"BJetness_avip1dsig");
}
void EventInfoSelector::Initialise(){
  //Event quantities
  EVENT_event_      = -9999;
  EVENT_run_        = -9999; 
  EVENT_lumiBlock_  = -9999;
  EVENT_genWeight_  = -9999;
  EVENT_genWeights_.clear();
  EVENT_genWeightsSource_.clear();
  EVENT_genHT       = -9999;
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
  Flag_HBHENoiseFilter                    = -9999;
  Flag_HBHENoiseIsoFilter                 = -9999;
  Flag_CSCTightHaloFilter                 = -9999;
  Flag_CSCTightHaloTrkMuUnvetoFilter      = -9999;
  Flag_CSCTightHalo2015Filter             = -9999;
  Flag_HcalStripHaloFilter                = -9999;
  Flag_hcalLaserEventFilter               = -9999;
  Flag_EcalDeadCellTriggerPrimitiveFilter = -9999;
  Flag_EcalDeadCellBoundaryEnergyFilter   = -9999;
  Flag_goodVertices                       = -9999;
  Flag_eeBadScFilter                      = -9999;
  Flag_ecalLaserCorrFilter                = -9999;
  Flag_trkPOGFilters                      = -9999;
  Flag_chargedHadronTrackResolutionFilter = -9999;
  Flag_muonBadTrackFilter                 = -9999;
  Flag_trkPOG_manystripclus53X            = -9999;
  Flag_trkPOG_toomanystripclus53X         = -9999;
  Flag_trkPOG_logErrorTooManyClusters     = -9999;
  Flag_METFilters                         = -9999;
  //BJetness variables
  BJetness_jetpt0          = -9999;
  BJetness_numleps         = -9999;
  BJetness_npvTrkOVcollTrk = -9999;
  BJetness_avip3dval       = -9999;
  BJetness_avip3dsig       = -9999;
  BJetness_avsip3dsig      = -9999;
  BJetness_avip1dsig       = -9999;  
}
