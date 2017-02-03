#include "BSMFramework/BSM3G_TNT_Maker/interface/PileupReweight.h"
PileupReweight::PileupReweight(std::string name, TTree* tree, bool debug, const pset& iConfig, edm::ConsumesCollector && ic):
  baseTree(name,tree,debug)
{
  if(debug) std::cout<<"in PileupReweight constructor"<<std::endl;
  _MiniAODv2 = iConfig.getParameter<bool>("MiniAODv2");
  _is_data   = iConfig.getParameter<bool>("is_data");
  PUInfo_         = ic.consumes<std::vector< PileupSummaryInfo> >(edm::InputTag("slimmedAddPileupInfo"));
  PUReweightfile_ = iConfig.getParameter<edm::FileInPath>("PUReweightfile");
  // Get data distribution from file
  const char *filePath = PUReweightfile_.fullPath().c_str();
  TFile file(filePath, "READ");
  TH1* h = NULL;
  file.GetObject("pileup",h);
  if( h == NULL ) {
    std::cerr << "\n\nERROR in PUWeight: Histogram 'pileup' does not exist in file \n.";
    throw std::exception();
  }
  h->SetDirectory(0);
  file.Close();

  // Computing weights
  // Store probabilites for each pu bin
  unsigned int nPUMax = 0;
  double *npuProbs = 0;
  
  //------------ 2015 PU SCENARIO ------------//
  //nPUMax = 52;
  //double npuWinter15_25ns[52] = {4.8551E-07,1.74806E-06,3.30868E-06,1.62972E-05,4.95667E-05,0.000606966,0.003307249,0.010340741,0.022852296,0.041948781,0.058609363,0.067475755,0.072817826,0.075931405,0.076782504,0.076202319,0.074502547,0.072355135,0.069642102,0.064920999,0.05725576,0.047289348,0.036528446,0.026376131,0.017806872,0.011249422,0.006643385,0.003662904,0.001899681,0.00095614,0.00050028,0.000297353,0.000208717,0.000165856,0.000139974,0.000120481,0.000103826,8.88868E-05,7.53323E-05,6.30863E-05,5.21356E-05,4.24754E-05,3.40876E-05,2.69282E-05,2.09267E-05,1.5989E-05,4.8551E-06,2.42755E-06,4.8551E-07,2.42755E-07,1.21378E-07,4.8551E-08};
  //npuProbs = npuWinter15_25ns;
  //------------------------------------------//

  //------------ 2016 PU SCENARIO ------------//
  //nPUMax = 50;
  //double npu_mix_2016_25ns_SpringMC_PUScenarioV1[nPUMax] = {0.000829312873542,0.00124276120498,0.00339329181587,0.00408224735376,0.00383036590008,0.00659159288946,0.00816022734493,0.00943640833116,0.0137777376066,0.017059392038,0.0213193035468,0.0247343174676,0.0280848773878,0.0323308476564,0.0370394341409,0.0456917721191,0.0558762890594,0.0576956187107,0.0625325287017,0.0591603758776,0.0656650815128,0.0678329011676,0.0625142146389,0.0548068448797,0.0503893295063,0.040209818868,0.0374446988111,0.0299661572042,0.0272024759921,0.0219328403791,0.0179586571619,0.0142926728247,0.00839941654725,0.00522366397213,0.00224457976761,0.000779274977993,0.000197066585944,7.16031761328e-05,1.0e-08,1.0e-08,1.0e-08,1.0e-08,1.0e-08,1.0e-08,1.0e-08,1.0e-08,1.0e-08,1.0e-08,1.0e-08,1.0e-08};//0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  //npuProbs = npu_mix_2016_25ns_SpringMC_PUScenarioV1;
  //------------------------------------------//

  //------------ Moriond17 PU SCENARIO ------------//
  nPUMax =75;
  double npu_Moriond17Scenario[nPUMax] = {1.78653e-05 ,2.56602e-05 ,5.27857e-05 ,8.88954e-05 ,0.000109362 ,0.000140973 ,0.000240998 ,0.00071209 ,0.00130121 ,0.00245255 ,0.00502589 ,0.00919534 ,0.0146697 ,0.0204126 ,0.0267586 ,0.0337697 ,0.0401478 ,0.0450159 ,0.0490577 ,0.0524855 ,0.0548159 ,0.0559937 ,0.0554468 ,0.0537687 ,0.0512055 ,0.0476713 ,0.0435312 ,0.0393107 ,0.0349812 ,0.0307413 ,0.0272425 ,0.0237115 ,0.0208329 ,0.0182459 ,0.0160712 ,0.0142498 ,0.012804 ,0.011571 ,0.010547 ,0.00959489 ,0.00891718 ,0.00829292 ,0.0076195 ,0.0069806 ,0.0062025 ,0.00546581 ,0.00484127 ,0.00407168 ,0.00337681 ,0.00269893 ,0.00212473 ,0.00160208 ,0.00117884 ,0.000859662 ,0.000569085 ,0.000365431 ,0.000243565 ,0.00015688 ,9.88128e-05 ,6.53783e-05 ,3.73924e-05 ,2.61382e-05 ,2.0307e-05 ,1.73032e-05 ,1.435e-05 ,1.36486e-05 ,1.35555e-05 ,1.37491e-05 ,1.34255e-05 ,1.33987e-05 ,1.34061e-05 ,1.34211e-05 ,1.34177e-05 ,1.32959e-05 ,1.33287e-05};
  npuProbs = npu_Moriond17Scenario;
  //------------------------------------------//

  // Check that binning of data-profile matches MC scenario
  if( nPUMax != static_cast<unsigned int>(h->GetNbinsX()) ) {
    std::cerr << "\n\nERROR number of bins (" << h->GetNbinsX() << ") in data PU-profile does not match number of bins (" << nPUMax << ") in MC scenario " << std::endl;
    throw std::exception();
  }

  std::vector<double> result(nPUMax,0.);
  double s = 0.;
  for(unsigned int npu = 0; npu < nPUMax; ++npu) {
    const double npu_estimated = h->GetBinContent(h->GetXaxis()->FindBin(npu));
    result[npu] = npu_estimated / npuProbs[npu];
    s += npu_estimated;
  }
  // normalize weights such that the total sum of weights over thw whole sample is 1.0, i.e., sum_i  result[i] * npu_probs[i] should be 1.0 (!)
  for(unsigned int npu = 0; npu < nPUMax; ++npu) {
    result[npu] /= s;
  }

  puWeigths_ = result;
  nPUMax_ = puWeigths_.size();

  // Clean up
  delete h;

  SetBranches();
}
PileupReweight::~PileupReweight(){
  delete tree_;
}

void PileupReweight::Fill(const edm::Event& iEvent){
  if(debug_) std::cout<<"getting PileupReweight info"<<std::endl;
  double w = 1.;
  if(!_is_data) {
    Handle<std::vector< PileupSummaryInfo > >  PUInfo;
    iEvent.getByToken(PUInfo_, PUInfo);
    std::vector<PileupSummaryInfo>::const_iterator PVI;
    float nPU = -1;
    for(PVI = PUInfo->begin(); PVI != PUInfo->end(); ++PVI) {
      int BX = PVI->getBunchCrossing();
      if(BX == 0) { 
	nPU = PVI->getTrueNumInteractions();
	continue;
      }
    }
    if( nPU >= nPUMax_ ) {
      //std::cerr << "WARNING: Number of PU vertices = " << nPU << " out of histogram binning." << std::endl;
      // In case nPU is out-of data-profile binning,
      // use weight from last bin
      w = puWeigths_.back();
    } else {
      w = puWeigths_.at(nPU);
    }
  }
  PUWeight=w;
  if(debug_) std::cout<<"got PileupReweight info"<<std::endl;
}
void PileupReweight::SetBranches(){
  if(debug_) std::cout<<"setting branches: calling AddBranch of PileupReweight"<<std::endl;
  AddBranch(&PUWeight,"PUWeight");
}
