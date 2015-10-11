/**
This Macro   
1. Plots variables in stack mode 

Need to specify
1. See Declare Constants
2. Pay attention if you have to use the proper path/name of the tree (TNT/BOOM or BOOM or whatelse)
3. Remember that if doasym is true -> v = n-1: where v is in bin[v] and n is in asymbin[n]
   Adjust also the range
*/
/////
//   To run: root -l StackPlots.cc+  
/////
/////
//   Prepare Root and Roofit
/////
#include "TFile.h"
#include "TH1F.h"
#include "THStack.h"
#include "TLegend.h"
#include "TTree.h"
#include "TTreePlayer.h"
#include "TStyle.h"
#include "TGaxis.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TEfficiency.h"
#include "TBranch.h"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <iostream>
using namespace std;
/////
//   Declare constants
/////
//Path and root files
const string path           = "/afs/cern.ch/work/f/fromeo/CMSSW_7_4_7/src/BSMFramework/BSM3G_TNT_Maker/macros/";
const char *samples[]       = {"WW", "WZ"};
const string selection      = "_10K";//Must write as _NameSelection
const string specsel        = "";
const string channel        = "";
//Corrections
const double xSecs[]        = {831.76, 4895., 61526.7, 63.21, 22.82, 10.32, 35.6, 35.6};//pb Signal:0.0007
const double evtsRead[]     = {10000., 10000.,10000.,10000., 10000.,10000.,10000., 10000.};//, 995600., 1000000., 7361970., 989608., 996920.,998848.};
const double Luminosity     = 1000.; //pb^-1
const bool LumiNorm         = true;
const bool PUcorr           = false;
const bool GenWgtcorr       = false;
//Plots
const bool show_ratio       = false; //true only with data
const bool save_plots       = false;
//Variables
const int    numvar      = 2;
const int    col_size    = 500; //500 = highest bin
const bool   doasym      = false; //If true, set bin[numvar] for the asymbinvar = num bin you choose in asym configuration
const double asymbin[11] = {0,50,100,150,200,250,300,400,600,900,1500};
//Impact parameter
const int inivar  = 0;
const int finvar  = 1;
const int posvtcr = 0;//leading 0, sublead 1, etc
//Kinematics
const char *variables[]         = {"Muon_pt","Muon_eta"};
const char *vartitleXaxis[]     = {"Mu pT (GeV)","Mu eta"};
//1=int,2=double,3=vec_int,4=vec_double //2,3 has to be implemented yet
const int    vartype[numvar]    = {   4,   4 };
const int    bin[numvar]        = { 100,  25 };
const double inRange[numvar]    = {   0,-2.5 };
const double endRange[numvar]   = {1000, 2.5 };
const int    ylogscale[numvar]  = {   0,   0 };
const double setminimum[numvar] = {   0,   0 };
const double setmaximum[numvar] = {  50,  10 };
const bool   show_title         = true;
/////
//   Declare functions 
/////
TFile* Call_TFile(string rootpla);
TH1F* get_emptyhisto(string rootpla, int v);
TH1F* get_h_var(unsigned int v, string var, string vaT, uint r, string rootplas, double err_AllBkg[][col_size], double ent_AllBkg[][col_size]);
void draw_plot(TH1F* h_data_var, THStack *hstack, TH1F* h_sum_var, TLegend *leg, int v, string var, int rootplas_size, double err_AllBkg[][col_size]);
void draw_plot_withratio(TCanvas* c1, TH1F* h_data_var, THStack *hstack, TH1F* h_sum_var, TLegend *leg, int v, string var, string vartitlexaxis, int rootplas_size, double err_AllBkg[][col_size]);
void setTDRStyle();
/////
//   Main function
/////
void StackPlots_2(){
 setTDRStyle();
 //Loop over all variables
 vector<string> variablesv(variables, variables + sizeof(variables)/sizeof(variables[0]));
 vector<string> vartitleXaxisv(vartitleXaxis, vartitleXaxis + sizeof(vartitleXaxis)/sizeof(vartitleXaxis[0]));
 for(int v=inivar; v<finvar; v++){
  //Legend  
  TLegend *leg = new TLegend(0.75, 0.25, 1.0, 0.9);
  leg->SetHeader("Samples");
  leg->SetBorderSize(0);
  leg->SetTextSize(0.05);
  //MC
  THStack *hstack = new THStack("hstack","hstack");
  double bkgstackintegral = 0.;
  TH1F *h_sum_var;
  if(variablesv[v]=="asymbinvar" && doasym) h_sum_var = new TH1F("","",bin[v],asymbin);
  else                                      h_sum_var = new TH1F("","",bin[v],inRange[v],endRange[v]);
  h_sum_var->Sumw2();
  //Data
  TH1F* h_data_var = get_emptyhisto(variablesv[v],v);
  //Loop over samples
  vector<string> rootplas(samples, samples + sizeof(samples)/sizeof(samples[0]));
  const uint rootplas_size = rootplas.size();
  double err_AllBkg[rootplas_size][col_size];
  double ent_AllBkg[rootplas_size][col_size];
  for(uint r=0; r<rootplas_size; r++) for(int b=0; b<bin[v]; b++) err_AllBkg[r][b] = 0.;
  for(uint r=0; r<rootplas_size; r++) for(int b=0; b<bin[v]; b++) ent_AllBkg[r][b] = 0.;
  for(uint r=0; r<rootplas_size; r++){
   //Declare histograms for variables
   TH1F *h_var;
   if(variablesv[v]=="asymbinvar" && doasym)  h_var = new TH1F(variablesv[v].c_str(),variablesv[v].c_str(),bin[v],asymbin);
   else                                       h_var = new TH1F(variablesv[v].c_str(),variablesv[v].c_str(),bin[v],inRange[v],endRange[v]);
   //Get histo
   if(rootplas[r]!="data") h_var      = get_h_var(v,variablesv[v],vartitleXaxisv[v],r,rootplas[r],err_AllBkg,ent_AllBkg);
   else                    h_data_var = get_h_var(v,variablesv[v],vartitleXaxisv[v],r,rootplas[r],err_AllBkg,ent_AllBkg);
   if(rootplas[r]!="data"){ 
    //cout << "entro in questo if?" << endl;
    //Put histos in the hstack 
    int col = -1;
    if(rootplas[r].substr(0,11)=="ADD2400")       col = 0;//Leave as an example for signal samples
    if(rootplas[r]=="TT")                         col = 1;
    if(rootplas[r]=="WJets")                      col = 2;
    if(rootplas[r]=="QCD")                        col = 3;
    if(rootplas[r]=="DY")                         col = 4;
    if(rootplas[r]=="WW")                         col = 5;
    if(rootplas[r]=="WZ")                         col = 6;
    if(rootplas[r]=="ZZ")                         col = 7;
    if(rootplas[r]=="tW" || rootplas[r]=="atW") col = 28;
    h_var->SetFillColor(col+3);
    h_var->SetLineColor(col+3);
    hstack->Add(h_var);
    leg->AddEntry(h_var,rootplas[r].c_str(),"L");
    //Sum them for the error
    h_sum_var->Add(h_sum_var,h_var); 
    cout<<"Evt "<<rootplas[r]<<" "<<h_var->Integral()<<endl;
    
    //cout<<variablesv[v]<<endl;

     //h_var->Draw();
     //hstack->Draw();


   }else{
    cout<<"Evt "<<rootplas[r]<<" "<<h_data_var->Integral()<<endl;
   }
   if(rootplas[r].substr(0,11)!="SSMToTauTau") bkgstackintegral += h_var->Integral();//Leave as an example for signal samples
  }//End rootplas
  cout<<"Tot "<<bkgstackintegral<<endl;  
  hstack->Draw("");
  hstack->Draw("same");
  /*
  //Find highest bin value
  double higherBinVal = 0;
  for(int h=1; h<=h_data_var->GetNbinsX(); h++) if(h_data_var->GetBinContent(h)>higherBinVal) higherBinVal=h_data_var->GetBinContent(h);
  if(hstack->GetMaximum()>higherBinVal)                                                       higherBinVal=hstack->GetMaximum();
  if(ylogscale[v]==1) h_data_var->SetMaximum((higherBinVal+0.25*higherBinVal)*10);
  else                h_data_var->SetMaximum(higherBinVal+0.25*higherBinVal);
  if(ylogscale[v]==1) h_data_var->SetMinimum(0.1);
  else                h_data_var->SetMinimum(0);
  if(ylogscale[v]==1) hstack->SetMaximum((higherBinVal+0.25*higherBinVal)*10);
  else                hstack->SetMaximum(higherBinVal+0.25*higherBinVal);
  if(ylogscale[v]==1) hstack->SetMinimum(0.1);
  else                hstack->SetMinimum(0);
  /////
  //   Draw
  /////
  //Canvas
  TCanvas* c1 = new TCanvas(variablesv[v].c_str(),variablesv[v].c_str(),200,200,700,600);
  if(ylogscale[v]==1) c1->SetLogy();

  //if(show_ratio)  draw_plot_withratio(c1, h_data_var,hstack,h_sum_var,leg,v,variablesv[v],vartitleXaxisv[v],rootplas_size,err_AllBkg);
  //if(!show_ratio) draw_plot(h_data_var,hstack,h_sum_var,leg,v,variablesv[v],rootplas_size,err_AllBkg);
  //Line
  TLine* line1 = new TLine(-0.75,0,-0.75,higherBinVal+0.25*higherBinVal);
  line1->SetLineColor(kRed);
  line1->SetLineWidth(3);
  //line1->Draw("same");
  //Save plots
  c1->Update();
  string namefile = specsel+variablesv[v]+selection+".png";;
  if(save_plots) c1->SaveAs(namefile.c_str());
  */
 }
}
/////
//   Call TFile to be read
/////
TFile* Call_TFile(string rootpla){
 string file_name = path+specsel+rootpla+selection+".root";
 //cout<<"file_name "<<file_name<<endl;
 TFile* f = new TFile(file_name.c_str(),"update");
 return f;
}
/////
//   Get an empty histogram
/////
TH1F* get_emptyhisto(string name, int v){
 TH1F *hist;
 if(name=="asymbinvar" && doasym) hist = new TH1F(name.c_str(),name.c_str(),bin[v],asymbin);
 else                                      hist = new TH1F(name.c_str(),name.c_str(),bin[v],inRange[v],endRange[v]); 
 hist->SetTitle(0);
 hist->SetMarkerStyle(8);
 hist->SetMarkerColor(1);
 hist->SetLineColor(1);
 vector<string> vartitleXaxisv(vartitleXaxis, vartitleXaxis + sizeof(vartitleXaxis)/sizeof(vartitleXaxis[0]));
 if(!show_ratio) hist->GetXaxis()->SetTitle(vartitleXaxisv[v].c_str());
 char bin_size_c[col_size]; float bin_size_f = ((endRange[v]-inRange[v])/bin[v]); sprintf(bin_size_c,"%.2f",bin_size_f); 
 string titleYaxis;
 //if(variablesv[v]=="asymbinvar" && doasym){
    titleYaxis = "Events";
 //}else{
 // titleYaxis = "Events/"+(string) bin_size_c;
 //}
 hist->GetYaxis()->SetTitle(titleYaxis.c_str());
return hist;
}
/////
//   Fill histo with double type
/////
TH1F* get_h_var(unsigned int v, string var, string varT, uint r, string rootplas, double err_AllBkg[][col_size], double ent_AllBkg[][col_size]){
 //Take histo
 TH1F *hist;
 if(var=="asymbinvar" && doasym) hist = new TH1F("hist","hist",bin[v],asymbin);
 else                            hist = new TH1F("hist","hist",bin[v],inRange[v],endRange[v]);
 TH1F *hist_err;
 if(var=="asymbinvar" && doasym) hist_err = new TH1F("hist_err","hist_err",bin[v],asymbin);
 else                            hist_err = new TH1F("hist_err","hist_err",bin[v],inRange[v],endRange[v]);
 hist_err->Sumw2();
 //Call tree and variables
 TFile* f = Call_TFile(rootplas.c_str()); TTree* tree; f->GetObject("TNT/BOOM",tree);
 Int_t curr_var_i;
 TBranch *b_curr_var_i;
 vector <double> * curr_var_vd;
 TBranch *b_curr_var_vd;
 if(vartype[v]==1){
  curr_var_i = 0;
  tree->SetBranchAddress(var.c_str(), &curr_var_i, &b_curr_var_i);
 }else if(vartype[v]==4){
  curr_var_vd = 0;

  //cout<<b_curr_var_vd->at(0)<<endl;

  tree->SetBranchAddress(var.c_str(), &curr_var_vd, &b_curr_var_vd);
 }
 //Run over entries
 vector<double> xsecs(xSecs, xSecs + sizeof(xSecs)/sizeof(xSecs[0]));
 vector<double> evtsread(evtsRead, evtsRead + sizeof(evtsRead)/sizeof(evtsRead[0]));
 //for(int en=0; en<tree->GetEntries(); en++)
 for(int en=0; en<10; en++)
 {
  Long64_t tentry = tree->LoadTree(en);
  //Get weights 
  double evt_wgt = 1.;
  if(rootplas!="data"){
   //Luminosity normalisation
   //(tough this is a global sample property, you may want to keep it here to avoid mistakes when changin the samples)
   if(LumiNorm){
    double lumi_wgt = xsecs[r]/evtsread[r]*Luminosity;
    evt_wgt = evt_wgt*lumi_wgt;
    cout << "evt_wgt" << evt_wgt << endl;
   }
   if(PUcorr){//PU correction (not included yet)
    double evt_pu = 1.;
    evt_wgt = evt_wgt*evt_pu;
   }
   if(GenWgtcorr){//GenWeight correction for aMC@NLO samples (not included yet)
    double evt_genwg = 1.;
    evt_wgt = evt_wgt*evt_genwg;
   }
  }
  //Get values
  if(vartype[v]==1){
   b_curr_var_i->GetEntry(tentry);
   hist->Fill(curr_var_i,evt_wgt);
   hist_err->Fill(curr_var_i,evt_wgt*evt_wgt);
  }else if(vartype[v]==4){
   
//    double lumi_wgt = xsecs[r]/evtsread[r]*Luminosity
//    //evt_wgt*lumi_wgt;
  //  evt_wgt = 1;

 
   b_curr_var_vd->GetEntry(tentry);
   if(curr_var_vd->size()>0){

    hist->Fill(curr_var_vd->at(posvtcr),evt_wgt);
    hist_err->Fill(curr_var_vd->at(posvtcr),evt_wgt*evt_wgt);

//cout<<"quale pt considero? "<<posvtcr<<endl;

   }else{
 //cout<<"Variable size of "<<var<<" is 0"<<endl;
   }
  }
 }
 if(rootplas!="data") for(int b=0; b<bin[v]; b++){
  err_AllBkg[r][b] = sqrt(hist_err->GetBinContent(b+1));
  ent_AllBkg[r][b] = hist->GetBinContent(b+1);
 }
 delete tree;
 return hist;
}
/////
//   Draw plot
/////
void draw_plot(TH1F* h_data_var, THStack *hstack, TH1F* h_sum_var, TLegend *leg, int v, string var, int rootplas_size, double err_AllBkg[][col_size]){
 //Data and bkg
 if(h_data_var->GetEntries()==0) gStyle->SetOptStat(0);
 TGaxis::SetMaxDigits(4);
 h_data_var->Draw("P");
 //hstack->Draw("textsame");
 hstack->Draw("same");//same
/// h_data_var->Draw("PEsame");
 if(show_title) h_data_var->SetTitle("#scale[0.90]{     Work in progress preliminary,   #sqrt{s} = 13 TeV, L = 47.9 pb^{-1}}");
 h_data_var->GetXaxis()->SetRangeUser(inRange[v],endRange[v]);
 h_data_var->GetXaxis()->SetLimits(inRange[v],endRange[v]);
 //h_data_var->Draw("text45same");
 gPad->RedrawAxis();
 if(!(h_data_var->GetEntries()==0)) leg->AddEntry(h_data_var,"data","L");
 leg->Draw();
 //Error for bkg
 //Bkg err
 double all_bkg_statErr_x[bin[v]]; double all_bkg_statErr_y[bin[v]]; double all_bkg_statErr_xerr[bin[v]]; double all_bkg_statErr_yerr[bin[v]];
 for(int b=0; b<bin[v]; b++){
  all_bkg_statErr_x[b] = 0; all_bkg_statErr_y[b] = 0; all_bkg_statErr_xerr[b] = 0; all_bkg_statErr_yerr[b] = 0;
  all_bkg_statErr_x[b] = h_sum_var->GetBinCenter(b+1);
  if(var=="asymbinvar" && doasym){
   all_bkg_statErr_xerr[b] = (asymbin[b+1]-asymbin[b])*0.5;
  }else{
   all_bkg_statErr_xerr[b] = ((endRange[v]-inRange[v])/bin[v])*0.5;
  }
  for(int r=0; r<rootplas_size; r++) all_bkg_statErr_yerr[b] += err_AllBkg[r][b]*err_AllBkg[r][b];
  all_bkg_statErr_y[b] = h_sum_var->GetBinContent(b+1); all_bkg_statErr_yerr[b] = sqrt(all_bkg_statErr_yerr[b]);
 }
 TGraphErrors *all_bkg_statErr = new TGraphErrors(bin[v], all_bkg_statErr_x, all_bkg_statErr_y, all_bkg_statErr_xerr, all_bkg_statErr_yerr);
 all_bkg_statErr->SetLineColor(kGray+3);
 all_bkg_statErr->SetFillStyle(3004);
 all_bkg_statErr->SetFillColor(kGray+3);
 all_bkg_statErr->Draw("E2same");
}




void draw_plot_withratio(TCanvas* c1, TH1F* h_data_var, THStack *hstack, TH1F* h_sum_var, TLegend *leg, int v, string var, string vartitlexaxis, int rootplas_size, double err_AllBkg[][col_size]){
 //Bottom plot
 TPad *c1_1 = new TPad("c1_1", "newpad",0.01,0.01,0.99,0.32);
 c1_1->Draw();
 c1_1->cd();
 c1_1->SetTopMargin(0.01);
 c1_1->SetBottomMargin(0.4);
 c1_1->SetRightMargin(0.01);
 c1_1->SetLeftMargin(0.125);
 //c1_1->SetFillStyle(0);
 double dataSUmc_x[bin[v]]; double dataSUmc_y[bin[v]]; double dataSUmc_xerr[bin[v]]; double dataSUmc_yerr[bin[v]];
 for(int j=0; j<bin[v]; j++){
  dataSUmc_x[j] = 0; dataSUmc_y[j] = 0; dataSUmc_xerr[j] = 0; dataSUmc_yerr[j] = 0;
  dataSUmc_x[j] = h_sum_var->GetBinCenter(j+1);  dataSUmc_xerr[j] = 0; 
  double mc_err = 0;
  for(int r=0; r<rootplas_size; r++) mc_err += err_AllBkg[r][j]*err_AllBkg[r][j];
  if(h_sum_var->GetBinContent(j+1)!=0){
   double rd = h_data_var->GetBinContent(j+1);
   double mc = h_sum_var->GetBinContent(j+1);
   dataSUmc_y[j]    = rd/mc; 
   dataSUmc_yerr[j] = sqrt(pow(sqrt(rd)/mc,2) + pow((rd*sqrt(mc_err))/(mc*mc),2));
  }else{
   dataSUmc_y[j]    = -1000000; 
   dataSUmc_yerr[j] = 0;
  }
 }
 if(var=="asymbinvar"){
  cout<<"const double sf[bin] = {";
  for(int ar=0; ar<bin[v]; ar++) cout<<dataSUmc_y[ar]<<",";
  cout<<"};"<<endl;
  cout<<"const double sferr[bin] = {";
  for(int ar=0; ar<bin[v]; ar++) cout<<dataSUmc_yerr[ar]<<",";
  cout<<"};"<<endl;
 }
 TGraphErrors *dataSUmc = new TGraphErrors(bin[v], dataSUmc_x, dataSUmc_y, dataSUmc_xerr, dataSUmc_yerr);
 dataSUmc->SetTitle(0);
 //dataSUmc->SetTitleSize(10);
 dataSUmc->SetMarkerStyle(7);
 //dataSUmc->SetMarkerColor(1);
 //dataSUmc->SetLineColor(1);
 dataSUmc->GetXaxis()->SetRangeUser(inRange[v],endRange[v]);
 dataSUmc->GetXaxis()->SetTitle(vartitlexaxis.c_str()); 
 dataSUmc->GetXaxis()->SetTitleSize(0.2);
 dataSUmc->GetYaxis()->SetTitle("Data/MC");
 dataSUmc->GetYaxis()->SetLabelSize(0.075);
 dataSUmc->GetYaxis()->SetTitleSize(0.15);
 dataSUmc->GetYaxis()->SetTitleOffset(0.35);
 dataSUmc->SetMinimum(0);  
 dataSUmc->SetMaximum(2);
 dataSUmc->GetXaxis()->SetRangeUser(inRange[v],endRange[v]);
 dataSUmc->GetXaxis()->SetLimits(inRange[v],endRange[v]);
 dataSUmc->Draw("APZ");
 dataSUmc->GetXaxis()->SetRangeUser(inRange[v],endRange[v]);
 dataSUmc->GetXaxis()->SetLimits(inRange[v],endRange[v]);
 TLine* line = new TLine(inRange[v],1,endRange[v],1);
 line->SetLineColor(kRed);
 line->SetLineWidth(2);
 line->Draw("same");
 //Top plots
 c1->cd();
 TPad *c1_2 = new TPad("c1_2", "newpad",0.01,0.33,0.99,0.99);
 if(ylogscale[v]==1) c1_2->SetLogy();
 c1_2->Draw();
 c1_2->cd();
 c1_2->SetTopMargin(0.065);
 c1_2->SetBottomMargin(0.1);
 c1_2->SetRightMargin(0.01);
 c1_2->SetLeftMargin(0.125);
 //c1_2->SetFillStyle(0);
}



/////
//   Set setTDRStyle_modified (from link https://twiki.cern.ch/twiki/pub/CMS/TRK10001/setTDRStyle_modified.C)
/////
void setTDRStyle(){
  TStyle *tdrStyle = new TStyle("tdrStyle","Style for P-TDR");
  // For the canvas:
  tdrStyle->SetCanvasBorderMode(0);
  tdrStyle->SetCanvasColor(kWhite);
  tdrStyle->SetCanvasDefH(600); //Height of canvas
  tdrStyle->SetCanvasDefW(600); //Width of canvas
  tdrStyle->SetCanvasDefX(0);   //POsition on screen
  tdrStyle->SetCanvasDefY(0);
  // For the Pad:
  tdrStyle->SetPadBorderMode(0);
  // tdrStyle->SetPadBorderSize(Width_t size = 1);
  tdrStyle->SetPadColor(kWhite);
  tdrStyle->SetPadGridX(false);
  tdrStyle->SetPadGridY(false);
  tdrStyle->SetGridColor(0);
  tdrStyle->SetGridStyle(3);
  tdrStyle->SetGridWidth(1);
  // For the frame:
  tdrStyle->SetFrameBorderMode(0);
  tdrStyle->SetFrameBorderSize(1);
  tdrStyle->SetFrameFillColor(0);
  tdrStyle->SetFrameFillStyle(0);
  tdrStyle->SetFrameLineColor(1);
  tdrStyle->SetFrameLineStyle(1);
  tdrStyle->SetFrameLineWidth(1);

  // For the histo:
  tdrStyle->SetHistFillColor(0);
  // tdrStyle->SetHistFillStyle(0);
  tdrStyle->SetHistLineColor(1);
  tdrStyle->SetHistLineStyle(0);
  tdrStyle->SetHistLineWidth(1);
  // tdrStyle->SetLegoInnerR(Float_t rad = 0.5);
  // tdrStyle->SetNumberContours(Int_t number = 20);
//  tdrStyle->SetEndErrorSize(0);
  tdrStyle->SetErrorX(0.);
//  tdrStyle->SetErrorMarker(20);

  tdrStyle->SetMarkerStyle(20);

  //For the fit/function:
  tdrStyle->SetOptFit(1);
  tdrStyle->SetFitFormat("5.4g");
  //tdrStyle->SetFuncColor(1);
  tdrStyle->SetFuncStyle(1);
  tdrStyle->SetFuncWidth(1);

  //For the date:
  tdrStyle->SetOptDate(0);
  // tdrStyle->SetDateX(Float_t x = 0.01);
  // tdrStyle->SetDateY(Float_t y = 0.01);

  // For the statistics box:
  tdrStyle->SetOptFile(0);
  tdrStyle->SetOptStat(""); // To display the mean and RMS:   SetOptStat("mr");
  //tdrStyle->SetStatColor(kWhite);
  tdrStyle->SetStatColor(kGray);
  tdrStyle->SetStatFont(42);

  tdrStyle->SetTextSize(11);
  tdrStyle->SetTextAlign(11);

  tdrStyle->SetStatTextColor(1);
  tdrStyle->SetStatFormat("6.4g");
  tdrStyle->SetStatBorderSize(0);
  tdrStyle->SetStatX(1.); //Starting position on X axis
  tdrStyle->SetStatY(1.); //Starting position on Y axis
  tdrStyle->SetStatFontSize(0.025); //Vertical Size
  tdrStyle->SetStatW(0.25); //Horizontal size 
  // tdrStyle->SetStatStyle(Style_t style = 1001);

  // Margins:
  tdrStyle->SetPadTopMargin(0.095);
  tdrStyle->SetPadBottomMargin(0.125);
  tdrStyle->SetPadLeftMargin(0.14);
  tdrStyle->SetPadRightMargin(0.1);

  // For the Global title:

  //  tdrStyle->SetOptTitle(0);
  tdrStyle->SetTitleFont(42);
  tdrStyle->SetTitleColor(1);
  tdrStyle->SetTitleTextColor(1);
  tdrStyle->SetTitleFillColor(10);
  tdrStyle->SetTitleFontSize(0.05);
  tdrStyle->SetTitleH(0.045); // Set the height of the title box
  //tdrStyle->SetTitleW(0); // Set the width of the title box
  tdrStyle->SetTitleX(0.20); // Set the position of the title box
  tdrStyle->SetTitleY(1.0); // Set the position of the title box
  // tdrStyle->SetTitleStyle(Style_t style = 1001);
  tdrStyle->SetTitleBorderSize(0);

  // For the axis titles:
  tdrStyle->SetTitleColor(1, "XYZ");
  tdrStyle->SetTitleFont(42, "XYZ");
  tdrStyle->SetTitleSize(0.06, "XYZ");
  // tdrStyle->SetTitleXSize(Float_t size = 0.02); // Another way to set the size?
  // tdrStyle->SetTitleYSize(Float_t size = 0.02);
  tdrStyle->SetTitleXOffset(0.9);
  tdrStyle->SetTitleYOffset(1.0);
  // tdrStyle->SetTitleOffset(1.1, "Y"); // Another way to set the Offset

  // For the axis labels:

  tdrStyle->SetLabelColor(1, "XYZ");
  tdrStyle->SetLabelFont(42, "XYZ");
  tdrStyle->SetLabelOffset(0.007, "XYZ");
  tdrStyle->SetLabelSize(0.05, "XYZ");

  // For the axis:

  tdrStyle->SetAxisColor(1, "XYZ");
  tdrStyle->SetStripDecimals(kTRUE);
  tdrStyle->SetTickLength(0.03, "XYZ");
  tdrStyle->SetNdivisions(510, "XYZ");
  tdrStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
  tdrStyle->SetPadTickY(1);

  // Change for log plots:
  tdrStyle->SetOptLogx(0);
  tdrStyle->SetOptLogy(0);
  tdrStyle->SetOptLogz(0);
  // Postscript options:
  // tdrStyle->SetPaperSize(15.,15.);
  // tdrStyle->SetLineScalePS(Float_t scale = 3);
  // tdrStyle->SetLineStyleString(Int_t i, const char* text);
  // tdrStyle->SetHeaderPS(const char* header);
  // tdrStyle->SetTitlePS(const char* pstitle);
  // tdrStyle->SetBarOffset(Float_t baroff = 0.5);
  // tdrStyle->SetBarWidth(Float_t barwidth = 0.5);
  // tdrStyle->SetPaintTextFormat(const char* format = "g");
  // tdrStyle->SetPalette(Int_t ncolors = 0, Int_t* colors = 0);
  // tdrStyle->SetTimeOffset(Double_t toffset);
  // tdrStyle->SetHistMinimumZero(kTRUE);

  tdrStyle->cd();
}
//const double xSecs[]    = {6025.2,    831.76};//4895. 
//const double evtsRead[] = {9002488., 4992231.};      
