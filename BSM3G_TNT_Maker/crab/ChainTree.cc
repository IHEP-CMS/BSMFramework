/**
This Macro
1. Uses a chain to merge the content of the trees that have been output by the analyzer and save it on a new rootpla 

Need to specify
0. See Declare constants
1. Do "voms-proxy-init --voms cms" if you read remote files
SingleMu: lcg-ls srm://stormfe1.pi.infn.it:8444/srm/managerv2?SFN=/cms/store/user/fromeo/SingleMuon
DrellYan: 
          lcg-ls srm://stormfe1.pi.infn.it:8444/srm/managerv2?SFN=/cms/store/user/fromeo/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8
         lcg-ls srm://stormfe1.pi.infn.it:8444/srm/managerv2?SFN=/cms/store/user/fromeo/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8
          lcg-ls srm://stormfe1.pi.infn.it:8444/srm/managerv2?SFN=/cms/store/user/fromeo/ZToMuMu_NNPDF30_13TeV-powheg_M_50_120
TTJets:   lcg-ls srm://stormfe1.pi.infn.it:8444/srm/managerv2?SFN=/cms/store/user/fromeo/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8        
tbarW:    lcg-ls srm://stormfe1.pi.infn.it:8444/srm/managerv2?SFN=/cms/store/user/fromeo/ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1
tW:       lcg-ls srm://stormfe1.pi.infn.it:8444/srm/managerv2?SFN=/cms/store/user/fromeo/ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1
WJets:    lcg-ls srm://stormfe1.pi.infn.it:8444/srm/managerv2?SFN=/cms/store/user/fromeo/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8
WW:       lcg-ls srm://stormfe1.pi.infn.it:8444/srm/managerv2?SFN=/cms/store/user/fromeo/WW_TuneCUETP8M1_13TeV-pythia8
WZ:       lcg-ls srm://stormfe1.pi.infn.it:8444/srm/managerv2?SFN=/cms/store/user/fromeo/WZ_TuneCUETP8M1_13TeV-pythia8
ZZ:       lcg-ls srm://stormfe1.pi.infn.it:8444/srm/managerv2?SFN=/cms/store/user/fromeo/ZZ_TuneCUETP8M1_13TeV-pythia8
Bari:     lcg-ls srm://storm-se-01.ba.infn.it:8444/srm/managerv2?SFN=/cms/store/user/fromeo/
TTJets:   lcg-ls srm://stormfe1.pi.infn.it:8444/srm/managerv2?SFN=/cms/store/user/fromeo/TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8
DrellYan: lcg-ls srm://stormfe1.pi.infn.it:8444/srm/managerv2?SFN=/cms/store/user/fromeo/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8
2. Do 0,$s/\/cms/root:\/\/xrootd-cms.infn.it\//g in the file specified in name_file
*/
/////
//   To run: root -l ChainTree.cc+  
/////
/////
//   Prepare Root and Roofit
/////
#include "TTree.h"
#include "TTreePlayer.h"
#include "TFile.h"
#include "TChain.h"
#include "TFileCollection.h"
#include <iostream>
using namespace std;
////
//   Declare constants
/////
const string tree_name    = "TNT/BOOM"; //The name of the tree you have defined in your rootplas
const string sample       = "TTJets";
const string name_file    = sample+".txt"; //List with the path of the rootfiles
const string name_rootple = sample+".root"; //The name of the new rootpla
/////
//   Main function
/////
void ChainTree(){
 //Put here the tree
 TTree *tree = new TTree("BOOM","BOOM");
 tree->SetMaxTreeSize(99000000000);
 //Create new file
 TFile *newfile = new TFile(name_rootple.c_str(),"recreate");
 newfile->cd();
 //Create chain, merge trees 
 TChain * chain      = new TChain(tree_name.c_str(),"");
 TFileCollection* fc = new TFileCollection("list", "list",name_file.c_str());
 chain->AddFileInfoList((TCollection*)fc->GetList());
 //Save it
 tree    = chain->CopyTree("");
 newfile = tree->GetCurrentFile();
 tree    = NULL;
 newfile->Write();
 newfile->Close();
 delete newfile; 
}
