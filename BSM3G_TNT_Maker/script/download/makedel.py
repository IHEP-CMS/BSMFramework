import os
import re
import sys
import commands
import subprocess

#Path you run this script
workpath = "/publicfs/cms/user/libh/jobs"


#Create the CopyFiles
##MC##
CopyScriptMCName = "delmc.sh"
CopyScriptMC      = file(CopyScriptMCName,"w")
##data##
CopyScriptDataName = "deldata.sh"
CopyScriptData      = file(CopyScriptDataName,"w")

#For Tier 2 Data
#PathT2 = '/pnfs/ihep.ac.cn/data/cms/store/user/fromeo'
PathT2 = '/pnfs/ihep.ac.cn/data/cms/store/user/binghuan'

#For Tier 3 Directory
mclocaldir = "/publicfs/cms/data/TopQuark/cms13TeV/FullMoriond2017/mc"
datalocaldir = "/publicfs/cms/data/TopQuark/cms13TeV/FullMoriond2017/data"
#mclocaldir = "/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/mc"
#datalocaldir = "/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/data"

mcchannels = [
##MC##
#Sig
#'ttHToNonbb_M125_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/crab_FullMorV1_ttHnobb/170216_105930/0000/',#'FullMorV1_ttHnobb', 54
##Bkg
#'TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/crab_FullMorV1_amcTTWJetsToLNuext2/170216_135729/0000/',#'FullMorV1_amcTTWJetsToLNuext2', 25
#'TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/crab_FullMorV1_amcTTWJetsToLNuext1/170216_140120/0000/',#'FullMorV1_amcTTWJetsToLNuext1', 17
#'TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8/crab_FullMorV1_amcTTZToLLNuNu_M-10_ext1/170216_140702/0000/',#'FullMorV1_amcTTZToLLNuNu_M-10_ext1', 18
#
#'TTZToLL_M-1to10_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_FullMorV1_TTZToLL_M1to10/170418_172802/0000/',#'FullMorV1_TTZToLL_M1to10',19
#'WGToLNuG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_FullMorV1_WGToLNuG_ext1/170418_173025/0000/',#'FullMorV1_WGToLNuG_ext1',145
#'ZGTo2LG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_FullMorV1_ZGTo2LG_ext1/170418_173222/0000/',#'FullMorV1_ZGTo2LG_ext1',210
#'TGJets_TuneCUETP8M1_13TeV_amcatnlo_madspin_pythia8/crab_FullMorV1_TGJets/170418_173441/0000/',#'FullMorV1_TGJets',19
#'TGJets_TuneCUETP8M1_13TeV_amcatnlo_madspin_pythia8/crab_FullMorV1_TGJets_ext1/170418_173659/0000/',#'FullMorV1_TGJets_ext1',37
#'TTGJets_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/crab_FullMorV1_TTGJets/170418_173858/0000/',#'FullMorV1_TTGJets', 80
#'WpWpJJ_EWK-QCD_TuneCUETP8M1_13TeV-madgraph-pythia8/crab_FullMorV1_WpWpJJ_EWK-QCD/170418_174133/0000/',#'FullMorV1_WpWpJJ_EWK-QCD',4
#'WWTo2L2Nu_DoubleScattering_13TeV-pythia8/crab_FullMorV1_WWTo2L2Nu_DS/170418_174405/0000/',#'FullMorV1_WWTo2L2Nu_DS',13
#'WWW_4F_TuneCUETP8M1_13TeV-amcatnlo-pythia8/crab_FullMorV1_WWW_4F/170418_174642/0000/',#'FullMorV1_WWW_4F',4
#'WWZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8/crab_FullMorV1_WWZ/170418_174841/0000/',#'FullMorV1_WWZ',4
#'WZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8/crab_FullMorV1_WZZ/170418_175041/0000/',#'FullMorV1_WZZ',5
#'ZZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8/crab_FullMorV1_ZZZ/170418_175243/0000/',#'FullMorV1_ZZZ',4
#'tZq_ll_4f_13TeV-amcatnlo-pythia8/crab_FullMorV1_tZq_ext1/170418_175440/0000/',#'FullMorV1_tZq_ext1',303
#'TTTT_TuneCUETP8M1_13TeV-amcatnlo-pythia8/crab_FullMorV1_TTTT/170418_175656/0000/',#'FullMorV1_TTTT',18
#'TTJets_SingleLeptFromTbar_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_FullMorV1_TTJets_sinLepTbar/170418_175902/0000/',#'FullMorV1_TTJets_sinLepTbar',196
#'TTJets_SingleLeptFromTbar_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_FullMorV1_TTJets_sinLepTbar_ext1/170418_180113/0000/',#'FullMorV1_TTJets_sinLepTbar_ext1',559
#'TTJets_SingleLeptFromT_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_FullMorV1_TTJets_sinLepT/170418_180402/0000/',#'FullMorV1_TTJets_sinLepT',167
#'TTJets_SingleLeptFromT_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_FullMorV1_TTJets_sinLepT_ext1/170418_180638/0000/',#'FullMorV1_TTJets_sinLepT_ext1',629
#'TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_FullMorV1_TTJets_diLep/170418_180842/0000/',#'FullMorV1_TTJets_diLep',83
#'TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_FullMorV1_TTJets_diLep_ext1/170418_181053/0000/',#'FullMorV1_TTJets_diLep_ext1',342
#'ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/crab_FullMorV1_powhegST/170418_181251/0000/',#'FullMorV1_powhegST',161
#'ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/crab_FullMorV1_powhegSaT/170418_181515/0000/',#'FullMorV1_powhegSaT',144
#'ST_t-channel_top_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1/crab_FullMorV1_powhegSTt/170418_181718/0000/',#'FullMorV1_powhegSTt',723
#'ST_t-channel_antitop_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1/crab_FullMorV1_powhegSTat/170418_181925/0000/',#'FullMorV1_powhegSTat',406
#'ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1/crab_FullMorV1_STs/170216_132119/0000/',#'FullMorV1_STs', 11
#'DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_FullMorV1_DY_M10to50/170418_182134/0000/',#'FullMorV1_DY_M10to50',555
#'DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_FullMorV1_DY_M50_ext1-v2/170418_182331/0000/',#'FullMorV1_DY_M50_ext1-v2',477
#'WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_FullMorV1_amcWJets/170418_182529/0000/',#'FullMorV1_amcWJets',300
#'WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8/crab_FullMorV1_WZTo3LNu/170418_182731/0000/',#'FullMorV1_WZTo3LNu',39
#'WWTo2L2Nu_13TeV-powheg/crab_FullMorV1_WWTo2L2Nu/170418_182933/0000/',#'FullMorV1_WWTo2L2Nu',24
#'ZZTo4L_13TeV_powheg_pythia8/crab_FullMorV1_ZZTo4L/170418_183130/0000/',#'FullMorV1_ZZTo4L',96
#
########
#'ttHTobb_M125_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/crab_FullMorV1_ttHbb/170216_130523/0000/',#'FullMorV1_ttHbb', 60
#'TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/crab_FullMorV1_TT/170216_130900/0000/',#'FullMorV1_TT', 493
#'TTToSemilepton_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/crab_FullMorV1_TTToSemilepton/170216_131212/0000/',#'FullMorV1_TTToSemilepton', 1260
#'TTToSemilepton_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/crab_FullMorV1_TTToSemilepton/170216_131212/0001/',#'FullMorV1_TTToSemilepton',  1260
#'TTToSemilepton_ttbbFilter_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/crab_FullMorV1_TTToSemilepton_ttbbFilter/170217_064724/0000/',#'FullMorV1_TTToSemilepton_ttbbFilter', 196
#'TTTo2L2Nu_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/crab_FullMorV1_TTTo2L2Nu/170216_131744/0000/',#'FullMorV1_TTTo2L2Nu', 732
#'ST_t-channel_top_4f_inclusiveDecays_TuneCUETP8M2T4_13TeV-powhegV2-madspin/crab_FullMorV1_STt/170217_064944/0000/',#'FullMorV1_STt', 61
#'ST_t-channel_antitop_4f_inclusiveDecays_TuneCUETP8M2T4_13TeV-powhegV2-madspin/crab_FullMorV1_SaTt/170216_132428/0000/',#'FullMorV1_SaTt', 84
#'ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M2T4/crab_FullMorV1_ST/170217_065140/0000/',#'FullMorV1_ST', 7
#'ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M2T4/crab_FullMorV1_SaT/170216_133251/0000/',#'FullMorV1_SaT', 9
#'WW_TuneCUETP8M1_13TeV-pythia8/crab_FullMorV1_WW/170216_133749/0000/',#'FullMorV1_WW', 6
#'WW_TuneCUETP8M1_13TeV-pythia8/crab_FullMorV1_WWext1/170216_134055/0000/',#'FullMorV1_WWext1', 58
#'WZ_TuneCUETP8M1_13TeV-pythia8/crab_FullMorV1_WZ/170217_065443/0000/',#'FullMorV1_WZ', 16
#'WZ_TuneCUETP8M1_13TeV-pythia8/crab_FullMorV1_WZext1/170216_134410/0000/',#'FullMorV1_WZext1', 64
#'ZZ_TuneCUETP8M1_13TeV-pythia8/crab_FullMorV1_ZZ/170216_134718/0000/',#'FullMorV1_ZZ', 6
#'ZZ_TuneCUETP8M1_13TeV-pythia8/crab_FullMorV1_ZZext1/170216_135129/0000/',#'FullMorV1_ZZext1', 33
#'TTWJetsToQQ_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/crab_FullMorV1_amcTTWJetsToQQ/170217_065640/0000/',#'FullMorV1_amcTTWJetsToQQ', 8
#'TTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8/crab_FullMorV1_amcTTZToQQ/170216_141058/0000/',#'FullMorV1_amcTTZToQQ' 10
#
#
#'TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8/crab_FullMorV1_amcTTZToLLNuNu_M-10_ext2/170322_193534/0000/',# 'FullMorV1_amcTTZToLLNuNu_M-10_ext2', 47
#'WJetsToLNu_HT-70To100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_FullMorV1_WJetsToLNu_HT-70To100/170322_193733/0000/',# 'FullMorV1_WJetsToLNu_HT-70To100', 104
#'WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_FullMorV1_WJetsToLNu_HT-100To200/170322_193931/0000/',# 'FullMorV1_WJetsToLNu_HT-100To200', 67
#'WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_FullMorV1_WJetsToLNu_HT-100To200_ext1/170322_194129/0000/',# 'FullMorV1_WJetsToLNu_HT-100To200_ext1', 151
#'WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_FullMorV1_WJetsToLNu_HT-100To200_ext2/170322_194326/0000/',# 'FullMorV1_WJetsToLNu_HT-100To200_ext2', 263
#'WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_FullMorV1_WJetsToLNu_HT-200To400/170322_194524/0000/',# 'FullMorV1_WJetsToLNu_HT-200To400', 45
#'WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_FullMorV1_WJetsToLNu_HT-200To400_ext1/170322_194723/0000/',# 'FullMorV1_WJetsToLNu_HT-200To400_ext1', 93
#'WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_FullMorV1_WJetsToLNu_HT-200To400_ext2/170322_194934/0000/',# 'FullMorV1_WJetsToLNu_HT-200To400_ext2', 176
#'WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_FullMorV1_WJetsToLNu_HT-400To600/170322_195144/0000/',# 'FullMorV1_WJetsToLNu_HT-400To600', 14
#'WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_FullMorV1_WJetsToLNu_HT-400To600_ext1/170322_195356/0000/',# 'FullMorV1_WJetsToLNu_HT-400To600_ext1', 25
#'WJetsToLNu_HT-600To800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_FullMorV1_WJetsToLNu_HT-600To800/170322_195555/0000/',# 'FullMorV1_WJetsToLNu_HT-600To800', 14
#'WJetsToLNu_HT-600To800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_FullMorV1_WJetsToLNu_HT-600To800_ext1/170322_195754/0000/',# 'FullMorV1_WJetsToLNu_HT-600To800_ext1', 97
#'WJetsToLNu_HT-800To1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_FullMorV1_WJetsToLNu_HT-800To1200/170322_195957/0000/',# 'FullMorV1_WJetsToLNu_HT-800To1200', 14
#'WJetsToLNu_HT-800To1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_FullMorV1_WJetsToLNu_HT-800To1200_ext1/170322_200201/0000/',# 'FullMorV1_WJetsToLNu_HT-800To1200_ext1', 58
#'WJetsToLNu_HT-1200To2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_FullMorV1_WJetsToLNu_HT-1200To2500/170322_202733/0000/',# 'FullMorV1_WJetsToLNu_HT-1200To2500', 4
#'WJetsToLNu_HT-1200To2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_FullMorV1_WJetsToLNu_HT-1200To2500_ext1/170322_202933/0000/',# 'FullMorV1_WJetsToLNu_HT-1200To2500_ext1', 58
#'WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_FullMorV1_WJetsToLNu_HT-2500ToInf/170322_203139/0000/',# 'FullMorV1_WJetsToLNu_HT-2500ToInf', 4
#'WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_FullMorV1_WJetsToLNu_HT-2500ToInf_ext1/170322_203341/0000/',# 'FullMorV1_WJetsToLNu_HT-2500ToInf_ext1', 23
#'DYJetsToLL_M-5to50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_FullMorV1_DY_M-5to50_HT-100to200/170322_203602/0000/',# 'FullMorV1_DY_M-5to50_HT-100to200', 7
#'DYJetsToLL_M-5to50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_FullMorV1_DY_M-5to50_HT-100to200_ext1/170322_203817/0000/',# 'FullMorV1_DY_M-5to50_HT-100to200_ext1', 46  
#'DYJetsToLL_M-5to50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_FullMorV1_DY_M-5to50_HT-200to400/170322_204027/0000/',# 'FullMorV1_DY_M-5to50_HT-200to400',  13
#'DYJetsToLL_M-5to50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_FullMorV1_DY_M-5to50_HT-200to400_ext1/170322_204234/0000/',# 'FullMorV1_DY_M-5to50_HT-200to400_ext1', 27  
#'DYJetsToLL_M-5to50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_FullMorV1_DY_M-5to50_HT-400to600/170322_204442/0000/',# 'FullMorV1_DY_M-5to50_HT-400to600',  10
#'DYJetsToLL_M-5to50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_FullMorV1_DY_M-5to50_HT-400to600_ext1/170322_204652/0000/',# 'FullMorV1_DY_M-5to50_HT-400to600_ext1', 18  
#'DYJetsToLL_M-5to50_HT-600toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_FullMorV1_DY_M-5to50_HT-600toInf/170322_204902/0000/',# 'FullMorV1_DY_M-5to50_HT-600toInf',  9
#'DYJetsToLL_M-5to50_HT-600toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_FullMorV1_DY_M-5to50_HT-600toInf_ext1/170322_205117/0000/',# 'FullMorV1_DY_M-5to50_HT-100toInf_ext1', 24 //typo is 600toInf_ext1  
#'DYJetsToLL_M-50_HT-70to100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_FullMorV1_DY_M-50_HT-70to100/170322_205323/0000/',# 'FullMorV1_DY_M-50_HT-70to100',  56
#'DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_FullMorV1_DY_M-50_HT-100to200/170322_205536/0000/',# 'FullMorV1_DY_M-50_HT-100to200',  20
#'DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_FullMorV1_DY_M-50_HT-100to200_ext1/170322_210618/0000/',# 'FullMorV1_DY_M-50_HT-100to200_ext1', 68  
#'DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_FullMorV1_DY_M-50_HT-200to400/170322_210826/0000/',# 'FullMorV1_DY_M-50_HT-200to400',  20
#'DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_FullMorV1_DY_M-50_HT-200to400_ext1/170322_211035/0000/',# 'FullMorV1_DY_M-50_HT-200to400_ext1', 65  
#'DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_FullMorV1_DY_M-50_HT-400to600/170322_211246/0000/',# 'FullMorV1_DY_M-50_HT-400to600',  10
#'DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_FullMorV1_DY_M-50_HT-400to600_ext1/170322_211450/0000/',# 'FullMorV1_DY_M-50_HT-400to600_ext1', 83  
#'DYJetsToLL_M-50_HT-600to800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_FullMorV1_DY_M-50_HT-600to800_v2/170322_211713/0000/',# 'FullMorV1_DY_M-50_HT-600to800_v2',  176 
#'DYJetsToLL_M-50_HT-800to1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_FullMorV1_DY_M-50_HT-800to1200/170322_211917/0000/',# 'FullMorV1_DY_M-50_HT-800to1200',  41
#'DYJetsToLL_M-50_HT-1200to2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_FullMorV1_DY_M-50_HT-1200to2500/170322_212129/0000/',# 'FullMorV1_DY_M-50_HT-1200to2500',  17
#'DYJetsToLL_M-50_HT-2500toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_FullMorV1_DY_M-50_HT-2500toInf/170322_212341/0000/',# 'FullMorV1_DY_M-50_HT-2500toInf',  9
'BstarToTW_M-1200_LH_TuneCUETP8M1_13TeV-madgraph-pythia8',
'BstarToTW_M-1400_LH_TuneCUETP8M1_13TeV-madgraph-pythia8',
'BstarToTW_M-1600_LH_TuneCUETP8M1_13TeV-madgraph-pythia8',
'BstarToTW_M-1800_LH_TuneCUETP8M1_13TeV-madgraph-pythia8',
'BstarToTW_M-2000_LH_TuneCUETP8M1_13TeV-madgraph-pythia8',
'BstarToTW_M-2200_LH_TuneCUETP8M1_13TeV-madgraph-pythia8',
'BstarToTW_M-2400_LH_TuneCUETP8M1_13TeV-madgraph-pythia8',
'BstarToTW_M-2600_LH_TuneCUETP8M1_13TeV-madgraph-pythia8',
'BstarToTW_M-2800_LH_TuneCUETP8M1_13TeV-madgraph-pythia8',
'BstarToTW_M-3000_LH_TuneCUETP8M1_13TeV-madgraph-pythia8',
'BstarToTW_M-1200_RH_TuneCUETP8M1_13TeV-madgraph-pythia8',
'BstarToTW_M-1400_RH_TuneCUETP8M1_13TeV-madgraph-pythia8',
'BstarToTW_M-1600_RH_TuneCUETP8M1_13TeV-madgraph-pythia8',
'BstarToTW_M-1800_RH_TuneCUETP8M1_13TeV-madgraph-pythia8',
'BstarToTW_M-2000_RH_TuneCUETP8M1_13TeV-madgraph-pythia8',
'BstarToTW_M-2200_RH_TuneCUETP8M1_13TeV-madgraph-pythia8',
'BstarToTW_M-2400_RH_TuneCUETP8M1_13TeV-madgraph-pythia8',
'BstarToTW_M-2600_RH_TuneCUETP8M1_13TeV-madgraph-pythia8',
'BstarToTW_M-2800_RH_TuneCUETP8M1_13TeV-madgraph-pythia8',
'BstarToTW_M-3000_RH_TuneCUETP8M1_13TeV-madgraph-pythia8',
]

datachannels = [

## Data##

#'SingleElectron/crab_FullMorV1_SEleBlockB1/170216_125915/0000/', #'FullMorV1_SEleBlockB1', 1758
#'SingleElectron/crab_FullMorV1_SEleBlockB1/170216_125915/0001/', #'FullMorV1_SEleBlockB1', 1758
#'SingleElectron/crab_FullMorV1_SEleBlockC1/170216_130855/0000/', #'FullMorV1_SEleBlockC1', 580
#'SingleElectron/crab_FullMorV1_SEleBlockD1/170216_131209/0000/', #'FullMorV1_SEleBlockD1', 972
#'SingleElectron/crab_FullMorV1_SEleBlockE1/170216_131445/0000/', #'FullMorV1_SEleBlockE1', 826
#'SingleElectron/crab_FullMorV1_SEleBlockF1/170216_131713/0000/', #'FullMorV1_SEleBlockF1', 523
#'SingleElectron/crab_FullMorV1_SEleBlockF2/170216_132117/0000/', #'FullMorV1_SEleBlockF2', 81
#'SingleElectron/crab_FullMorV1_SEleBlockG1/170216_132450/0000/', #'FullMorV1_SEleBlockG1', 1423
#'SingleElectron/crab_FullMorV1_SEleBlockG1/170216_132450/0001/', #'FullMorV1_SEleBlockG1', 1423
#'SingleElectron/crab_FullMorV1_SEleBlockH1/170216_132719/0000/', #'FullMorV1_SEleBlockH1', 1541
#'SingleElectron/crab_FullMorV1_SEleBlockH1/170216_132719/0001/', #'FullMorV1_SEleBlockH1', 1541
#'SingleElectron/crab_FullMorV1_SEleBlockH2/170217_142948/0000/',#'FullMorV1_SEleBlockH2', 41
#'SingleMuon/crab_FullMorV1_SMuBlockB1/170216_133716/0000/', #'FullMorV1_SMuBlockB1',  1756
#'SingleMuon/crab_FullMorV1_SMuBlockB1/170216_133716/0001/', #'FullMorV1_SMuBlockB1',  1756
#'SingleMuon/crab_FullMorV1_SMuBlockC1/170217_102735/0000/', #'FullMorV1_SMuBlockC1',  580
#'SingleMuon/crab_FullMorV1_SMuBlockD1/170216_134614/0000/', #'FullMorV1_SMuBlockD1',  972
#'SingleMuon/crab_FullMorV1_SMuBlockE1/170217_081050/0000/', #'FullMorV1_SMuBlockE1',  826
#'SingleMuon/crab_FullMorV1_SMuBlockF1/170217_081305/0000/', #'FullMorV1_SMuBlockF1',  523
#'SingleMuon/crab_FullMorV1_SMuBlockF2/170216_135711/0000/', #'FullMorV1_SMuBlockF2',  81
#'SingleMuon/crab_FullMorV1_SMuBlockG1/170216_140118/0000/', #'FullMorV1_SMuBlockG1',  1423
#'SingleMuon/crab_FullMorV1_SMuBlockG1/170216_140118/0001/', #'FullMorV1_SMuBlockG1',  1423
#'SingleMuon/crab_FullMorV1_SMuBlockH1/170216_140404/0000/', #'FullMorV1_SMuBlockH1',  1541
#'SingleMuon/crab_FullMorV1_SMuBlockH1/170216_140404/0001/', #'FullMorV1_SMuBlockH1',  1541
#'SingleMuon/crab_FullMorV1_SMuBlockH2/170216_140830/0000/', #'FullMorV1_SMuBlockH2',  41





]

for mcchannel in mcchannels:

        mclocalFiledir = mclocaldir + "/"+mcchannel
        print >> CopyScriptMC, "env -i X509_USER_PROXY=/tmp/x509up_u$UID gfal-rm -r srm://srm.ihep.ac.cn:8443/srm/managerv2?SFN="+PathT2+"/"+mcchannel

for datachannel in datachannels:

        datalocalFiledir = datalocaldir + "/"+datachannel
        print >> CopyScriptData, "env -i X509_USER_PROXY=/tmp/x509up_u$UID gfal-rm -r srm://srm.ihep.ac.cn:8443/srm/managerv2?SFN="+PathT2+"/"+datachannel





os.popen('chmod +x '+CopyScriptMCName)
os.popen('chmod +x '+CopyScriptDataName)
