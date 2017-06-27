import os
import re
import sys
import commands
import subprocess
#import xmlrpclib

#os.command("grid-proxy-init")
print ''
print '###############    Remember to setup grid vo first'
print '###############    voms-proxy-init --voms cms'
print ''
##For Electron data and MC
#se = "cmsse02.na.infn.it"
#port = "8446"
#path = "/dpm/na.infn.it/home/cms/store/user/oiorio/2012/Fall12/"

#For Muon data
#se = "stormfe1.pi.infn.it"
#port = "8444"
#path = "/cms/store/user/mmerola/SingleTop/Analysis/2014/"


#dir ='"srm://'+se + ":" + port + "/srm/managerv2?SFN=" + path +'"'
#command_ls = 'lcg-ls -l -b -D srmv2 -T srmv2 '+ dir


#For Local data
#path = '/pnfs/ihep.ac.cn/data/cms/store/user/aspiezia'
path = '/pnfs/ihep.ac.cn/data/cms/store/user/binghuan'
command_ls = 'ls -ltr ' + path + '"'


localdir = "/publicfs/cms/data/TopQuark/cms13TeV/FullMoriond2017/mc"
#localdir = "/publicfs/cms/data/TopQuark/tWsemiLeptonic/july2013/"#
#localdir = "/publicfs/cms/data/TopQuark/tWsemiLeptonic/April2014"
#localdir = "/publicfs/cms/data/TopQuark/cms13TeV/Samples8X/mc"
#localdir = "/publicfs/cms/data/TopQuark/cms13TeV/Samples2607/data"

nSimultaneous = 35
#nSimultaneous 
channels = [

'TTZToLL_M-1to10_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_FullMorV1_TTZToLL_M1to10/170418_172802/0000/',#'FullMorV1_TTZToLL_M1to10',19
'WGToLNuG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_FullMorV1_WGToLNuG_ext1/170418_173025/0000/',#'FullMorV1_WGToLNuG_ext1',145
'ZGTo2LG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_FullMorV1_ZGTo2LG_ext1/170418_173222/0000/',#'FullMorV1_ZGTo2LG_ext1',210
'TGJets_TuneCUETP8M1_13TeV_amcatnlo_madspin_pythia8/crab_FullMorV1_TGJets/170418_173441/0000/',#'FullMorV1_TGJets',19
'TGJets_TuneCUETP8M1_13TeV_amcatnlo_madspin_pythia8/crab_FullMorV1_TGJets_ext1/170418_173659/0000/',#'FullMorV1_TGJets_ext1',37
'TTGJets_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/crab_FullMorV1_TTGJets/170418_173858/0000/',#'FullMorV1_TTGJets', 80
'WpWpJJ_EWK-QCD_TuneCUETP8M1_13TeV-madgraph-pythia8/crab_FullMorV1_WpWpJJ_EWK-QCD/170418_174133/0000/',#'FullMorV1_WpWpJJ_EWK-QCD',4
'WWTo2L2Nu_DoubleScattering_13TeV-pythia8/crab_FullMorV1_WWTo2L2Nu_DS/170418_174405/0000/',#'FullMorV1_WWTo2L2Nu_DS',13
'WWW_4F_TuneCUETP8M1_13TeV-amcatnlo-pythia8/crab_FullMorV1_WWW_4F/170418_174642/0000/',#'FullMorV1_WWW_4F',4
'WWZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8/crab_FullMorV1_WWZ/170418_174841/0000/',#'FullMorV1_WWZ',4
'WZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8/crab_FullMorV1_WZZ/170418_175041/0000/',#'FullMorV1_WZZ',5
'ZZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8/crab_FullMorV1_ZZZ/170418_175243/0000/',#'FullMorV1_ZZZ',4
'tZq_ll_4f_13TeV-amcatnlo-pythia8/crab_FullMorV1_tZq_ext1/170418_175440/0000/',#'FullMorV1_tZq_ext1',303
'TTTT_TuneCUETP8M1_13TeV-amcatnlo-pythia8/crab_FullMorV1_TTTT/170418_175656/0000/',#'FullMorV1_TTTT',18
'TTJets_SingleLeptFromTbar_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_FullMorV1_TTJets_sinLepTbar/170418_175902/0000/',#'FullMorV1_TTJets_sinLepTbar',196
'TTJets_SingleLeptFromTbar_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_FullMorV1_TTJets_sinLepTbar_ext1/170418_180113/0000/',#'FullMorV1_TTJets_sinLepTbar_ext1',559
'TTJets_SingleLeptFromT_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_FullMorV1_TTJets_sinLepT/170418_180402/0000/',#'FullMorV1_TTJets_sinLepT',167
'TTJets_SingleLeptFromT_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_FullMorV1_TTJets_sinLepT_ext1/170418_180638/0000/',#'FullMorV1_TTJets_sinLepT_ext1',629
'TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_FullMorV1_TTJets_diLep/170418_180842/0000/',#'FullMorV1_TTJets_diLep',83
'TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_FullMorV1_TTJets_diLep_ext1/170418_181053/0000/',#'FullMorV1_TTJets_diLep_ext1',342
'ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/crab_FullMorV1_powhegST/170418_181251/0000/',#'FullMorV1_powhegST',161
'ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/crab_FullMorV1_powhegSaT/170418_181515/0000/',#'FullMorV1_powhegSaT',144
'ST_t-channel_top_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1/crab_FullMorV1_powhegSTt/170418_181718/0000/',#'FullMorV1_powhegSTt',723
'ST_t-channel_antitop_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1/crab_FullMorV1_powhegSTat/170418_181925/0000/',#'FullMorV1_powhegSTat',406
'DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_FullMorV1_DY_M10to50/170418_182134/0000/',#'FullMorV1_DY_M10to50',555
'DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_FullMorV1_DY_M50_ext1-v2/170418_182331/0000/',#'FullMorV1_DY_M50_ext1-v2',477
'WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_FullMorV1_amcWJets/170418_182529/0000/',#'FullMorV1_amcWJets',300
'WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8/crab_FullMorV1_WZTo3LNu/170418_182731/0000/',#'FullMorV1_WZTo3LNu',39
'WWTo2L2Nu_13TeV-powheg/crab_FullMorV1_WWTo2L2Nu/170418_182933/0000/',#'FullMorV1_WWTo2L2Nu',24
'ZZTo4L_13TeV_powheg_pythia8/crab_FullMorV1_ZZTo4L/170418_183130/0000/',#'FullMorV1_ZZTo4L',96
##################          MC   ######################
#'ttHTobb_M125_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/crab_FullMorV1_ttHbb/170216_130523/0000/',#'FullMorV1_ttHbb', 60
#'ttHToNonbb_M125_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/crab_FullMorV1_ttHnobb/170216_105930/0000/',#'FullMorV1_ttHnobb', 54
#'TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/crab_FullMorV1_TT/170216_130900/0000/',#'FullMorV1_TT', 493
#'TTToSemilepton_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/crab_FullMorV1_TTToSemilepton/170216_131212/0000/',#'FullMorV1_TTToSemilepton', 1260
#'TTToSemilepton_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/crab_FullMorV1_TTToSemilepton/170216_131212/0001/',#'FullMorV1_TTToSemilepton',  1260
#'TTToSemilepton_ttbbFilter_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/crab_FullMorV1_TTToSemilepton_ttbbFilter/170217_064724/0000/',#'FullMorV1_TTToSemilepton_ttbbFilter', 196
#'TTTo2L2Nu_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/crab_FullMorV1_TTTo2L2Nu/170216_131744/0000/',#'FullMorV1_TTTo2L2Nu', 732
#'ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1/crab_FullMorV1_STs/170216_132119/0000/',#'FullMorV1_STs', 11
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
#'TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/crab_FullMorV1_amcTTWJetsToLNuext2/170216_135729/0000/',#'FullMorV1_amcTTWJetsToLNuext2', 25
#'TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/crab_FullMorV1_amcTTWJetsToLNuext1/170216_140120/0000/',#'FullMorV1_amcTTWJetsToLNuext1', 17
#'TTWJetsToQQ_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/crab_FullMorV1_amcTTWJetsToQQ/170217_065640/0000/',#'FullMorV1_amcTTWJetsToQQ', 8
#'TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8/crab_FullMorV1_amcTTZToLLNuNu_M-10_ext1/170216_140702/0000/',#'FullMorV1_amcTTZToLLNuNu_M-10_ext1', 18
#'TTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8/crab_FullMorV1_amcTTZToQQ/170216_141058/0000/',#'FullMorV1_amcTTZToQQ' 10
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
##'DYJetsToLL_M-5to50_HT-600toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_FullMorV1_DY_M-5to50_HT-100toInf_ext1/170322_205117/0000/',# 'FullMorV1_DY_M-5to50_HT-100toInf_ext1', 24 //typo is 600toInf_ext1  
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

















###Before 2017 Moriond ####
#'TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/crab_TTHLep8020_TTW_13Feb_V2/170213_145721/0000/',
#'ttHToNonbb_M125_13TeV_powheg_pythia8/crab_TTHLep8020_TTHnbb_13Feb_V2/170213_150859/0000/',
#'TT_TuneCUETP8M1_13TeV-powheg-pythia8/crab_TTHLep8017SIHEP_TT/161203_115309/0000/',
#'TT_TuneCUETP8M1_13TeV-powheg-pythia8/crab_TTHLep8017SIHEP_TT/161203_115309/0001/',
#'TT_TuneCUETP8M1_13TeV-powheg-pythia8/crab_TTHLep8017SIHEP_TT/161203_115309/0002/',
#'TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/crab_Full2202_TTW/160728_122529/0000',
#'TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8/crab_Full2202_TTZ/160728_122555/0000',
#'ttHTobb_M125_13TeV_powheg_pythia8/crab_First2016_TTHbb/160726_140559/0000',
#'TT_TuneCUETP8M1_13TeV-powheg-pythia8/crab_First2016_TT/160726_140659/0000',
#'TT_TuneCUETP8M1_13TeV-powheg-pythia8/crab_First2016_TT/160726_140659/0001',
#'TT_TuneCUETP8M1_13TeV-powheg-pythia8/crab_First2016_TT/160726_140659/0002',
#'ttHToNonbb_M125_13TeV_powheg_pythia8/crab_First2016_TTHnbb/160726_140630/0000',
#'WZ_TuneCUETP8M1_13TeV-pythia8/crab_First2016_WZ/160726_141004/0000',
#'DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_First2016_2_DY/160726_151834/0000',
#'DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_First2016_2_DY/160726_151834/0001',
#'ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/crab_First2016_2_ST/160726_151730/0000',
#'ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/crab_First2016_2_SaT/160726_151805/0000',
#'WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_First2016_2_WJets/160726_151923/0000',
#'ZZ_TuneCUETP8M1_13TeV-pythia8/crab_First2016_2_ZZ/160726_152149/0000',
#'WW_TuneCUETP8M1_13TeV-pythia8/crab_First2016_2_WW/160726_151952/0000',

#'ttHToNonbb_M125_13TeV_powheg_pythia8/crab_Full2202_TTHnbbMass/160517_140907/0000',
#'ZprimeToTauTau_M_5000_TuneCUETP8M1_tauola_13TeV_pythia8/crab_Full2202_Zp5000/160301_113101/0000',
#'ZprimeToTauTau_M_4500_TuneCUETP8M1_tauola_13TeV_pythia8/crab_Full2202_Zp4500/160301_113037/0000',
#'ZprimeToTauTau_M_4000_TuneCUETP8M1_tauola_13TeV_pythia8/crab_Full2202_Zp4000/160301_113013/0000',
#'ZprimeToTauTau_M_2500_TuneCUETP8M1_tauola_13TeV_pythia8/crab_Full2202_Zp2500/160301_112546/0000',
#'ZprimeToTauTau_M_2000_TuneCUETP8M1_tauola_13TeV_pythia8/crab_Full2202_Zp2000/160301_112522/0000',
#'ZprimeToTauTau_M_1500_TuneCUETP8M1_tauola_13TeV_pythia8/crab_Full2202_Zp1500/160301_112457/0000',
#'ZprimeToTauTau_M_1000_TuneCUETP8M1_tauola_13TeV_pythia8/crab_Full2202_Zp1000/160301_112426/0000',
#'ZprimeToTauTau_M_500_TuneCUETP8M1_tauola_13TeV_pythia8/crab_Full2202_Zp500/160301_112354/0000',
#"ttHTobb_M125_13TeV_powheg_pythia8/crab_Full2202_TTHbb/160222_223411/0000",
#"ttHToNonbb_M125_13TeV_powheg_pythia8/crab_Full2202_TTHnbb/160222_223430/0000",
#"ZZ_TuneCUETP8M1_13TeV-pythia8/crab_Full2202_ZZ/160222_224032/0000",
#"WZ_TuneCUETP8M1_13TeV-pythia8/crab_Full2202_WZ/160222_224007/0000",
#"ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/crab_Full2202_SaT/160222_223547/0000",
#"WW_TuneCUETP8M1_13TeV-pythia8/crab_Full2202_WW/160222_223947/0000",
#"ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/crab_Full2202_ST/160222_223524/0000",
#"ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1/crab_Full2202_STs/160222_223608/0000",
#"QCD_HT100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_Full2202_QCD_HT100to200/160222_224054/0000",
#"QCD_HT100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_Full2202_QCD_HT100to200/160222_224054/0001",
#"QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_Full2202_QCD_HT300to500/160222_224202/0000",
#"QCD_HT200to300_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_Full2202_QCD_HT200to300/160222_224121/0000",
#"QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_Full2202_QCD_HT500to700/160222_224225/0000",
#"QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_Full2202_QCD_HT700to1000/160222_224245/0000",
#"ZprimeToTauTau_M_1500_TuneCUETP8M1_tauola_13TeV_pythia8/crab_HN2Tau3101_Zp1500/160131_162759/0000",
#"ZprimeToTauTau_M_2000_TuneCUETP8M1_tauola_13TeV_pythia8/crab_HN2Tau3101_Zp2000/160131_162817/0000",
#"QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_Full2202_QCD_HT1000to1500/160222_224313/0000",
#"ST_t-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1/crab_Full2202_STtext1/160222_223658/0000",
#"QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_Full2202_QCD_HT1500to2000/160222_224332/0000",
#"QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_Full2202_QCD_HT2000toInf/160222_224406/0000",

#"DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_Full2202_DY/160222_223837/0000",
#"DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_Full2202_amcDY1050/160222_223718/0000",
#"DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_Full2202_amcDY/160222_223759/0000",
#"DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_Full2202_amcDY1050ext1/160222_223740/0000",
#"DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_Full2202_amcDY1050ext1/160222_223740/0001",
#"TT_TuneCUETP8M1_13TeV-powheg-pythia8/crab_Full2202_TT/160222_223457/0000",
#"TT_TuneCUETP8M1_13TeV-powheg-pythia8/crab_Full2202_TT/160222_223457/0001",
#"TT_TuneCUETP8M1_13TeV-powheg-pythia8/crab_Full2202_TT/160222_223457/0002",
#"ST_t-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1/crab_Full2202_STt/160222_223629/0000",
#"ST_t-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1/crab_Full2202_STtext1/160222_223658/0000",

#"WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_Full2202_WJets/160222_223915/0000",
#"WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_Full2202_WJets/160222_223915/0001",
#"WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_Full2202_amcWJets/160222_223856/0000",
###################  data  ##############
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



#'SingleMuon/crab_First2016_SMuB/160726_142254/0000',
#'SingleMuon/crab_First2016_SMuB/160726_142254/0001',
#'SingleMuon/crab_First2016_SMuC/160726_142334/0000',
#'SingleElectron/crab_Firts2016_SEleB/160726_142148/0000',
#'SingleElectron/crab_Firts2016_SEleB/160726_142148/0001',
#'SingleElectron/crab_Firts2016_SEleC/160726_142218/0000',

#'SingleElectron/crab_Full2202_SEle_16Dec2015_JsonGold/160225_143951/0000',
#'SingleElectron/crab_Full2202_SEle_16Dec2015_JsonGold/160225_143951/0001',
#'SingleMuon/crab_Full2202_SMu_16Dec2015S_JsonGold/160224_170442/0000',
#'SingleMuon/crab_Full2202_SMu_16Dec2015S_JsonGold/160224_170442/0001',
#'DoubleEG/crab_Full2202_DEG_16Dec2015S_JsonGold/160224_170506/0000',
#'DoubleEG/crab_Full2202_DEG_16Dec2015S_JsonGold/160224_170506/0001',
#'DoubleMuon/crab_Full2202_DMu_16Dec2015_JsonGold/160224_170549/0000',
#'DoubleMuon/crab_Full2202_DMu_16Dec2015_JsonGold/160224_170549/0001',
#'MuonEG/crab_Full2202_MuEG_16Dec2015_JsonGold/160224_170528/0000',
#'MuonEG/crab_Full2202_MuEG_16Dec2015_JsonGold/160224_170528/0001',
#"ExtendedWeakIsospinModel_mumujj_L5000_M500_CalcHEP/crab_Full2202_mumujj_L5000_M500/160222_225730/0000",
#"DoubleEG/crab_Full2202_DEG_16Dec2015S_JsonSilv/160222_224854/0000",
#"DoubleEG/crab_Full2202_DEG_16Dec2015S_JsonSilv/160222_224854/0001",
#"MuonEG/crab_Full2202_MuEG_16Dec2015_JsonSilv/160222_224932/0000",
#"MuonEG/crab_Full2202_MuEG_16Dec2015_JsonSilv/160222_224932/0001",
#"DoubleMuon/crab_Full2202_DMu_16Dec2015_JsonSilv/160222_224952/0000",
#"DoubleMuon/crab_Full2202_DMu_16Dec2015_JsonSilv/160222_224952/0001",
#"SingleElectron/crab_Full2202_SEle_16Dec2015_JsonSilv/160222_224812/0000",
#"SingleElectron/crab_Full2202_SEle_16Dec2015_JsonSilv/160222_224812/0001",
#"SingleMuon/crab_Full2202_SMu_16Dec2015S_JsonSilv/160222_224832/0000",
#"SingleMuon/crab_Full2202_SMu_16Dec2015S_JsonSilv/160222_224832/0001",
#"ExtendedWeakIsospinModel_eejj_L15000_M1500_CalcHEP/crab_Full2202_eejj_L15000_M1500/160222_225457/0000",
#"ExtendedWeakIsospinModel_eejj_L15000_M2500_CalcHEP/crab_Full2202_eejj_L15000_M2500/160222_225515/0000",
#"ExtendedWeakIsospinModel_eejj_L15000_M3500_CalcHEP/crab_Full2202_eejj_L15000_M3500/160222_225534/0000",
#"ExtendedWeakIsospinModel_eejj_L15000_M4500_CalcHEP/crab_Full2202_eejj_L15000_M4500/160222_225553/0000",
#"ExtendedWeakIsospinModel_eejj_L15000_M500_CalcHEP/crab_Full2202_eejj_L15000_M500/160222_225438/0000",
#"ExtendedWeakIsospinModel_eejj_L25000_M2500_CalcHEP/crab_Full2202_eejj_L25000_M2500/160222_225631/0000",
#"ExtendedWeakIsospinModel_eejj_L25000_M3500_CalcHEP/crab_Full2202_eejj_L25000_M3500/160222_225651/0000",
#"ExtendedWeakIsospinModel_eejj_L25000_M4500_CalcHEP/crab_Full2202_eejj_L25000_M4500/160222_225711/0000",
#"ExtendedWeakIsospinModel_eejj_L25000_M500_CalcHEP/crab_Full2202_eejj_L25000_M500/160222_225612/0000",
#"ExtendedWeakIsospinModel_eejj_L5000_M1500_CalcHEP/crab_Full2202_eejj_L5000_M1500/160222_225319/0000",
#"ExtendedWeakIsospinModel_eejj_L5000_M2500_CalcHEP/crab_Full2202_eejj_L5000_M2500/160222_225339/0000",
#"ExtendedWeakIsospinModel_eejj_L5000_M3500_CalcHEP/crab_Full2202_eejj_L5000_M3500/160222_225358/0000",
#"ExtendedWeakIsospinModel_eejj_L5000_M4500_CalcHEP/crab_Full2202_eejj_L5000_M4500/160222_225417/0000",
#"ExtendedWeakIsospinModel_eejj_L5000_M500_CalcHEP/crab_Full2202_eejj_L5000_M500/160222_225300/0000",
##"ExtendedWeakIsospinModel_mumujj_L15000_M1500_CalcHEP/crab_Full2202_mumujj_L15000_M1500/160222_225924/0000",
#"ExtendedWeakIsospinModel_mumujj_L15000_M2500_CalcHEP/crab_Full2202_mumujj_L15000_M2500/160222_225943/0000",
#"ExtendedWeakIsospinModel_mumujj_L15000_M3500_CalcHEP/crab_Full2202_mumujj_L15000_M3500/160222_230002/0000",
#"ExtendedWeakIsospinModel_mumujj_L15000_M4500_CalcHEP/crab_Full2202_mumujj_L15000_M4500/160222_230024/0000",
#"ExtendedWeakIsospinModel_mumujj_L15000_M500_CalcHEP/crab_HN1502bj_mumujj_L15000_M500/160216_125922/0000",
#"ExtendedWeakIsospinModel_mumujj_L25000_M1500_CalcHEP/crab_Full2202_mumujj_L25000_M1500/160222_230328/0000",
#"ExtendedWeakIsospinModel_mumujj_L25000_M4500_CalcHEP/crab_Full2202_mumujj_L25000_M4500/160222_230433/0000",
#"ExtendedWeakIsospinModel_mumujj_L25000_M500_CalcHEP/crab_Full2202_mumujj_L25000_M500/160222_230043/0000",
#"ExtendedWeakIsospinModel_mumujj_L5000_M1500_CalcHEP/crab_Full2202_mumujj_L5000_M1500/160222_225749/0000",
#"ExtendedWeakIsospinModel_mumujj_L5000_M2500_CalcHEP/crab_Full2202_mumujj_L5000_M2500/160222_225808/0000",
#"ExtendedWeakIsospinModel_mumujj_L5000_M3500_CalcHEP/crab_Full2202_mumujj_L5000_M3500/160222_225827/0000",
#"ExtendedWeakIsospinModel_mumujj_L5000_M4500_CalcHEP/crab_Full2202_mumujj_L5000_M4500/160222_225846/0000",
#"ExtendedWeakIsospinModel_mumujj_L5000_M500_CalcHEP/crab_Full2202_mumujj_L5000_M500/160222_225730/0000",

#"MC/SChannel",
#"MC/SbarChannel",
#"MC/TChannel",
#"MC/TbarChannel",
#"MC/TWChannel",
#"MC/TbarWChannel",
#"MC/W1Jet",
#"MC/W2Jets",
#"MC/W3Jets",
#"MC/W4Jets",
#"MC/TTBarSemiLep",
#"MC/TTBarFullLep",
#"MC/ZJets",
#"MC/QCDMu",
#"MC/QCDMuBig",
#"MC/WW",
#"MC/WZ",
#"MC/ZZ",
#"MC/QCD_Pt_80to170_BCtoE",
#"MC/QCD_Pt_170to250_BCtoE",
#"MC/QCD_Pt_20to30_BCtoE",
#"MC/QCD_Pt_30to80_BCtoE",
#"MC/QCD_Pt_80to170_EMEnriched",
#"MC/QCD_Pt_20to30_EMEnriched",
#"MC/QCD_Pt_30to80_EMEnriched",
#"MC/QCD_Pt_170to250_EMEnriched",
#
#"Data/Mu_A_22Jan",
#"Data/Mu_B1_22Jan",
#"Data/Mu_B2_22Jan",
#"Data/Mu_B3_22Jan",
#"Data/Mu_B4_22Jan",
#"Data/Mu_C1_22Jan",
#"Data/Mu_C2_22Jan",
#"Data/Mu_C3_22Jan",
#"Data/Mu_C4_22Jan",
#"Data/Mu_C5_22Jan",
#"Data/Mu_C6_22Jan",
#"Data/Mu_D1_22Jan",
#"Data/Mu_D2_22Jan",
#"Data/Mu_D3_22Jan",
#"Data/Mu_D4_22Jan",
#"Data/Mu_D5_22Jan",
#"Data/Mu_D6_22Jan",
#
#"MC/Systematics/TWChannelFullLep_Q2Up",
#"MC/Systematics/TWChannelFullLep_Q2Down",
#"MC/Systematics/TbarWChannelFullLep_Q2Up",
#"MC/Systematics/TbarWChannelFullLep_Q2Down",
#
#"MC/Systematics/TWChannelThadWlep_Q2Up",
#"MC/Systematics/TWChannelThadWlep_Q2Down",
#"MC/Systematics/TWChannelTlepWhad_Q2Up",
#"MC/Systematics/TWChannelTlepWhad_Q2Down",
#"MC/Systematics/TbarWChannelThadWlep_Q2Up",
#"MC/Systematics/TbarWChannelThadWlep_Q2Down",
#"MC/Systematics/TbarWChannelTlepWhad_Q2Up",
#"MC/Systematics/TbarWChannelTlepWhad_Q2Down",
#
#"Ferdos/Systematics/TTBar_MassUp",
#"Ferdos/Systematics/TTBar_MassDown",
#"Ferdos/Systematics/TTBar_MatchingUp",
#"Ferdos/Systematics/TTBar_MatchingDown",
#"Ferdos/Systematics/TTBar_Q2Up",
##"Ferdos/Systematics/TTBar_Q2Down",
#
#"Ferdos/Systematics/TChannel_MassUp",
#"Ferdos/Systematics/TChannel_MassDown",
#"Ferdos/Systematics/TbarChannel_MassUp",
#"Ferdos/Systematics/TbarChannel_MassDown",
#"Ferdos/Systematics/TChannel_Q2Up",
#"Ferdos/Systematics/TChannel_Q2Down",
#"Ferdos/Systematics/TbarChannel_Q2Up",
#"Ferdos/Systematics/TbarChannel_Q2Down#",
#
#
#"Ferdos/Systematics/SChannel_MassUp",
#"Ferdos/Systematics/SChannel_MassDown",
#"Ferdos/Systematics/SbarChannel_MassUp",
#"Ferdos/Systematics/SbarChannel_MassDown",
#"Ferdos/SbarChannel_MassUp",
#"Ferdos/SbarChannel_MassDown",
#"Ferdos/Systematics/SChannel_Q2Up",
#"Ferdos/Systematics/SChannel_Q2Down",
#"Ferdos/Systematics/SbarChannel_Q2Up",
#"Ferdos/Systematics/SbarChannel_Q2Down",
#
#
#"crab_Full2701_DEG_05Oct/160127_115152/0000",
#"crab_Full2701_DEG_PR/160127_115307/0000",
#"crab_Full2701_TTv1/160127_114324/0000",
#"crab_Full2701_TTv2/160127_114346/0000",
#"crab_Full2701_TTv2/160127_114346/0001",
#"crab_Full2701_TTv2/160127_114346/0002",
#
#
#"Ferdos/Systematics/WJets_Q2Up",
#"Ferdos/Systematics/WJets_Q2Down",
#"Ferdos/Systematics/WJets_MatchingUp",
#"Ferdos/Systematics/WJets_MatchingDown",
#"0000",
#"0001",
#"0002", 
#"ttHTobb_M125_13TeV_powheg_pythia8/crab_Partial0812_TTHbb/151208_142607/0000", 
#"TT_TuneCUETP8M1_13TeV-powheg-pythia8/crab_Partial0812_TTv2/151208_142642/0000",
#"TT_TuneCUETP8M1_13TeV-powheg-pythia8/crab_Partial0812_TTv2/151208_142642/0001",
#"TT_TuneCUETP8M1_13TeV-powheg-pythia8/crab_Partial0812_TTv2/151208_142642/0002",
#"TT_TuneCUETP8M1_13TeV-powheg-pythia8/crab_Partial0812_TTv1/151208_142624/0000",
  ]

nstart = 0
nRetries = 1 
nRetries = nRetries+1 #actually the range is +1

for attempt in range(1,nRetries):
    nstart = 0;
    print "attempt number " +str(attempt)
    for channel in channels:

	localFiledir = localdir + "/"+channel
	if not os.path.isdir(localFiledir):
	    commands.getoutput("mkdir -p "+localFiledir)
        command_ls_channel = command_ls[:-1] + "/"+ channel 
        command_ls_channel_output = commands.getoutput(command_ls_channel)
        files = command_ls_channel_output.split('\n')
        
        print command_ls_channel
        print command_ls_channel_output
        print "channel " + channel + " files: "
        size = len(files)
	print size
	#if size==0:
	#    continue
        for iline in files[1:]:
	    iline_tuple = iline.split()
	    file = iline_tuple[8]
	    print file
            gridfile_size = iline_tuple[4]
            if "Merged" in channel:
                #filename = str(re.findall(path + "*\.root" , file))[(int(len(path))+2):-2]
                filename = str(re.findall(".*\.root" , file))[(int(len(path))+int(len(channel))+3):-2]
            else: filename = str(re.findall("OutTree.*\.root" , file))[2:-2]
	    print filename
            if filename=='':
		continue;
            command_check_doubles = 'ls -ltr ' + localFiledir +"/"+ filename
            doubles_check = commands.getstatusoutput(command_check_doubles)
	    #print filename
	    #print command_check_doubles
	    #print doubles_check
            #if nstart==size: print " I think I'm the last one " + filename 
	    equalSize = False
            if doubles_check[0]==0:
		#print doubles_check
		local_file_size = doubles_check[1].split()[4]
		#print local_file_size, gridfile_size
		if local_file_size==gridfile_size:
		    equalSize = True
		else:
		    print "not equal size, remove local one and redownload from grid"
		    command_rm = 'rm -f ' + localFiledir + "/"+filename
		    print command_rm
		    os.system(command_rm)
            if doubles_check[0]==512 or (not equalSize):
                nstart = nstart +1
                #print " is multiple " +str(nstart)
                command_cp = ""
                if (nstart % nSimultaneous ==0):
	 	   # command_cp = 'lcg-cp -b -D srmv2 -T srmv2 "srm://'+se + ":" + port + "/srm/managerv2?SFN=" + file + '" "'+ localFiledir + "/" + filename +'"'
                   #command_cp = 'cp ' + path + "/" + channel + "/" + file + " " + localFiledir 
                    command_cp = 'srmcp --debug srm://srm.ihep.ac.cn'+ path + "/" + channel + "/" + filename + " " + 'file:///' + localFiledir 
                    print "test once "+str(nstart)
                else:
                    #command_cp = 'lcg-cp -b -D srmv2 -T srmv2 "srm://'+se + ":" + port + "/srm/managerv2?SFN=" + file + '" "'+ localFiledir + "/" + filename +'" &'
                    #command_cp = 'cp ' + path + "/" + channel + "/" + file + " " + localFiledir  
                    command_cp = 'srmcp --debug srm://srm.ihep.ac.cn'+ path + "/" + channel + "/" + filename + " " + 'file:///' + localFiledir 
                #If it's not found try to get it
                print "check  "+str(doubles_check)
		if int(gridfile_size)!=0:
                    os.system(command_cp)
                print command_cp
	    print gridfile_size
	    if int(gridfile_size)==0:  ##### remove 0 size file to avoid problem....
		command_rm = 'rm -f ' + localFiledir + "/"+filename
		print "remove 0 size file ", command_rm
		os.system(command_rm)
            #print file
            #print filename

    

#print dir
#print command_ls
