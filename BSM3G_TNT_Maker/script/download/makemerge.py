import os
import re
import sys
import commands
import subprocess

#Path you run this script
workpath = "/publicfs/cms/user/libh/jobs"


#Create the MergeFiles
##MC##
MergeScriptMCName = "Mergemc.sh"
MergeScriptMC      = file(MergeScriptMCName,"w")
##data##
MergeScriptDataName = "Mergedata.sh"
MergeScriptData      = file(MergeScriptDataName,"w")

#For Tier 2 Data
PathT2 = '/pnfs/ihep.ac.cn/data/cms/store/user/fromeo'

#For Tier 3 Directory
mclocaldir = "/publicfs/cms/data/TopQuark/cms13TeV/Samples2607/mc"
datalocaldir = "/publicfs/cms/data/TopQuark/cms13TeV/Samples2607/data"
mergelocaldir = "/publicfs/cms/data/TopQuark/cms13TeV/Samples2607/merged"

mcchannels = [

##MC##
'ttHTobb_M125_13TeV_powheg_pythia8/crab_First2016_TTHbb/160726_140559/0000',
'TT_TuneCUETP8M1_13TeV-powheg-pythia8/crab_First2016_TT/160726_140659/0000',
'TT_TuneCUETP8M1_13TeV-powheg-pythia8/crab_First2016_TT/160726_140659/0001',
'TT_TuneCUETP8M1_13TeV-powheg-pythia8/crab_First2016_TT/160726_140659/0002',
'ttHToNonbb_M125_13TeV_powheg_pythia8/crab_First2016_TTHnbb/160726_140630/0000',
#'WZ_TuneCUETP8M1_13TeV-pythia8/crab_First2016_WZ/160726_141004/0000',
'DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_First2016_2_DY/160726_151834/0000',
'DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_First2016_2_DY/160726_151834/0001',
'ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/crab_First2016_2_ST/160726_151730/0000',
'ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/crab_First2016_2_SaT/160726_151805/0000',
'WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_First2016_2_WJets/160726_151923/0000',
'ZZ_TuneCUETP8M1_13TeV-pythia8/crab_First2016_2_ZZ/160726_152149/0000',
'WW_TuneCUETP8M1_13TeV-pythia8/crab_First2016_2_WW/160726_151952/0000',
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
#
#"ZprimeToTauTau_M_1500_TuneCUETP8M1_tauola_13TeV_pythia8/crab_HN2Tau3101_Zp1500/160131_162759/0000",
#"ZprimeToTauTau_M_2000_TuneCUETP8M1_tauola_13TeV_pythia8/crab_HN2Tau3101_Zp2000/160131_162817/0000",
#"QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_Full2202_QCD_HT1000to1500/160222_224313/0000",

]

datachannels = [

## Data##
'SingleMuon/crab_First2016_SMuB/160726_142254/0000',
'SingleMuon/crab_First2016_SMuB/160726_142254/0001',
'SingleMuon/crab_First2016_SMuC/160726_142334/0000',
'SingleElectron/crab_Firts2016_SEleB/160726_142148/0000',
'SingleElectron/crab_Firts2016_SEleB/160726_142148/0001',
'SingleElectron/crab_Firts2016_SEleC/160726_142218/0000',
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
#
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
#"ExtendedWeakIsospinModel_mumujj_L15000_M1500_CalcHEP/crab_Full2202_mumujj_L15000_M1500/160222_225924/0000",
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






]

for mcchannel in mcchannels:

        mclocalFiledir = mclocaldir + "/"+mcchannel
        print >> MergeScriptMC,"hadd -f "+mclocalFiledir+"/Merged.root "+mclocalFiledir+"/*.root"

for datachannel in datachannels:

        datalocalFiledir = datalocaldir + "/"+datachannel
        print >> MergeScriptData, "hadd -f "+datalocalFiledir+"/Merged.root "+datalocalFiledir+"/*.root"




