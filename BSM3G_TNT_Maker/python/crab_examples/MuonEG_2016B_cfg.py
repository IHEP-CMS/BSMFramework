from WMCore.Configuration import Configuration
import os
config = Configuration()

config.section_("General")
config.General.requestName = "MuonEG_2016B"
config.General.workArea = "grid"
config.General.transferOutputs=True

config.section_("JobType")
config.JobType.pluginName = "Analysis"
config.JobType.psetName = "/afs/cern.ch/work/f/fromeo/CMSSW_8_0_23_BTagAna/src/RecoBTag/PerformanceMeasurements/test/runBTagAnalyzer_cfg.py"
config.JobType.disableAutomaticOutputCollection = False
config.JobType.pyCfgParams = ['useTTbarFilter=True', 'miniAOD=True', 'maxEvents=-1', 'runOnData=True']

config.section_("Data")
config.Data.inputDataset = "/MuonEG/Run2016B-23Sep2016-v3/MINIAOD"
config.Data.inputDBS = "global"
config.Data.splitting = "LumiBased"
config.Data.unitsPerJob = 100
config.Data.lumiMask = '/afs/cern.ch/work/f/fromeo/CMSSW_8_0_23_BTagAna/src/RecoBTag/PerformanceMeasurements/test/ttbar/data/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'
config.Data.publication = True
config.Data.ignoreLocality = False
config.Data.outLFNDirBase = '/store/group/phys_btag/Commissioning/TTbar2//TTbarCrabNt8023/'

config.section_("Site")
config.Site.storageSite = "T2_CH_CERN"
