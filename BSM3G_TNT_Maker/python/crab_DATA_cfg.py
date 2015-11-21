from WMCore.Configuration import Configuration
config = Configuration()

config.section_('General')
config.General.requestName = 'SingleElectronRunDv4'
config.General.workArea    = 'SingleElectronRunDv4'

config.section_('JobType')
config.JobType.pluginName = 'Analysis'
#config.JobType.psetName   = 'miniAOD_RunC.py'
config.JobType.psetName   = 'miniAOD_RunD.py'
config.section_('Data')
#config.Data.inputDataset  = '/SingleElectron/Run2015C-PromptReco-v1/MINIAOD'
#config.Data.inputDataset  = '/SingleMuon/Run2015C-PromptReco-v1/MINIAOD'
#config.Data.inputDataset  = '/SingleElectron/Run2015D-PromptReco-v3/MINIAOD'
#config.Data.inputDataset  = '/SingleMuon/Run2015D-PromptReco-v3/MINIAOD'
config.Data.inputDataset  = '/SingleElectron/Run2015D-PromptReco-v4/MINIAOD'
#config.Data.inputDataset  = '/SingleMuon/Run2015D-PromptReco-v4/MINIAOD'
config.Data.inputDBS      = 'global'
config.Data.splitting     = 'LumiBased'
config.Data.unitsPerJob   = 40
config.Data.lumiMask      = '/afs/cern.ch/work/a/aspiezia/BSM/Francesco/v2/CMSSW_7_4_12/src/BSMFramework/BSM3G_TNT_Maker/data/JSON/Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON.txt'
config.Data.outLFNDirBase = '/store/user/aspiezia/'

config.section_('Site')
config.Site.storageSite = 'T2_CN_Beijing'
