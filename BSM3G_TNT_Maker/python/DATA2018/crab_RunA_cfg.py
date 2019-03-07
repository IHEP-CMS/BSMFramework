from WMCore.Configuration import Configuration
config = Configuration()

config.section_('General')
config.General.requestName = 'RunA'
config.General.workArea    = '2018v01'

config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.psetName   = 'miniAOD_RunAll.py'
config.JobType.sendExternalFolder = True
config.section_('Data')
config.Data.inputDataset  = '/MET/Run2018A-17Sep2018-v1/MINIAOD'
config.Data.inputDBS      = 'global'
config.Data.splitting     = 'Automatic'
config.Data.lumiMask      = '/afs/cern.ch/work/a/aspiezia/Ntuplizer/2018/CMSSW_10_2_5/src/BSMFramework/BSM3G_TNT_Maker/data/JSON/Cert_314472-325175_13TeV_PromptReco_Collisions18_JSON.txt'
config.Data.outLFNDirBase = '/store/user/aspiezia/'

config.section_('Site')
config.Site.storageSite = 'T2_CN_Beijing'

