from WMCore.Configuration import Configuration
config = Configuration()

config.section_('General')
config.General.requestName = 'RunC'
config.General.workArea    = '2017v01'

config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.psetName   = 'miniAOD_RunC.py'
config.JobType.sendExternalFolder = True
config.section_('Data')
config.Data.inputDataset  = '/MET/Run2017C-31Mar2018-v1/MINIAOD'
config.Data.inputDBS      = 'global'
config.Data.splitting     = 'Automatic'
config.Data.lumiMask      = '/afs/cern.ch/work/a/aspiezia/Ntuplizer/2017/new/CMSSW_9_4_10/src/BSMFramework/BSM3G_TNT_Maker/data/JSON/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt'
config.Data.outLFNDirBase = '/store/user/aspiezia/'

config.section_('Site')
config.Site.storageSite = 'T2_CN_Beijing'

