from WMCore.Configuration import Configuration
config = Configuration()

config.section_('General')
config.General.requestName = '2016RunC'
config.General.workArea    = '2016v01'

config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.psetName   = 'miniAOD_RunBCD.py'
config.JobType.sendExternalFolder = True
config.section_('Data')
config.Data.inputDataset  = '/MET/Run2016C-17Jul2018-v1/MINIAOD'
config.Data.inputDBS      = 'global'
config.Data.splitting     = 'Automatic'
config.Data.lumiMask      = '/afs/cern.ch/work/a/aspiezia/Ntuplizer/2017/new/CMSSW_9_4_10/src/BSMFramework/BSM3G_TNT_Maker/data/JSON/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'
config.Data.outLFNDirBase = '/store/user/aspiezia/'

config.section_('Site')
config.Site.storageSite = 'T2_CN_Beijing'

