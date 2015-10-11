from WMCore.Configuration import Configuration
config = Configuration()

config.section_('General')
config.General.requestName = 'TTHbb_TTv2_6jets3b_coll16_Var_10_10'
config.General.workArea    = 'TTHbb_TTv2_6jets3b_coll16_Var_10_10'

config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.psetName   = '/afs/cern.ch/work/f/fromeo/CMSSW_7_4_7/src/BSMFramework/BSM3G_TNT_Maker/python/miniAOD_TTHbb.py'
config.section_('Data')
#config.Data.inputDataset  = '/SUSYGluGluToHToTauTau_M-160_TuneCUETP8M1_13TeV-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/MINIAODSIM'
config.Data.inputDataset  = '/TT_TuneCUETP8M1_13TeV-powheg-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2/MINIAODSIM'
#config.Data.inputDataset  = '/TT_TuneCUETP8M1_13TeV-powheg-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9_ext3-v1/MINIAODSIM'
#config.Data.inputDataset  = '/ttHTobb_M125_13TeV_powheg_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/MINIAODSIM'
config.Data.inputDBS      = 'global'
config.Data.splitting     = 'FileBased'
config.Data.totalUnits    = 2400
config.Data.unitsPerJob   = 6
config.Data.outLFNDirBase = '/store/user/fromeo/'

config.section_('Site')
config.Site.storageSite = 'T2_IT_Bari'
#'T2_IT_Bari' 
#'T2_IT_Pisa' 
