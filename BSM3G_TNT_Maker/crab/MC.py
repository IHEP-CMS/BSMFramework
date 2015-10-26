from WMCore.Configuration import Configuration
config = Configuration()

config.section_('General')
config.General.requestName = 'ZprimeToTauTau_M2000_MiniAODv1_4'
config.General.workArea    = 'ZprimeToTauTau_M2000_MiniAODv1_4'

config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.psetName   = '/afs/cern.ch/work/f/fromeo/CMSSW_7_4_12_patch4/src/BSMFramework/BSM3G_TNT_Maker/python/miniAOD.py'

config.section_('Data')
config.Data.inputDataset  = '/ZprimeToTauTau_M_2000_TuneCUETP8M1_tauola_13TeV_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/MINIAODSIM'
#config.Data.inputDataset  = '/ZprimeToTauTau_M_2000_TuneCUETP8M1_tauola_13TeV_pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/MINIAODSIM'
#config.Data.inputDataset  = '/SUSYGluGluToHToTauTau_M-160_TuneCUETP8M1_13TeV-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/MINIAODSIM'
#config.Data.inputDataset  = '/TT_TuneCUETP8M1_13TeV-powheg-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2/MINIAODSIM'
#config.Data.inputDataset  = '/TT_TuneCUETP8M1_13TeV-powheg-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9_ext3-v1/MINIAODSIM'
#config.Data.inputDataset  = '/ttHTobb_M125_13TeV_powheg_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/MINIAODSIM'
config.Data.inputDBS      = 'global'
config.Data.splitting     = 'FileBased'
config.Data.totalUnits    = 2500 #With 'FileBased' splitting tells how many files to analyse
config.Data.unitsPerJob   = 1    #Tells you how many units per job
config.Data.outLFNDirBase = '/store/user/fromeo/'

config.section_('Site')
config.Site.storageSite = 'T2_IT_Bari'
#'T2_IT_Bari' 'T2_IT_Pisa' 
