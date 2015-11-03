from WMCore.Configuration import Configuration
config = Configuration()

config.section_('General')
config.General.requestName = 'DYJetsToLL_HT100To200'
config.General.workArea    = 'DYJetsToLL_HT100To200'

config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.psetName   = 'miniAOD_MC.py'
config.section_('Data')
#config.Data.inputDataset  = '/TprimeBToTZ_M-700_LH_TuneCUETP8M1_13TeV-madgraph-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/MINIAODSIM'
#config.Data.inputDataset  = '/TprimeBToTZ_M-1000_LH_TuneCUETP8M1_13TeV-madgraph-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/MINIAODSIM'
#config.Data.inputDataset  = '/TprimeBToTZ_M-1400_LH_TuneCUETP8M1_13TeV-madgraph-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/MINIAODSIM'
#config.Data.inputDataset  = '/TprimeBToTZ_M-1700_LH_TuneCUETP8M1_13TeV-madgraph-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/MINIAODSIM'
config.Data.inputDataset  = '/DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2/MINIAODSIM'
#config.Data.inputDataset  = '/DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2/MINIAODSIM'
#config.Data.inputDataset  = '/DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2/MINIAODSIM'
#config.Data.inputDataset  = '/DYJetsToLL_M-50_HT-600toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2/MINIAODSIM'
#config.Data.inputDataset  = '/TT_TuneCUETP8M1_13TeV-powheg-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2/MINIAODSIM'
#config.Data.inputDataset  = '/ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/MINIAODSIM'
#config.Data.inputDataset  = '/ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/MINIAODSIM'
config.Data.inputDBS      = 'global'
config.Data.splitting     = 'FileBased'
config.Data.totalUnits    = 3000
config.Data.unitsPerJob   = 4
config.Data.outLFNDirBase = '/store/user/aspiezia/'

config.section_('Site')
config.Site.storageSite = 'T2_CN_Beijing'
