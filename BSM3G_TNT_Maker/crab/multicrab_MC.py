if __name__ == '__main__':
 #####
 ##   Multicrab configuration
 #####
 from CRABClient.UserUtilities import config, getUsernameFromSiteDB
 config = config()
 from CRABAPI.RawCommand import crabCommand
 from CRABClient.ClientExceptions import ClientException
 from httplib import HTTPException
 config.General.workArea = 'Crab_projects'
 
 def submit(config):
  try:
   crabCommand('submit', config = config)
  except HTTPException as hte:
   print "Failed submitting task: %s" % (hte.headers)
  except ClientException as cle:
   print "Failed submitting task: %s" % (cle)
 #####
 ##   Crab configuration
 #####
 datasetnames  = [
#'First2016_TTHbb',
#'First2016_TTHnbb',
#'First2016_TT',
'First2016_2_ST',
'First2016_2_SaT',
'First2016_2_DY',
'First2016_2_WJets',
'First2016_2_WW',
'First2016_2_WZ',
'First2016_2_ZZ'
                 ]
 datasetinputs = [
#'/ttHTobb_M125_13TeV_powheg_pythia8/RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/MINIAODSIM',
#'/ttHToNonbb_M125_13TeV_powheg_pythia8/RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/MINIAODSIM',
#'/TT_TuneCUETP8M1_13TeV-powheg-pythia8/RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14_ext3-v1/MINIAODSIM',
'/ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v2/MINIAODSIM',
'/ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/MINIAODSIM',
'/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0_ext1-v1/MINIAODSIM',
'/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/MINIAODSIM',
'/WW_TuneCUETP8M1_13TeV-pythia8/RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/MINIAODSIM',
'/WZ_TuneCUETP8M1_13TeV-pythia8/RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/MINIAODSIM',
'/ZZ_TuneCUETP8M1_13TeV-pythia8/RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/MINIAODSIM'
                 ]
 for d in range(0,len(datasetnames)):
  config.section_('General')
  config.General.requestName = datasetnames[d]
  config.General.workArea    = datasetnames[d]
  config.section_('JobType')
  config.JobType.pluginName  = 'Analysis'
  config.JobType.psetName    = '/afs/cern.ch/work/f/fromeo/CMSSW_8_0_12_FW/src/BSMFramework/BSM3G_TNT_Maker/python/miniAOD_MC.py'
  config.section_('Data')
  config.Data.inputDataset   = datasetinputs[d]
  config.Data.inputDBS       = 'global'
  config.Data.splitting      = 'FileBased'
  config.Data.totalUnits     = 2500 #With 'FileBased' splitting tells how many files to analyse
  config.Data.unitsPerJob    = 1    #Tells you how many units per job
  config.Data.outLFNDirBase  = '/store/user/fromeo/'
  config.section_('Site')
  config.Site.storageSite    = 'T2_CN_Beijing'
  submit(config)
