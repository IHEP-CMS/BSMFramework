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
 datasetnames  = ['TTH', 'TTv1']
 datasetinputs = ['/ttHToNonbb_M125_13TeV_powheg_pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/MINIAODSIM', '/TT_TuneCUETP8M1_13TeV-powheg-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/MINIAODSIM']
 #numdatasets   = len(len(datasetnames))
 #for d in range(0,numdatasets):
 for d in range(0,2):
  config.section_('General')
  config.General.requestName = 'datasetnames[d]'
  config.General.workArea    = 'datasetnames[d]'
  config.section_('JobType')
  config.JobType.pluginName = 'Analysis'
  config.JobType.psetName   = '/afs/cern.ch/work/f/fromeo/CMSSW_7_4_12_patch4/src/BSMFramework/BSM3G_TNT_Maker/python/miniAOD_TTHLep_SR.py'
  config.section_('Data')
  config.Data.inputDataset  = 'datasetinputs[d]'
  config.Data.inputDBS      = 'global'
  config.Data.splitting     = 'FileBased'
  config.Data.totalUnits    = 2500 #With 'FileBased' splitting tells how many files to analyse
  config.Data.unitsPerJob   = 1    #Tells you how many units per job
  config.Data.outLFNDirBase = '/store/user/fromeo/'
  config.section_('Site')
  config.Site.storageSite = 'T2_CN_Beijing'
  submit(config)
