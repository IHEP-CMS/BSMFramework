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
 datasetnames  = ['TTHtest', 'TTtest']
 datasetinputs = ['/ttHTobb_M125_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/RunIISpring16MiniAODv2-premix_withHLT_80X_mcRun2_asymptotic_v14-v1/MINIAODSIM', '/TT_TuneCUETP8M1_13TeV-powheg-pythia8/RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14_ext3-v1/MINIAODSIM']

 #numdatasets   = len(len(datasetnames))
 #for d in range(0,numdatasets):
 for d in range(0,2):
  config.section_('General')
  config.General.requestName = 'datasetnames[d]'
  config.General.workArea    = 'datasetnames[d]'
  config.section_('JobType')
  config.JobType.pluginName = 'Analysis'
  config.JobType.psetName   = '/afs/cern.ch/work/j/jthomasw/CMSSW_8_0_24_patch1/src/BSMFramework/BSM3G_TNT_Maker/python/miniAOD_MC.py'
  config.section_('Data')
  config.Data.inputDataset  = 'datasetinputs[d]'
  config.Data.inputDBS      = 'global'
  config.Data.splitting     = 'FileBased'
  config.Data.totalUnits    = 2500 #With 'FileBased' splitting tells how many files to analyse
  config.Data.unitsPerJob   = 1    #Tells you how many units per job
  config.Data.outLFNDirBase = '/store/user/jthomasw/'

  config.section_('Site')
  # Where the output files will be transmitted. Check one has permission to write here.
  config.Site.storageSite = 'T2_CN_Beijing'
submit(config)
