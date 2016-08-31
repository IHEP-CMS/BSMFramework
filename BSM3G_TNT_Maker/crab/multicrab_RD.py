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
'Firts2016_2_SEleB',
'Firts2016_2_SEleC',
'First2016_2_SMuB',
'First2016_2_SMuC'
                 ]
 datasetinputs = [
'/SingleElectron/Run2016B-PromptReco-v2/MINIAOD',
'/SingleElectron/Run2016C-PromptReco-v2/MINIAOD',
'/SingleMuon/Run2016B-PromptReco-v2/MINIAOD',
'/SingleMuon/Run2016C-PromptReco-v2/MINIAOD'
                 ]
 for d in range(0,len(datasetnames)):
  config.section_('General')
  config.General.requestName = datasetnames[d]
  config.General.workArea    = datasetnames[d]
  config.section_('JobType')
  config.JobType.pluginName  = 'Analysis'
  config.JobType.psetName    = '/afs/cern.ch/work/f/fromeo/CMSSW_8_0_12_FW/src/BSMFramework/BSM3G_TNT_Maker/python/miniAOD_RD.py'
  config.section_('Data')
  config.Data.inputDataset   = datasetinputs[d]
  config.Data.inputDBS       = 'global'
  config.Data.splitting      = 'LumiBased'
  config.Data.unitsPerJob    = 30
  #Golden
  config.Data.lumiMask       = '/afs/cern.ch/work/f/fromeo/CMSSW_8_0_12_FW/src/BSMFramework/BSM3G_TNT_Maker/data/JSON/Cert_271036-276097_13TeV_PromptReco_Collisions16_JSON_NoL1T_v2.txt'
  config.Data.outLFNDirBase  = '/store/user/fromeo/'
  config.section_('Site')
  config.Site.storageSite    = 'T2_CN_Beijing'
  submit(config)
