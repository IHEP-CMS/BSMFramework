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
'TTHbb_SEleB',
'TTHbb_SEleC',
'TTHbb_SEleD',
'TTHbb_SEleE',
'TTHbb_SEleF',
'TTHbb_SEleG',
'TTHbb_SMuB',
'TTHbb_SMuC',
'TTHbb_SMuD',
'TTHbb_SMuE',
'TTHbb_SMuF',
'TTHbb_SMuG'
                 ]
 datasetinputs = [
'/SingleElectron/Run2016B-23Sep2016-v3/MINIAOD',
'/SingleElectron/Run2016C-23Sep2016-v1/MINIAOD',
'/SingleElectron/Run2016D-23Sep2016-v1/MINIAOD',
'/SingleElectron/Run2016E-23Sep2016-v1/MINIAOD',
'/SingleElectron/Run2016F-23Sep2016-v1/MINIAOD',
'/SingleElectron/Run2016G-23Sep2016-v1/MINIAOD',
'/SingleMuon/Run2016B-23Sep2016-v3/MINIAOD',
'/SingleMuon/Run2016C-23Sep2016-v1/MINIAOD',
'/SingleMuon/Run2016D-23Sep2016-v1/MINIAOD',
'/SingleMuon/Run2016E-23Sep2016-v1/MINIAOD',
'/SingleMuon/Run2016F-23Sep2016-v1/MINIAOD',
'/SingleMuon/Run2016G-23Sep2016-v1/MINIAOD'
                 ]
 for d in range(0,len(datasetnames)):
  config.section_('General')
  config.General.requestName = datasetnames[d]
  config.General.workArea    = datasetnames[d]
  config.section_('JobType')
  config.JobType.pluginName  = 'Analysis'
  config.JobType.psetName    = '/afs/cern.ch/work/f/fromeo/CMSSW_8_0_20_FW/src/BSMFramework/BSM3G_TNT_Maker/python/miniAOD_RDbj.py'
  config.JobType.allowUndistributedCMSSW = True
  config.section_('Data')
  config.Data.inputDataset   = datasetinputs[d]
  config.Data.inputDBS       = 'global'
  config.Data.splitting      = 'LumiBased'
  config.Data.unitsPerJob    = 70
  #Golden
  config.Data.lumiMask       = '/afs/cern.ch/work/f/fromeo/CMSSW_8_0_20_FW/src/BSMFramework/BSM3G_TNT_Maker/data/JSON/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'
  #config.Data.outLFNDirBase  = '/store/user/fromeo/'
  config.Data.outLFNDirBase = '/store/group/cmst3/user/fromeo/TTHbb/'
  config.section_('Site')
  config.Site.storageSite    = 'T2_CH_CERN'#T2_CN_Beijing'
  submit(config)
