if __name__ == '__main__':
 #####
 ##   Multicrab configuration
 #####
 import sys
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
'TTHbb_SEleBlockB1'
#'TTHbb_SEleBlockC1',
#'TTHbb_SEleBlockD1',
#'TTHbb_SEleBlockE1',
#'TTHbb_SEleBlockF1',
#'TTHbb_SEleBlockF2',
#'TTHbb_SEleBlockG1',
#'TTHbb_SEleBlockH1',
#'TTHbb_SEleBlockH2',
#'TTHbb_SMuBlockB1',
#'TTHbb_SMuBlockC1',
#'TTHbb_SMuBlockD1',
#'TTHbb_SMuBlockE1',
#'TTHbb_SMuBlockF1',
#'TTHbb_SMuBlockF2',
#'TTHbb_SMuBlockG1',
#'TTHbb_SMuBlockH1',
#'TTHbb_SMuBlockH2',
#'TTHbb_DblEGBlockB1',
#'TTHbb_DblEGBlockC1',
#'TTHbb_DblEGBlockD1',
#'TTHbb_DblEGBlockE1',
#'TTHbb_DblEGBlockF1',
#'TTHbb_DblEGBlockF2',
#'TTHbb_DblEGBlockG1',
#'TTHbb_DblEGBlockH1',
#'TTHbb_DblEGBlockH2',
#'TTHbb_DblMuBlockB1',
#'TTHbb_DblMuBlockC1',
#'TTHbb_DblMuBlockD1',
#'TTHbb_DblMuBlockE1',
#'TTHbb_DblMuBlockF1',
#'TTHbb_DblMuBlockF2',
#'TTHbb_DblMuBlockG1',
#'TTHbb_DblMuBlockH1',
#'TTHbb_DblMuBlockH2',
#'TTHbb_MuEGBlockB1',
#'TTHbb_MuEGBlockC1',
#'TTHbb_MuEGBlockD1',
#'TTHbb_MuEGBlockE1',
#'TTHbb_MuEGBlockF1',
#'TTHbb_MuEGBlockF2',
#'TTHbb_MuEGBlockG1',
#'TTHbb_MuEGBlockH1',
#'TTHbb_MuEGBlockH2'
                 ]
 datasetinputs = [
 # SingleElectron dataset : AT LEAST 1 high-energy electron in the event.
 '/SingleElectron/Run2016B-23Sep2016-v3/MINIAOD'
 #'/SingleElectron/Run2016C-23Sep2016-v1/MINIAOD',
 #'/SingleElectron/Run2016D-23Sep2016-v1/MINIAOD',
 #'/SingleElectron/Run2016E-23Sep2016-v1/MINIAOD',
 #'/SingleElectron/Run2016F-23Sep2016-v1/MINIAOD',
 #'/SingleElectron/Run2016F-23Sep2016-v1/MINIAOD',
 #'/SingleElectron/Run2016G-23Sep2016-v1/MINIAOD',
 #'/SingleElectron/Run2016H-PromptReco-v2/MINIAOD',
 #'/SingleElectron/Run2016H-PromptReco-v3/MINIAOD,',
 # SingleMuon dataset : AT LEAST 1 high-energy muon in the event.
 #'/SingleMuon/Run2016B-23Sep2016-v3/MINIAOD',
 #'/SingleMuon/Run2016C-23Sep2016-v1/MINIAOD',
 #'/SingleMuon/Run2016D-23Sep2016-v1/MINIAOD',
 #'/SingleMuon/Run2016E-23Sep2016-v1/MINIAOD',
 #'/SingleMuon/Run2016F-23Sep2016-v1/MINIAOD',
 #'/SingleMuon/Run2016F-23Sep2016-v1/MINIAOD',
 #'/SingleMuon/Run2016G-23Sep2016-v1/MINIAOD',
 #'/SingleMuon/Run2016H-PromptReco-v2/MINIAOD',
 #'/SingleMuon/Run2016H-PromptReco-v3/MINIAOD',
 #'/DoubleEG/Run2016B-23Sep2016-v3/MINIAOD',
 #'/DoubleEG/Run2016C-23Sep2016-v1/MINIAOD',
 #'/DoubleEG/Run2016D-23Sep2016-v1/MINIAOD',
 #'/DoubleEG/Run2016E-23Sep2016-v1/MINIAOD',
 #'/DoubleEG/Run2016F-23Sep2016-v1/MINIAOD',
 #'/DoubleEG/Run2016F-23Sep2016-v1/MINIAOD',
 #'/DoubleEG/Run2016G-23Sep2016-v1/MINIAOD',
 #'/DoubleEG/Run2016H-PromptReco-v2/MINIAOD',
 #'/DoubleEG/Run2016H-PromptReco-v3/MINIAOD',
 #'/MuonEG/Run2016B-23Sep2016-v3/MINIAOD',
 #'/MuonEG/Run2016C-23Sep2016-v1/MINIAOD',
 #'/MuonEG/Run2016D-23Sep2016-v1/MINIAOD',
 #'/MuonEG/Run2016E-23Sep2016-v1/MINIAOD',
 #'/MuonEG/Run2016F-23Sep2016-v1/MINIAOD',
 #'/MuonEG/Run2016F-23Sep2016-v1/MINIAOD',
 #'/MuonEG/Run2016G-23Sep2016-v1/MINIAOD',
 #'/MuonEG/Run2016H-PromptReco-v2/MINIAOD',
 #'/MuonEG/Run2016H-PromptReco-v3/MINIAOD',
 #'/DoubleMuon/Run2016B-23Sep2016-v3/MINIAOD',
 #'/DoubleMuon/Run2016C-23Sep2016-v1/MINIAOD',
 #'/DoubleMuon/Run2016D-23Sep2016-v1/MINIAOD',
 #'/DoubleMuon/Run2016E-23Sep2016-v1/MINIAOD',
 #'/DoubleMuon/Run2016F-23Sep2016-v1/MINIAOD',
 #'/DoubleMuon/Run2016F-23Sep2016-v1/MINIAOD',
 #'/DoubleMuon/Run2016G-23Sep2016-v1/MINIAOD',
 #'/DoubleMuon/Run2016H-PromptReco-v2/MINIAOD',
 #'/DoubleMuon/Run2016H-PromptReco-v3/MINIAOD'
                 ]

JECBlockBCD = [
'BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_L1FastJet_AK4PFchs.txt',
'BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_L2Relative_AK4PFchs.txt',
'BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_L3Absolute_AK4PFchs.txt',
'BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_L2L3Residual_AK4PFchs.txt',
'BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_Uncertainty_AK4PFchs.txt',
'BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_L1FastJet_AK4PFPuppi.txt',
'BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_L2Relative_AK4PFPuppi.txt',
'BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_L3Absolute_AK4PFPuppi.txt',
'BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_L2L3Residual_AK4PFPuppi.txt',
'BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_Uncertainty_AK4PFPuppi.txt',
'BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_L1FastJet_AK8PFchs.txt',
'BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_L2Relative_AK8PFchs.txt',
'BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_L3Absolute_AK8PFchs.txt',
'BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_L2L3Residual_AK8PFchs.txt',
'BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_Uncertainty_AK8PFchs.txt'
]

JECBlockEF = [
'BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer16_23Sep2016EFV4_DATA/Summer16_23Sep2016EFV4_DATA_L1FastJet_AK4PFchs.txt',
'BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer16_23Sep2016EFV4_DATA/Summer16_23Sep2016EFV4_DATA_L2Relative_AK4PFchs.txt',
'BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer16_23Sep2016EFV4_DATA/Summer16_23Sep2016EFV4_DATA_L3Absolute_AK4PFchs.txt',
'BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer16_23Sep2016EFV4_DATA/Summer16_23Sep2016EFV4_DATA_L2L3Residual_AK4PFchs.txt',
'BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer16_23Sep2016EFV4_DATA/Summer16_23Sep2016EFV4_DATA_Uncertainty_AK4PFchs.txt',
'BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer16_23Sep2016EFV4_DATA/Summer16_23Sep2016EFV4_DATA_L1FastJet_AK4PFPuppi.txt',
'BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer16_23Sep2016EFV4_DATA/Summer16_23Sep2016EFV4_DATA_L2Relative_AK4PFPuppi.txt',
'BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer16_23Sep2016EFV4_DATA/Summer16_23Sep2016EFV4_DATA_L3Absolute_AK4PFPuppi.txt',
'BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer16_23Sep2016EFV4_DATA/Summer16_23Sep2016EFV4_DATA_L2L3Residual_AK4PFPuppi.txt',
'BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer16_23Sep2016EFV4_DATA/Summer16_23Sep2016EFV4_DATA_Uncertainty_AK4PFPuppi.txt',
'BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer16_23Sep2016EFV4_DATA/Summer16_23Sep2016EFV4_DATA_L1FastJet_AK8PFchs.txt',
'BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer16_23Sep2016EFV4_DATA/Summer16_23Sep2016EFV4_DATA_L2Relative_AK8PFchs.txt',
'BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer16_23Sep2016EFV4_DATA/Summer16_23Sep2016EFV4_DATA_L3Absolute_AK8PFchs.txt',
'BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer16_23Sep2016EFV4_DATA/Summer16_23Sep2016EFV4_DATA_L2L3Residual_AK8PFchs.txt',
'BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer16_23Sep2016EFV4_DATA/Summer16_23Sep2016EFV4_DATA_Uncertainty_AK8PFchs.txt'
]

JECBlockG = [
'BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer16_23Sep2016GV4_DATA/Summer16_23Sep2016GV4_DATA_L1FastJet_AK4PFchs.txt',
'BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer16_23Sep2016GV4_DATA/Summer16_23Sep2016GV4_DATA_L2Relative_AK4PFchs.txt',
'BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer16_23Sep2016GV4_DATA/Summer16_23Sep2016GV4_DATA_L3Absolute_AK4PFchs.txt',
'BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer16_23Sep2016GV4_DATA/Summer16_23Sep2016GV4_DATA_L2L3Residual_AK4PFchs.txt',
'BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer16_23Sep2016GV4_DATA/Summer16_23Sep2016GV4_DATA_Uncertainty_AK4PFchs.txt',
'BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer16_23Sep2016GV4_DATA/Summer16_23Sep2016GV4_DATA_L1FastJet_AK4PFPuppi.txt',
'BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer16_23Sep2016GV4_DATA/Summer16_23Sep2016GV4_DATA_L2Relative_AK4PFPuppi.txt',
'BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer16_23Sep2016GV4_DATA/Summer16_23Sep2016GV4_DATA_L3Absolute_AK4PFPuppi.txt',
'BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer16_23Sep2016GV4_DATA/Summer16_23Sep2016GV4_DATA_L2L3Residual_AK4PFPuppi.txt',
'BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer16_23Sep2016GV4_DATA/Summer16_23Sep2016GV4_DATA_Uncertainty_AK4PFPuppi.txt',
'BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer16_23Sep2016GV4_DATA/Summer16_23Sep2016GV4_DATA_L1FastJet_AK8PFchs.txt',
'BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer16_23Sep2016GV4_DATA/Summer16_23Sep2016GV4_DATA_L2Relative_AK8PFchs.txt',
'BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer16_23Sep2016GV4_DATA/Summer16_23Sep2016GV4_DATA_L3Absolute_AK8PFchs.txt',
'BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer16_23Sep2016GV4_DATA/Summer16_23Sep2016GV4_DATA_L2L3Residual_AK8PFchs.txt',
'BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer16_23Sep2016GV4_DATA/Summer16_23Sep2016GV4_DATA_Uncertainty_AK8PFchs.txt'
]

JECBlockH = [
'BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer16_23Sep2016HV4_DATA/Summer16_23Sep2016HV4_DATA_L1FastJet_AK4PFchs.txt',
'BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer16_23Sep2016HV4_DATA/Summer16_23Sep2016HV4_DATA_L2Relative_AK4PFchs.txt',
'BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer16_23Sep2016HV4_DATA/Summer16_23Sep2016HV4_DATA_L3Absolute_AK4PFchs.txt',
'BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer16_23Sep2016HV4_DATA/Summer16_23Sep2016HV4_DATA_L2L3Residual_AK4PFchs.txt',
'BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer16_23Sep2016HV4_DATA/Summer16_23Sep2016HV4_DATA_Uncertainty_AK4PFchs.txt',
'BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer16_23Sep2016HV4_DATA/Summer16_23Sep2016HV4_DATA_L1FastJet_AK4PFPuppi.txt',
'BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer16_23Sep2016HV4_DATA/Summer16_23Sep2016HV4_DATA_L2Relative_AK4PFPuppi.txt',
'BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer16_23Sep2016HV4_DATA/Summer16_23Sep2016HV4_DATA_L3Absolute_AK4PFPuppi.txt',
'BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer16_23Sep2016HV4_DATA/Summer16_23Sep2016HV4_DATA_L2L3Residual_AK4PFPuppi.txt',
'BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer16_23Sep2016HV4_DATA/Summer16_23Sep2016HV4_DATA_Uncertainty_AK4PFPuppi.txt',
'BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer16_23Sep2016HV4_DATA/Summer16_23Sep2016HV4_DATA_L1FastJet_AK8PFchs.txt',
'BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer16_23Sep2016HV4_DATA/Summer16_23Sep2016HV4_DATA_L2Relative_AK8PFchs.txt',
'BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer16_23Sep2016HV4_DATA/Summer16_23Sep2016HV4_DATA_L3Absolute_AK8PFchs.txt',
'BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer16_23Sep2016HV4_DATA/Summer16_23Sep2016HV4_DATA_L2L3Residual_AK8PFchs.txt',
'BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer16_23Sep2016HV4_DATA/Summer16_23Sep2016HV4_DATA_Uncertainty_AK8PFchs.txt'
]

goodRunsLists = [
'/afs/cern.ch/work/j/jthomasw/private/IHEP/CMSSW/CMSSW_8_0_26_patch2/src/BSMFramework/BSM3G_TNT_Maker/data/JSON/runF_split_JSON/Cert_276831-278801_13TeV_23Sep2016ReReco_Collisions16_JSON.txt',
'/afs/cern.ch/work/j/jthomasw/private/IHEP/CMSSW/CMSSW_8_0_26_patch2/src/BSMFramework/BSM3G_TNT_Maker/data/JSON/runF_split_JSON/Cert_278802-280385_13TeV_23Sep2016ReReco_Collisions16_JSON.txt',
'/afs/cern.ch/work/j/jthomasw/private/IHEP/CMSSW/CMSSW_8_0_26_patch2/src/BSMFramework/BSM3G_TNT_Maker/data/JSON/GOLDEN/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'
]

#for d in range(0,len(datasetnames)):
for d in range(0,1):
    print 'multicrab.py: Running datasetname: ', datasetnames[d]
    JECFiles = []
    tempJSON = ''
    if 'BlockB' in datasetnames[d] or 'BlockC' in datasetnames[d] or 'BlockD' in datasetnames[d]:
        print 'multicrab.py: Run Block B/C/D'
        JECFiles = JECBlockBCD
        tempJSON = goodRunsLists[2]
    if 'BlockE' in datasetnames[d]:
        print 'multicrab.py: Run Block E'
        JECFiles = JECBlockEF
        tempJSON = goodRunsLists[2]
    if 'BlockG' in datasetnames[d]:
        print 'multicrab.py: Run Block G'
        JECFiles = JECBlockG
        tempJSON = goodRunsLists[2]
    if 'BlockH' in datasetnames[d]:
        print 'multicrab.py: Run Block H'
        JECFiles = JECBlockH
        tempJSON = goodRunsLists[2]
    if 'BlockF1' in datasetnames[d]:
        print 'multicrab.py: Run Block F1'
        JECFiles = JECBlockEF
        tempJSON = goodRunsLists[0]
    if 'BlockF2' in datasetnames[d]:
        print 'multicrab.py: Run Block F2'
        JECFiles = JECBlockEF
        tempJSON = goodRunsLists[1]

    print 'multicrab.py: JSON File = ', tempJSON
    try:
        testJECFiles = JECFiles[14]
    except(IndexError):
        print 'multicrab.py: Failed to access JEC list element.'
        print 'multicrab.py: Not eneough JEC files proivided.'
        sys.exit()
    try:
        testJSON = goodRunsLists[2]
    except(IndexError):
        print 'multicrab.py: Failed to access JSON list element.'
        print 'multicrab.py: Not eneough JSON files proivided.'
        sys.exit()

    nameJECAK4PFchsDATA1 = "optionJECAK4PFchsDATA1="+JECFiles[0]
    nameJECAK4PFchsDATA2 = "optionJECAK4PFchsDATA2="+JECFiles[1]
    nameJECAK4PFchsDATA3 = "optionJECAK4PFchsDATA3="+JECFiles[2]
    nameJECAK4PFchsDATA4 = "optionJECAK4PFchsDATA4="+JECFiles[3]
    nameJECAK4PFchsDATAUnc = "optionJECAK4PFchsDATAUnc="+JECFiles[4]
    nameJECAK4PFPuppiDATA1 = "optionJECAK4PFPuppiDATA1="+JECFiles[5]
    nameJECAK4PFPuppiDATA2 = "optionJECAK4PFPuppiDATA2="+JECFiles[6]
    nameJECAK4PFPuppiDATA3 = "optionJECAK4PFPuppiDATA3="+JECFiles[7]
    nameJECAK4PFPuppiDATA4 = "optionJECAK4PFPuppiDATA4="+JECFiles[8]
    nameJECAK4PFPuppiDATAUnc = "optionJECAK4PFPuppiDATAUnc="+JECFiles[9]
    nameJECAK8PFchsDATA1 = "optionJECAK8PFchsDATA1="+JECFiles[10]
    nameJECAK8PFchsDATA2 = "optionJECAK8PFchsDATA2="+JECFiles[11]
    nameJECAK8PFchsDATA3 = "optionJECAK8PFchsDATA3="+JECFiles[12]
    nameJECAK8PFchsDATA4 = "optionJECAK8PFchsDATA4="+JECFiles[13]
    nameJECAK8PFchsDATAUnc = "optionJECAK8PFchsDATAUnc="+JECFiles[14]

    config.section_('General')
    config.General.requestName = datasetnames[d]
    config.General.workArea    = datasetnames[d]

    config.section_('JobType')
    config.JobType.pluginName  = 'Analysis'
    # List of parameters to pass to CMSSW parameter-set configuration file:
    config.JobType.psetName    = '/afs/cern.ch/work/j/jthomasw/private/IHEP/CMSSW/CMSSW_8_0_26_patch2/src/BSMFramework/BSM3G_TNT_Maker/python/miniAOD_RD_new.py'
    config.JobType.allowUndistributedCMSSW = True
    config.JobType.sendExternalFolder = True
    ofParam = 'ofName=' + datasetnames[d]
    config.JobType.pyCfgParams = [nameJECAK4PFchsDATA1,
                                  nameJECAK4PFchsDATA2,
                                  nameJECAK4PFchsDATA3,
                                  nameJECAK4PFchsDATA4,
                                  nameJECAK4PFchsDATAUnc,
                                  nameJECAK4PFPuppiDATA1,
                                  nameJECAK4PFPuppiDATA2,
                                  nameJECAK4PFPuppiDATA3,
                                  nameJECAK4PFPuppiDATA4,
                                  nameJECAK4PFPuppiDATAUnc,
                                  nameJECAK8PFchsDATA1,
                                  nameJECAK8PFchsDATA2,
                                  nameJECAK8PFchsDATA3,
                                  nameJECAK8PFchsDATA4,
                                  nameJECAK8PFchsDATAUnc,
                                  ofParam
                                  ]

    config.section_('Data')
    config.Data.inputDataset   = datasetinputs[d]
    config.Data.inputDBS       = 'global'
    config.Data.splitting      = 'LumiBased'
    config.Data.unitsPerJob    = 30
    # Golden
    config.Data.lumiMask       = tempJSON
    config.Data.outLFNDirBase = '/store/user/jthomasw/TTHbb/BSMFramework/output/'
    print 'multicrab.py: outLFNDirBase = /store/user/jthomasw/TTHbb/BSMFramework/output/'
    #config.Data.publication = True

    config.section_('Site')
    config.Site.storageSite    = 'T2_CN_Beijing'#'T2_CH_CERN'
    print 'multicrab.py: Submitting Jobs'
    submit(config)
