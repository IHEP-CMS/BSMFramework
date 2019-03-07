import FWCore.ParameterSet.Config as cms

############   Initial standard configs   ############
process = cms.Process("Demo")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10000
process.load('Configuration.Geometry.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.Services_cff')
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Geometry.CaloEventSetup.CaloTowerConstituents_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
process.GlobalTag.globaltag = '102X_dataRun2_Sep2018Rereco_v1'
process.prefer("GlobalTag")
######################################################


##################   Input files   ###################
process.source = cms.Source("PoolSource",
  fileNames = cms.untracked.vstring(
    '/store/data/Run2018A/MET/MINIAOD/17Sep2018-v1/60000/F83780CE-9F43-AC48-AAE8-EB567EA13D82.root'
  ),
  skipEvents = cms.untracked.uint32(0)
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )
######################################################


##################   Output file   ###################
process.TFileService = cms.Service("TFileService",
  fileName = cms.string("OutTree.root")
)
######################################################


###############   ELECTRON ID SECTION   ##############
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV70_cff',
                 'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Fall17_94X_V2_cff']
from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
setupEgammaPostRecoSeq(process,
                       runEnergyCorrections=False,
                       era='2018-Prompt')
######################################################


#################   MET FILTERS   ####################
baddetEcallist = cms.vuint32(
    [872439604,872422825,872420274,872423218,872423215,872416066,872435036,872439336,872420273,872436907,872420147,872439731,
     872436657,872420397,872439732,872439339,872439603,872422436,872439861,872437051,872437052,872420649,872422436,872421950,
     872437185,872422564,872421566,872421695,872421955,872421567,872437184,872421951,872421694,872437056,872437057,872437313])
process.ecalBadCalibReducedMINIAODFilter = cms.EDFilter(
    "EcalBadCalibFilter",
    EcalRecHitSource = cms.InputTag("reducedEgamma:reducedEERecHits"),
    ecalMinEt        = cms.double(50.),
    baddetEcal    = baddetEcallist, 
    taggingMode = cms.bool(True),
    debug = cms.bool(False)
    )
######################################################


####################  TAU v2   #######################
from BSMFramework.BSM3G_TNT_Maker.runTauIdMVA import *
na = TauIDEmbedder(process, cms,
        debug=True,
        toKeep = ["2017v2", "newDM2017v2", "dR0p32017v2"]
        )
na.runTauID()
######################################################


###############   Analysis parameters   ##############
process.TNT = cms.EDAnalyzer("BSM3G_TNT_Maker",
  ifevtriggers      = cms.bool(True), # True means you want to require the triggers
  maxtriggerversion = cms.double(50), # please leave it as a double
  evtriggers        = cms.vstring(
    'HLT_PFMET120_PFMHT120_IDTight_v',
    'HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v',
  ),
  # Choose which information you want to use
  fillBoostedJetinfo    = cms.bool(True),
  fillelectronpatinfo   = cms.bool(True),
  filleventinfo         = cms.bool(True),
  fillgeninfo           = cms.bool(False),
  filljetinfo           = cms.bool(True),
  fillMETinfo           = cms.bool(True),
  fillmuoninfo          = cms.bool(True),
  fillPVinfo            = cms.bool(True),
  fillPileupReweight    = cms.bool(True),
  filltauinfo           = cms.bool(True),
  filltauboostedinfo    = cms.bool(False),
  filltriggerinfo       = cms.bool(True),
  fillphotoninfo        = cms.bool(False), 
  # Choose format 
  reduced   = cms.bool(True),
  is_data   = cms.bool(True),
  debug_    = cms.bool(False),
  PuppiVar  = cms.bool(False),
  MC2016    = cms.bool(False),
  # Input tags 
  bits                = cms.InputTag("TriggerResults","","HLT"),
  prescales           = cms.InputTag("patTrigger"),
  objects             = cms.InputTag("selectedPatTrigger"),  
  vertices            = cms.InputTag("offlineSlimmedPrimaryVertices"),
  beamSpot            = cms.InputTag("offlineBeamSpot"),
  muons               = cms.InputTag("slimmedMuons"),
  patElectrons        = cms.InputTag("slimmedElectrons"),
  ebRecHits = cms.InputTag("reducedEgamma","reducedEBRecHits"),
  taus                = cms.InputTag("NewTauIDsEmbedded"),
  tausboosted         = cms.InputTag("slimmedTausBoosted"),
  jets                = cms.InputTag("slimmedJets"),
  jetsPUPPI           = cms.InputTag("slimmedJetsPuppi"),
  fatjets             = cms.InputTag("slimmedJetsAK8"),
  mets                = cms.InputTag("slimmedMETs"),
  metsPUPPI           = cms.InputTag("slimmedMETsPuppi"),
  metFilterBits       = cms.InputTag("TriggerResults", "", "RECO"),
  photons             = cms.InputTag("slimmedPhotons"),
  packedPFCandidates  = cms.InputTag("packedPFCandidates"), 
  pruned              = cms.InputTag("prunedGenParticles"),
  jerAK4PFchs     =  cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JER/Fall17_V3_MC_PtResolution_AK4PFchs.txt"),
  jerAK4PFchsSF   =  cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JER/Fall17_V3_MC_SF_AK4PFchs.txt"),
  jerAK4PFPuppi   =  cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JER/Fall17_V3_MC_PtResolution_AK4PFPuppi.txt"),
  jerAK4PFPuppiSF =  cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JER/Fall17_V3_MC_SF_AK4PFPuppi.txt"),
  jerAK8PFchs     =  cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JER/Fall17_V3_MC_PtResolution_AK8PFchs.txt"),
  jerAK8PFchsSF   =  cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JER/Fall17_V3_MC_SF_AK8PFchs.txt"),
  jerAK8PFPuppi   =  cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JER/Fall17_V3_MC_PtResolution_AK8PFPuppi.txt"),
  jerAK8PFPuppiSF =  cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JER/Fall17_V3_MC_SF_AK8PFPuppi.txt"),
  jecPayloadNamesAK4PFchsMC1   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Fall17_17Nov2017_V32_MC/Fall17_17Nov2017_V32_MC_L1FastJet_AK4PFchs.txt"),
  jecPayloadNamesAK4PFchsMC2   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Fall17_17Nov2017_V32_MC/Fall17_17Nov2017_V32_MC_L2Relative_AK4PFchs.txt"),
  jecPayloadNamesAK4PFchsMC3   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Fall17_17Nov2017_V32_MC/Fall17_17Nov2017_V32_MC_L3Absolute_AK4PFchs.txt"),
  jecPayloadNamesAK4PFchsMCUnc = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Fall17_17Nov2017_V32_MC/Fall17_17Nov2017_V32_MC_Uncertainty_AK4PFchs.txt"),
  jecPayloadNamesAK4PFchsDATA1   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Fall17_17Nov2017B_V32_DATA/Fall17_17Nov2017B_V32_DATA_L1FastJet_AK4PFchs.txt"),
  jecPayloadNamesAK4PFchsDATA2   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Fall17_17Nov2017B_V32_DATA/Fall17_17Nov2017B_V32_DATA_L2Relative_AK4PFchs.txt"),
  jecPayloadNamesAK4PFchsDATA3   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Fall17_17Nov2017B_V32_DATA/Fall17_17Nov2017B_V32_DATA_L3Absolute_AK4PFchs.txt"),
  jecPayloadNamesAK4PFchsDATA4   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Fall17_17Nov2017B_V32_DATA/Fall17_17Nov2017B_V32_DATA_L2L3Residual_AK4PFchs.txt"),
  jecPayloadNamesAK4PFchsDATAUnc = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Fall17_17Nov2017B_V32_DATA/Fall17_17Nov2017B_V32_DATA_Uncertainty_AK4PFchs.txt"),
  jecPayloadNamesAK4PFPuppiMC1   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Fall17_17Nov2017_V32_MC/Fall17_17Nov2017_V32_MC_L1FastJet_AK4PFPuppi.txt"),
  jecPayloadNamesAK4PFPuppiMC2   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Fall17_17Nov2017_V32_MC/Fall17_17Nov2017_V32_MC_L2Relative_AK4PFPuppi.txt"),
  jecPayloadNamesAK4PFPuppiMC3   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Fall17_17Nov2017_V32_MC/Fall17_17Nov2017_V32_MC_L3Absolute_AK4PFPuppi.txt"),
  jecPayloadNamesAK4PFPuppiMCUnc = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Fall17_17Nov2017_V32_MC/Fall17_17Nov2017_V32_MC_Uncertainty_AK4PFPuppi.txt"),
  jecPayloadNamesAK4PFPuppiDATA1   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Fall17_17Nov2017B_V32_DATA/Fall17_17Nov2017B_V32_DATA_L1FastJet_AK4PFPuppi.txt"),
  jecPayloadNamesAK4PFPuppiDATA2   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Fall17_17Nov2017B_V32_DATA/Fall17_17Nov2017B_V32_DATA_L2Relative_AK4PFPuppi.txt"),
  jecPayloadNamesAK4PFPuppiDATA3   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Fall17_17Nov2017B_V32_DATA/Fall17_17Nov2017B_V32_DATA_L3Absolute_AK4PFPuppi.txt"),
  jecPayloadNamesAK4PFPuppiDATA4   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Fall17_17Nov2017B_V32_DATA/Fall17_17Nov2017B_V32_DATA_L2L3Residual_AK4PFPuppi.txt"),
  jecPayloadNamesAK4PFPuppiDATAUnc = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Fall17_17Nov2017B_V32_DATA/Fall17_17Nov2017B_V32_DATA_Uncertainty_AK4PFPuppi.txt"),
  jecPayloadNamesAK8PFchsMC1   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Fall17_17Nov2017_V32_MC/Fall17_17Nov2017_V32_MC_L1FastJet_AK8PFchs.txt"),
  jecPayloadNamesAK8PFchsMC2   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Fall17_17Nov2017_V32_MC/Fall17_17Nov2017_V32_MC_L2Relative_AK8PFchs.txt"),
  jecPayloadNamesAK8PFchsMC3   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Fall17_17Nov2017_V32_MC/Fall17_17Nov2017_V32_MC_L3Absolute_AK8PFchs.txt"),
  jecPayloadNamesAK8PFchsMCUnc = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Fall17_17Nov2017_V32_MC/Fall17_17Nov2017_V32_MC_Uncertainty_AK8PFchs.txt"),
  jecPayloadNamesAK8PFchsDATA1   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Fall17_17Nov2017B_V32_DATA/Fall17_17Nov2017B_V32_DATA_L1FastJet_AK8PFchs.txt"),
  jecPayloadNamesAK8PFchsDATA2   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Fall17_17Nov2017B_V32_DATA/Fall17_17Nov2017B_V32_DATA_L2Relative_AK8PFchs.txt"),
  jecPayloadNamesAK8PFchsDATA3   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Fall17_17Nov2017B_V32_DATA/Fall17_17Nov2017B_V32_DATA_L3Absolute_AK8PFchs.txt"),
  jecPayloadNamesAK8PFchsDATA4   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Fall17_17Nov2017B_V32_DATA/Fall17_17Nov2017B_V32_DATA_L2L3Residual_AK8PFchs.txt"),
  jecPayloadNamesAK8PFchsDATAUnc = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Fall17_17Nov2017B_V32_DATA/Fall17_17Nov2017B_V32_DATA_Uncertainty_AK8PFchs.txt"),
  PUReweightfile      = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/PUReweight/PileUpReweighting2017.root"),
  PuppiWeightFilePath = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/PUPPI/puppiCorr.root"),
  Pvtx_ndof_min   = cms.double(4.),
  Pvtx_vtx_max    = cms.double(24.),
  Pvtx_vtxdxy_max = cms.double(24.),
  vtx_ndof_min        = cms.int32(4),
  vtx_rho_max         = cms.int32(2),
  vtx_position_z_max  = cms.double(24.),
  Muon_pt_min         = cms.double(10.),
  Muon_eta_max        = cms.double(50),
  patElectron_pt_min  = cms.double(10.),
  patElectron_eta_max = cms.double(50),
  Tau_pt_min          = cms.double(15.),
  Tau_eta_max         = cms.double(50.),
  Jet_pt_min = cms.double(15.),
  Photon_pt_min   = cms.double(5.0),
  Photon_eta_max  = cms.double(5.0),
)
######################################################


#############   Dump gen particle list   #############
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.printGenParticleList = cms.EDAnalyzer("ParticleListDrawer",
  maxEventsToPrint = cms.untracked.int32(-1),
  printVertex = cms.untracked.bool(True),
  src = cms.InputTag("prunedGenParticles")
)
######################################################


####################   PROCESS   #####################
process.p = cms.Path(
process.ecalBadCalibReducedMINIAODFilter*
process.egammaPostRecoSeq *
process.rerunMvaIsolationSequence *
process.NewTauIDsEmbedded *
process.TNT
)
######################################################
