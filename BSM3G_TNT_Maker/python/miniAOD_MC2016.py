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
process.GlobalTag.globaltag = '94X_mcRun2_asymptotic_v3'
process.prefer("GlobalTag")
######################################################


##################   Input files   ###################
process.source = cms.Source("PoolSource",
  fileNames = cms.untracked.vstring(
    '/store/mc/RunIISummer16MiniAODv3/TprimeBToTZ_M-900_LH_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3-v2/10000/2C0646E3-650D-E911-AA52-003048FFD772.root'
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
                       runVID=True,
                       applyEnergyCorrections=False,
                       era='2016-Legacy')
######################################################


####################  TAU v2   #######################
from BSMFramework.BSM3G_TNT_Maker.runTauIdMVA import *
na = TauIDEmbedder(process, cms,
        debug=True,
        toKeep = ["2017v2", "newDM2017v2", "dR0p32017v2"]
        )
na.runTauID()
######################################################


###################  PREFIRING   #####################
process.prefiringweight = cms.EDProducer("L1ECALPrefiringWeightProducer",
                                 ThePhotons = cms.InputTag("slimmedPhotons"),
	                         TheJets = cms.InputTag("slimmedJets"),
                                 L1Maps = cms.string("L1PrefiringMaps_new.root"),
                                 DataEra = cms.string("2016BtoH"), #Use 2016BtoH for 2016
                                 UseJetEMPt = cms.bool(False), #can be set to true to use jet prefiring maps parametrized vs pt(em) instead of pt
	                         PrefiringRateSystematicUncty = cms.double(0.2) #Minimum relative prefiring uncty per object
)
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
  fillgeninfo           = cms.bool(True),
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
  is_data   = cms.bool(False),
  debug_    = cms.bool(False),
  PuppiVar  = cms.bool(False),
  MC2016    = cms.bool(True),
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
  metFilterBits       = cms.InputTag("TriggerResults", "", "PAT"),
  photons             = cms.InputTag("slimmedPhotons"),
  packedPFCandidates  = cms.InputTag("packedPFCandidates"), 
  pruned              = cms.InputTag("prunedGenParticles"),
  jerAK4PFchs     =  cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JER/Summer16_25nsV1_MC_PtResolution_AK4PFchs.txt"),
  jerAK4PFchsSF   =  cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JER/Summer16_25nsV1_MC_SF_AK4PFchs.txt"),
  jerAK4PFPuppi   =  cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JER/Summer16_25nsV1_MC_PtResolution_AK4PFPuppi.txt"),
  jerAK4PFPuppiSF =  cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JER/Summer16_25nsV1_MC_SF_AK4PFPuppi.txt"),
  jerAK8PFchs     =  cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JER/Summer16_25nsV1_MC_PtResolution_AK8PFchs.txt"),
  jerAK8PFchsSF   =  cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JER/Summer16_25nsV1_MC_SF_AK8PFchs.txt"),
  jerAK8PFPuppi   =  cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JER/Summer16_25nsV1_MC_PtResolution_AK8PFPuppi.txt"),
  jerAK8PFPuppiSF =  cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JER/Summer16_25nsV1_MC_SF_AK8PFPuppi.txt"),
  jecPayloadNamesAK4PFchsMC1   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Summer16_07Aug2017_V11_MC/Summer16_07Aug2017_V11_MC_L1FastJet_AK4PFchs.txt"),
  jecPayloadNamesAK4PFchsMC2   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Summer16_07Aug2017_V11_MC/Summer16_07Aug2017_V11_MC_L2Relative_AK4PFchs.txt"),
  jecPayloadNamesAK4PFchsMC3   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Summer16_07Aug2017_V11_MC/Summer16_07Aug2017_V11_MC_L3Absolute_AK4PFchs.txt"),
  jecPayloadNamesAK4PFchsMCUnc = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Summer16_07Aug2017_V11_MC/Summer16_07Aug2017_V11_MC_Uncertainty_AK4PFchs.txt"),
  jecPayloadNamesAK4PFchsDATA1   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer16_07Aug2017BCD_V11_DATA/Summer16_07Aug2017BCD_V11_DATA_L1FastJet_AK4PFchs.txt"),
  jecPayloadNamesAK4PFchsDATA2   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer16_07Aug2017BCD_V11_DATA/Summer16_07Aug2017BCD_V11_DATA_L2Relative_AK4PFchs.txt"),
  jecPayloadNamesAK4PFchsDATA3   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer16_07Aug2017BCD_V11_DATA/Summer16_07Aug2017BCD_V11_DATA_L3Absolute_AK4PFchs.txt"),
  jecPayloadNamesAK4PFchsDATA4   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer16_07Aug2017BCD_V11_DATA/Summer16_07Aug2017BCD_V11_DATA_L2L3Residual_AK4PFchs.txt"),
  jecPayloadNamesAK4PFchsDATAUnc = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer16_07Aug2017BCD_V11_DATA/Summer16_07Aug2017BCD_V11_DATA_Uncertainty_AK4PFchs.txt"),
  jecPayloadNamesAK4PFPuppiMC1   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Summer16_07Aug2017_V11_MC/Summer16_07Aug2017_V11_MC_L1FastJet_AK4PFPuppi.txt"),
  jecPayloadNamesAK4PFPuppiMC2   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Summer16_07Aug2017_V11_MC/Summer16_07Aug2017_V11_MC_L2Relative_AK4PFPuppi.txt"),
  jecPayloadNamesAK4PFPuppiMC3   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Summer16_07Aug2017_V11_MC/Summer16_07Aug2017_V11_MC_L3Absolute_AK4PFPuppi.txt"),
  jecPayloadNamesAK4PFPuppiMCUnc = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Summer16_07Aug2017_V11_MC/Summer16_07Aug2017_V11_MC_Uncertainty_AK4PFPuppi.txt"),
  jecPayloadNamesAK4PFPuppiDATA1   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer16_07Aug2017BCD_V11_DATA/Summer16_07Aug2017BCD_V11_DATA_L1FastJet_AK4PFPuppi.txt"),
  jecPayloadNamesAK4PFPuppiDATA2   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer16_07Aug2017BCD_V11_DATA/Summer16_07Aug2017BCD_V11_DATA_L2Relative_AK4PFPuppi.txt"),
  jecPayloadNamesAK4PFPuppiDATA3   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer16_07Aug2017BCD_V11_DATA/Summer16_07Aug2017BCD_V11_DATA_L3Absolute_AK4PFPuppi.txt"),
  jecPayloadNamesAK4PFPuppiDATA4   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer16_07Aug2017BCD_V11_DATA/Summer16_07Aug2017BCD_V11_DATA_L2L3Residual_AK4PFPuppi.txt"),
  jecPayloadNamesAK4PFPuppiDATAUnc = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer16_07Aug2017BCD_V11_DATA/Summer16_07Aug2017BCD_V11_DATA_Uncertainty_AK4PFPuppi.txt"),
  jecPayloadNamesAK8PFchsMC1   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Summer16_07Aug2017_V11_MC/Summer16_07Aug2017_V11_MC_L1FastJet_AK8PFchs.txt"),
  jecPayloadNamesAK8PFchsMC2   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Summer16_07Aug2017_V11_MC/Summer16_07Aug2017_V11_MC_L2Relative_AK8PFchs.txt"),
  jecPayloadNamesAK8PFchsMC3   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Summer16_07Aug2017_V11_MC/Summer16_07Aug2017_V11_MC_L3Absolute_AK8PFchs.txt"),
  jecPayloadNamesAK8PFchsMCUnc = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Summer16_07Aug2017_V11_MC/Summer16_07Aug2017_V11_MC_Uncertainty_AK8PFchs.txt"),
  jecPayloadNamesAK8PFchsDATA1   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer16_07Aug2017BCD_V11_DATA/Summer16_07Aug2017BCD_V11_DATA_L1FastJet_AK8PFchs.txt"),
  jecPayloadNamesAK8PFchsDATA2   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer16_07Aug2017BCD_V11_DATA/Summer16_07Aug2017BCD_V11_DATA_L2Relative_AK8PFchs.txt"),
  jecPayloadNamesAK8PFchsDATA3   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer16_07Aug2017BCD_V11_DATA/Summer16_07Aug2017BCD_V11_DATA_L3Absolute_AK8PFchs.txt"),
  jecPayloadNamesAK8PFchsDATA4   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer16_07Aug2017BCD_V11_DATA/Summer16_07Aug2017BCD_V11_DATA_L2L3Residual_AK8PFchs.txt"),
  jecPayloadNamesAK8PFchsDATAUnc = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer16_07Aug2017BCD_V11_DATA/Summer16_07Aug2017BCD_V11_DATA_Uncertainty_AK8PFchs.txt"),
  PUReweightfile      = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/PUReweight/PileUpReweighting2016.root"),
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
process.egammaPostRecoSeq *
process.rerunMvaIsolationSequence *
process.NewTauIDsEmbedded *
process.prefiringweight *
process.TNT
)
######################################################
