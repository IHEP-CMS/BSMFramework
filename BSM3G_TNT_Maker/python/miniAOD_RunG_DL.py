import FWCore.ParameterSet.Config as cms
#####
##   Initial standard configs
#####
process = cms.Process("Demo")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10000
#process.load("Configuration.StandardSequences.Geometry_cff")
##process.load('Configuration.Geometry.GeometryIdeal_cff')
process.load('Configuration.Geometry.GeometryRecoDB_cff')
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
#process.load("Configuration.StandardSequences.MagneticField_38T_PostLS1_cff")
#process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")

#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = '80X_dataRun2_Prompt_v15'
process.prefer("GlobalTag")
process.load('Configuration.StandardSequences.Services_cff')
#process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
#####
##   Input files
#####
process.source = cms.Source("PoolSource",
  fileNames = cms.untracked.vstring(

    #Data, e e
    '/store/data/Run2016G/DoubleEG/MINIAOD/23Sep2016-v1/100000/0608426D-3F8E-E611-A52D-00237DF28460.root',
    '/store/data/Run2016G/DoubleEG/MINIAOD/23Sep2016-v1/100000/06FC0D6B-F88A-E611-A24D-00215E2EB6EE.root',

    #Data, e mu
    '/store/data/Run2016G/MuonEG/MINIAOD/23Sep2016-v1/100000/005AB7E9-0B93-E611-AC81-848F69FD2925.root',
    '/store/data/Run2016G/MuonEG/MINIAOD/23Sep2016-v1/100000/005FB5C1-ED8F-E611-BAE7-0025905A607E.root',

    #Data, mu mu
    '/store/data/Run2016G/DoubleMuon/MINIAOD/23Sep2016-v1/100000/00993A51-DF90-E611-A4EE-7845C4FC3650.root',
    '/store/data/Run2016G/DoubleMuon/MINIAOD/23Sep2016-v1/100000/00DD00F8-008C-E611-8CD0-00266CFFC9C4.root'

  ),
  skipEvents = cms.untracked.uint32(0)
)
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )
#####
##   BTAGGING WITH HIP MITIGATION
#####
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.Geometry.GeometryRecoDB_cff")
from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
updateJetCollection(
  process,
  jetSource = cms.InputTag('slimmedJets'),
  jetCorrections = ('AK4PFchs', cms.vstring(['L1FastJet','L2Relative','L3Absolute']), 'None'),
  btagDiscriminators = ['pfCombinedInclusiveSecondaryVertexV2BJetTags', 'pfCombinedMVAV2BJetTags', 'pfJetProbabilityBJetTags', 'pfCombinedCvsLJetTags', 'pfCombinedCvsBJetTags'],
  runIVF=True,
  btagPrefix = 'new' # optional, in case interested in accessing both the old and new discriminator values
)
process.options = cms.untracked.PSet( allowUnscheduled = cms.untracked.bool(True) )
#####
##   ELECTRON ID SECTION
#####
# Set up everything that is needed to compute electron IDs and
# add the ValueMaps with ID decisions into the event data stream
# Load tools and function definitions
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
#process.load("RecoEgamma.ElectronIdentification.ElectronMVAValueMapProducer_cfi")
#process.load("RecoEgamma.ElectronIdentification.egmGsfElectronIDs_cfi")
#process.egmGsfElectronIDs.physicsObjectSrc = cms.InputTag('slimmedElectrons')
#from PhysicsTools.SelectorUtils.centralIDRegistry import central_id_registry
#process.egmGsfElectronIDSequence = cms.Sequence(process.egmGsfElectronIDs)
switchOnVIDElectronIdProducer(process, DataFormat.MiniAOD)
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Spring15_25ns_V1_cff',
                 'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV60_cff',
                 'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring15_25ns_nonTrig_V1_cff',
                 'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring15_25ns_Trig_V1_cff'
                ]
# Add them to the VID producer
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)
"""
#####
##   JEC (to check if they need to be used in miniAOD)
#####
## JEC
#from JetMETCorrections.Configuration.JetCorrectionServices_cff import *
#process.ak4PFCHSL1Fastjet = cms.ESProducer(
#  'L1FastjetCorrectionESProducer',
#  level       = cms.string('L1FastJet'),
#  algorithm   = cms.string('AK4PFchs'),
#  srcRho      = cms.InputTag( 'fixedGridRhoFastjetAll' ),
#  useCondDB = cms.untracked.bool(True)
#)
#process.ak4PFchsL2Relative = ak4CaloL2Relative.clone( algorithm = 'AK4PFchs' )
#process.ak4PFchsL3Absolute = ak4CaloL3Absolute.clone( algorithm = 'AK4PFchs' )
#process.ak4PFchsResidual   = ak5PFResidual.clone( algorithm = 'AK4PFchs' )
#process.ak4PFchsL1L2L3     = cms.ESProducer("JetCorrectionESChain",
#  correctors = cms.vstring('ak4PFCHSL1Fastjet','ak4PFchsL2Relative','ak4PFchsL3Absolute'),#,'ak4PFchsResidual'),
#  useCondDB = cms.untracked.bool(True)
#)
## JEC (a la TTHLep)
#from RecoJets.Configuration.RecoJets_cff import *
#from RecoJets.Configuration.RecoPFJets_cff import *
#from JetMETCorrections.Configuration.JetCorrectionProducersAllAlgos_cff import *
#from JetMETCorrections.Configuration.JetCorrectionServicesAllAlgos_cff import *
#from JetMETCorrections.Configuration.JetCorrectionServices_cff import *
#process.ak4PFCHSL1Fastjet = cms.ESProducer(
#  'L1FastjetCorrectionESProducer',
#  level       = cms.string('L1FastJet'),
#  algorithm   = cms.string('AK4PFchs'),
#  srcRho      = cms.InputTag( 'fixedGridRhoFastjetCentralNeutral' ), #'fixedGridRhoFastjetAll' ),
#  useCondDB = cms.untracked.bool(True)
#  )
#process.ak4PFchsL2Relative  =  ak5PFL2Relative.clone( algorithm = 'AK4PFchs' )
#process.ak4PFchsL3Absolute  =  ak5PFL3Absolute.clone( algorithm = 'AK4PFchs' )
#process.ak4PFchsResidual    =  ak5PFResidual.clone(   algorithm = 'AK4PFchs' )
#process.ak4PFCHSL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
#  correctors = cms.vstring(
#  'ak4PFCHSL1Fastjet',
#  'ak4PFchsL2Relative',
#  'ak4PFchsL3Absolute',
#  'ak4PFchsResidual'),
#  useCondDB = cms.untracked.bool(True)
#)
#####
##   MET (to check if they need to be used in miniAOD)
#####
# filter out anomalous MET from detector noise, cosmic rays, and beam halo particles
#process.load("RecoMET.METFilters.metFilters_cff")
#process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
#  vertexCollection = cms.InputTag('offlineSlimmedPrimaryVertices'),
#  minimumNDOF = cms.uint32(4) ,
#  maxAbsZ = cms.double(24),
#  maxd0 = cms.double(2)
#)
###___________________________HCAL_Noise_Filter________________________________||
#process.load('CommonTools.RecoAlgos.HBHENoiseFilterResultProducer_cfi')
#process.HBHENoiseFilterResultProducer.minZeros = cms.int32(99999)
#process.ApplyBaselineHBHENoiseFilter = cms.EDFilter('BooleanFlagFilter',
#  inputLabel = cms.InputTag('HBHENoiseFilterResultProducer','HBHENoiseFilterResult'),
#  reverseDecision = cms.bool(False)
#)
"""
#####
##   For tt+X
#####
# Setting input particle collections to be used by the tools
genParticleCollection = "prunedGenParticles"
genJetCollection      = "slimmedGenJets"
jetFlavourInfos       = "genJetFlavourInfos"
jetAlgo               = "AntiKt"
rParam                = 0.4
genJetPtMin           = 20.
genJetAbsEtaMax       = 2.4
from PhysicsTools.JetMCAlgos.HadronAndPartonSelector_cfi import selectedHadronsAndPartons
process.selectedHadronsAndPartons = selectedHadronsAndPartons.clone(
  particles = genParticleCollection
)
from PhysicsTools.JetMCAlgos.AK4PFJetsMCFlavourInfos_cfi import ak4JetFlavourInfos
process.genJetFlavourInfos = ak4JetFlavourInfos.clone(
    jets         = genJetCollection,
    rParam       = cms.double(rParam),
    jetAlgorithm = cms.string(jetAlgo)
)
from PhysicsTools.JetMCAlgos.GenHFHadronMatcher_cff import matchGenBHadron
process.matchGenBHadron = matchGenBHadron.clone(
  genParticles = genParticleCollection,
  jetFlavourInfos = jetFlavourInfos
)
from PhysicsTools.JetMCAlgos.GenHFHadronMatcher_cff import matchGenCHadron
process.matchGenCHadron = matchGenCHadron.clone(
  genParticles = genParticleCollection,
  jetFlavourInfos = jetFlavourInfos
)
#####
##   Output file
#####
process.TFileService = cms.Service("TFileService",
  fileName = cms.string("Data_RunG_DL.root")
)
#####
##   Analysis parameters
#####
process.TNT = cms.EDAnalyzer("BSM3G_TNT_Maker",
  #### Running options
  # Choose which trigger you want (do NOT need to put * as it will consider all the versions by default)
  ifevtriggers      = cms.bool(True), # True means you want to require the triggers
  maxtriggerversion = cms.double(10), # please leave it as a double
  evtriggers        = cms.vstring(
    # ttH SL and DL triggers Jan 2017
    'HLT_IsoMu22_v',
    'HLT_IsoTkMu22_v',
    'HLT_IsoMu24_v',
    'HLT_IsoTkMu24_v',
    'HLT_Ele27_eta2p1_WPTight_Gsf_v',
    'HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v',
    'HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v',
    'HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v',
    'HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v',
    'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v',
    'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v',
    'HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v',
    'HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v',
    'HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v',
     #Electron
     'HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v',
     #Muon
     'HLT_IsoMu20_v',
     'HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v',
     'HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v',
    #TTHLep
     #Electron
     'HLT_Ele23_WPLoose_Gsf_v', #Data
     'HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v', #MC
     #Muon
     'HLT_IsoTkMu20_v',
     #Cross Ele-Mu
     'HLT_DiMu9_Ele9_CaloIdL_TrackIdL_v',
     'HLT_Mu8_DiEle12_CaloIdL_TrackIdL_v',
     'HLT_TripleMu_12_10_5_v',
     'HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v',
     #Other
     'HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v',
  ),
  # Choose which information you want to use
  fillgeninfo           = cms.bool(False), #F
  fillgenHFCategoryinfo = cms.bool(False), #F
  filleventinfo         = cms.bool(True),
  filltriggerinfo       = cms.bool(True),
  fillPVinfo            = cms.bool(True),
  fillmuoninfo          = cms.bool(True),
  fillelectronpatinfo   = cms.bool(True),
  filltauinfo           = cms.bool(True),
  filljetinfo           = cms.bool(True),
  filltthjetinfo        = cms.bool(False), #F
  fillBoostedJetinfo    = cms.bool(True),
  fillTopSubJetinfo     = cms.bool(False), #F
  fillTauJetnessinfo    = cms.bool(True),
  fillBJetnessinfo      = cms.bool(True),
  fillBJetnessFVinfo    = cms.bool(True),
  fillBTagReweight      = cms.bool(True),
  fillPileupReweight    = cms.bool(True),
  fillMETinfo           = cms.bool(True),
  fillphotoninfo        = cms.bool(False), #F
  # Choose format
  MiniAODv2 = cms.bool(True),
  is_data   = cms.bool(True),
  debug_    = cms.bool(False),
  super_TNT = cms.bool(False),
  AJVar     = cms.bool(False),
  tthlepVar = cms.bool(True),
  bjetnessselfilter = cms.bool(False),
  bjetnessproducer  = cms.bool(True),
  PuppiVar  = cms.bool(False),
  qglVar    = cms.bool(True),
  # Input tags
  bits                = cms.InputTag("TriggerResults","","HLT"),
  prescales           = cms.InputTag("patTrigger"),
  objects             = cms.InputTag("selectedPatTrigger"),
  vertices            = cms.InputTag("offlineSlimmedPrimaryVertices"),
  beamSpot            = cms.InputTag("offlineBeamSpot"),
  muons               = cms.InputTag("slimmedMuons"),
  patElectrons        = cms.InputTag("slimmedElectrons"),
  electronVetoIdMap   = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-veto"),
  electronLooseIdMap  = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-loose"),
  electronMediumIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-medium"),
  electronTightIdMap  = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-tight"),
  eleMVATrigIdMap        = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-Trig-V1-wp80"),
  eleMVAnonTrigIdMap     = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-nonTrig-V1-wp80"),
  eleMVATrigwp90IdMap    = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-Trig-V1-wp90"),
  eleMVAnonTrigwp90IdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-nonTrig-V1-wp90"),
  eleHEEPIdMap                 = cms.InputTag("egmGsfElectronIDs:heepElectronID-HEEPV60"),
  elemvaValuesMap_nonTrig      = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Values"),
  elemvaCategoriesMap_nonTrig  = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Categories"),
  elemvaValuesMap_Trig         = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15Trig25nsV1Values"),
  elemvaCategoriesMap_Trig     = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15Trig25nsV1Categories"),
  taus                = cms.InputTag("slimmedTaus"),
  #jets                = cms.InputTag("selectedPatJetsAK8PFCHS"),
  jets                = cms.InputTag("slimmedJets"),#selectedUpdatedPatJets #slimmedJets"),
  jetsPUPPI           = cms.InputTag("slimmedJetsPuppi"),
  fatjets             = cms.InputTag("slimmedJetsAK8"),
  topsubjets          = cms.InputTag("slimmedJetsCMSTopTagCHSPacked", "SubJets"),
  mets                = cms.InputTag("slimmedMETs"),
  metsPUPPI           = cms.InputTag("slimmedMETsPuppi"),
  metFilterBits       = cms.InputTag("TriggerResults", "", "RECO"),
  photons             = cms.InputTag("slimmedPhotons"),
  packedPFCandidates  = cms.InputTag("packedPFCandidates"),
  pruned              = cms.InputTag("prunedGenParticles"),
  # ======== JER (only applied to MC.) =========
  # Note that this is a config for Data.
  # The following MC JER/JEC will not be
  # used but are needed to prevent the
  # BSMFramework code from crashing.
 jerAK4PFchs     =  cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JER/JRDatabase/textFiles/Spring16_25nsV10_MC/Spring16_25nsV10_MC_PtResolution_AK4PFchs.txt"),
 jerAK4PFchsSF   =  cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JER/JRDatabase/textFiles/Spring16_25nsV10_MC/Spring16_25nsV10_MC_SF_AK4PFchs.txt"),
 jerAK4PFPuppi   =  cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JER/JRDatabase/textFiles/Spring16_25nsV10_MC/Spring16_25nsV10_MC_PtResolution_AK4PFPuppi.txt"),
 jerAK4PFPuppiSF =  cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JER/JRDatabase/textFiles/Spring16_25nsV10_MC/Spring16_25nsV10_MC_SF_AK4PFPuppi.txt"),
 jerAK8PFchs     =  cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JER/JRDatabase/textFiles/Spring16_25nsV10_MC/Spring16_25nsV10_MC_PtResolution_AK8PFchs.txt"),
 jerAK8PFchsSF   =  cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JER/JRDatabase/textFiles/Spring16_25nsV10_MC/Spring16_25nsV10_MC_SF_AK8PFchs.txt"),
 jerAK8PFPuppi   =  cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JER/JRDatabase/textFiles/Spring16_25nsV10_MC/Spring16_25nsV10_MC_PtResolution_AK8PFPuppi.txt"),
 jerAK8PFPuppiSF =  cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JER/JRDatabase/textFiles/Spring16_25nsV10_MC/Spring16_25nsV10_MC_SF_AK8PFchs.txt"),
  # ===========  JEC - CORRECTIONS ON FLY ===========
  #=== MC ===
  # L1FastJet
  # L2Relative
  # L3Absolute
  # MC_Uncertainty
  jecPayloadNamesAK4PFchsMC1   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Spring16_25nsV6_MC_L1FastJet_AK4PFchs.txt"),
  jecPayloadNamesAK4PFchsMC2   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Spring16_25nsV6_MC_L2Relative_AK4PFchs.txt"),
  jecPayloadNamesAK4PFchsMC3   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Spring16_25nsV6_MC_L3Absolute_AK4PFchs.txt"),
  jecPayloadNamesAK4PFchsMCUnc = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Spring16_25nsV6_MC_Uncertainty_AK4PFchs.txt"),
  jecPayloadNamesAK4PFPuppiMC1   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Spring16_25nsV6_MC_L1FastJet_AK4PFPuppi.txt"),
  jecPayloadNamesAK4PFPuppiMC2   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Spring16_25nsV6_MC_L2Relative_AK4PFPuppi.txt"),
  jecPayloadNamesAK4PFPuppiMC3   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Spring16_25nsV6_MC_L3Absolute_AK4PFPuppi.txt"),
  jecPayloadNamesAK4PFPuppiMCUnc = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Spring16_25nsV6_MC_Uncertainty_AK4PFPuppi.txt"),
  jecPayloadNamesAK8PFchsMC1   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Spring16_25nsV6_MC_L1FastJet_AK8PFchs.txt"),
  jecPayloadNamesAK8PFchsMC2   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Spring16_25nsV6_MC_L2Relative_AK8PFchs.txt"),
  jecPayloadNamesAK8PFchsMC3   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Spring16_25nsV6_MC_L3Absolute_AK8PFchs.txt"),
  jecPayloadNamesAK8PFchsMCUnc = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Spring16_25nsV6_MC_Uncertainty_AK8PFchs.txt"),
  #=== DATA ===
  jecPayloadNamesAK4PFchsDATA1   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Spring16_23Sep2016GV2_DATA/Spring16_23Sep2016GV2_DATA_L1FastJet_AK4PFchs.txt"),
  jecPayloadNamesAK4PFchsDATA2   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Spring16_23Sep2016GV2_DATA/Spring16_23Sep2016GV2_DATA_L2Relative_AK4PFchs.txt"),
  jecPayloadNamesAK4PFchsDATA3   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Spring16_23Sep2016GV2_DATA/Spring16_23Sep2016GV2_DATA_L3Absolute_AK4PFchs.txt"),
  jecPayloadNamesAK4PFchsDATA4   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Spring16_23Sep2016GV2_DATA/Spring16_23Sep2016GV2_DATA_L2L3Residual_AK4PFchs.txt"),
  jecPayloadNamesAK4PFchsDATAUnc = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Spring16_23Sep2016GV2_DATA/Spring16_23Sep2016GV2_DATA_Uncertainty_AK4PFchs.txt"),
  jecPayloadNamesAK4PFPuppiDATA1   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Spring16_23Sep2016GV2_DATA/Spring16_23Sep2016GV2_DATA_L1FastJet_AK4PFPuppi.txt"),
  jecPayloadNamesAK4PFPuppiDATA2   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Spring16_23Sep2016GV2_DATA/Spring16_23Sep2016GV2_DATA_L2Relative_AK4PFPuppi.txt"),
  jecPayloadNamesAK4PFPuppiDATA3   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Spring16_23Sep2016GV2_DATA/Spring16_23Sep2016GV2_DATA_L3Absolute_AK4PFPuppi.txt"),
  jecPayloadNamesAK4PFPuppiDATA4   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Spring16_23Sep2016GV2_DATA/Spring16_23Sep2016GV2_DATA_L2L3Residual_AK4PFPuppi.txt"),
  jecPayloadNamesAK4PFPuppiDATAUnc = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Spring16_23Sep2016GV2_DATA/Spring16_23Sep2016GV2_DATA_Uncertainty_AK4PFPuppi.txt"),
  jecPayloadNamesAK8PFchsDATA1   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Spring16_23Sep2016GV2_DATA/Spring16_23Sep2016GV2_DATA_L1FastJet_AK8PFchs.txt"),
  jecPayloadNamesAK8PFchsDATA2   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Spring16_23Sep2016GV2_DATA/Spring16_23Sep2016GV2_DATA_L2Relative_AK8PFchs.txt"),
  jecPayloadNamesAK8PFchsDATA3   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Spring16_23Sep2016GV2_DATA/Spring16_23Sep2016GV2_DATA_L3Absolute_AK8PFchs.txt"),
  jecPayloadNamesAK8PFchsDATA4   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Spring16_23Sep2016GV2_DATA/Spring16_23Sep2016GV2_DATA_L2L3Residual_AK8PFchs.txt"),
  jecPayloadNamesAK8PFchsDATAUnc = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Spring16_23Sep2016GV2_DATA/Spring16_23Sep2016GV2_DATA_Uncertainty_AK8PFchs.txt"),

  # PILEUP REWEIGHTING
  PUReweightfile      = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/PUReweight/PileUpReweighting2016.root"),
  # BTAG REWEIGHTING
  BTAGReweightfile1   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/BTAGReweight/csv_rwt_fit_hf_v2_final_2016_09_7test.root"),
  BTAGReweightfile2   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/BTAGReweight/csv_rwt_fit_lf_v2_final_2016_09_7test.root"),
  # Object selection
  # Primary vertex cuts
  Pvtx_ndof_min   = cms.double(4.),
  Pvtx_vtx_max    = cms.double(24.),
  Pvtx_vtxdxy_max = cms.double(24.),
  # Obj primary vertex cuts
  vtx_ndof_min        = cms.int32(4),
  vtx_rho_max         = cms.int32(2),
  vtx_position_z_max  = cms.double(24.),
  # Muon cuts
  Muon_pt_min         = cms.double(5.),
  Muon_eta_max        = cms.double(50),
  # Electron cuts
  patElectron_pt_min  = cms.double(5.),
  patElectron_eta_max = cms.double(50),
  # Tau cuts
  Tau_pt_min          = cms.double(15.),
  Tau_eta_max         = cms.double(50.),
  # Jet cuts
  Jet_pt_min = cms.double(10.),
  # Photon cuts
  Photon_pt_min   = cms.double(5.0),
  Photon_eta_max  = cms.double(5.0),
  # ttHFCategorization
  genJetPtMin               = cms.double(genJetPtMin),
  genJetAbsEtaMax           = cms.double(genJetAbsEtaMax),
  genJets                   = cms.InputTag(genJetCollection),
  genBHadJetIndex           = cms.InputTag("matchGenBHadron", "genBHadJetIndex"),
  genBHadFlavour            = cms.InputTag("matchGenBHadron", "genBHadFlavour"),
  genBHadFromTopWeakDecay   = cms.InputTag("matchGenBHadron", "genBHadFromTopWeakDecay"),
  genBHadPlusMothers        = cms.InputTag("matchGenBHadron", "genBHadPlusMothers"),
  genBHadPlusMothersIndices = cms.InputTag("matchGenBHadron", "genBHadPlusMothersIndices"),
  genBHadIndex              = cms.InputTag("matchGenBHadron", "genBHadIndex"),
  genBHadLeptonHadronIndex  = cms.InputTag("matchGenBHadron", "genBHadLeptonHadronIndex"),
  genBHadLeptonViaTau       = cms.InputTag("matchGenBHadron", "genBHadLeptonViaTau"),
  genCHadJetIndex           = cms.InputTag("matchGenCHadron", "genCHadJetIndex"),
  genCHadFlavour            = cms.InputTag("matchGenCHadron", "genCHadFlavour"),
  genCHadFromTopWeakDecay   = cms.InputTag("matchGenCHadron", "genCHadFromTopWeakDecay"),
  genCHadBHadronId          = cms.InputTag("matchGenCHadron", "genCHadBHadronId"),
)
#####
##   Dump gen particle list
#####
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.printGenParticleList = cms.EDAnalyzer("ParticleListDrawer",
  maxEventsToPrint = cms.untracked.int32(-1),
  printVertex = cms.untracked.bool(True),
  src = cms.InputTag("prunedGenParticles")
)
#process.p = cms.Path(process.printGenParticleList)
#BJetness producer
process.load('BJetnessTTHbb.BJetness.BJetness_cfi')
process.BJetness = cms.EDProducer('BJetness')
process.BJetness.is_data = cms.bool(True)
process.BJetness.vertices = cms.InputTag("offlineSlimmedPrimaryVertices")
process.BJetness.eleMVAnonTrigIdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-nonTrig-V1-wp80")
process.BJetness.patElectrons = cms.InputTag("slimmedElectrons")
process.BJetness.muons = cms.InputTag("slimmedMuons")
process.BJetness.jets = cms.InputTag("slimmedJets")#selectedUpdatedPatJets #slimmedJets"),
process.BJetness.jecPayloadNamesAK4PFchsMC1 = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Spring16_25nsV6_MC_L1FastJet_AK4PFchs.txt")
process.BJetness.jecPayloadNamesAK4PFchsMC2 = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Spring16_25nsV6_MC_L2Relative_AK4PFchs.txt")
process.BJetness.jecPayloadNamesAK4PFchsMC3 = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Spring16_25nsV6_MC_L3Absolute_AK4PFchs.txt")
process.BJetness.jecPayloadNamesAK4PFchsMCUnc = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Spring16_25nsV6_MC_Uncertainty_AK4PFchs.txt")
process.BJetness.jecPayloadNamesAK4PFchsDATA1 = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Spring16_23Sep2016GV2_DATA/Spring16_23Sep2016GV2_DATA_L1FastJet_AK4PFchs.txt")
process.BJetness.jecPayloadNamesAK4PFchsDATA2 = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Spring16_23Sep2016GV2_DATA/Spring16_23Sep2016GV2_DATA_L2Relative_AK4PFchs.txt")
process.BJetness.jecPayloadNamesAK4PFchsDATA3 = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Spring16_23Sep2016GV2_DATA/Spring16_23Sep2016GV2_DATA_L3Absolute_AK4PFchs.txt")
process.BJetness.jecPayloadNamesAK4PFchsDATA4 = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Spring16_23Sep2016GV2_DATA/Spring16_23Sep2016GV2_DATA_L2L3Residual_AK4PFchs.txt")
process.BJetness.jecPayloadNamesAK4PFchsDATAUnc = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Spring16_23Sep2016GV2_DATA/Spring16_23Sep2016GV2_DATA_Uncertainty_AK4PFchs.txt")
process.BJetness.jerAK4PFchs   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JER/Spring16_25nsV6_MC_PtResolution_AK4PFchs.txt")
process.BJetness.jerAK4PFchsSF = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JER/Spring16_25nsV6_MC_SF_AK4PFchs.txt")
#QG likelihood
process.load('BSMFramework.BSM3G_TNT_Maker.QGTagger_cfi')
process.QGTagger.srcJets       = cms.InputTag('slimmedJets')#selectedUpdatedPatJets #slimmedJets
process.QGTagger.jetsLabel     = cms.string('QGL_AK4PFchs')
#Run analysis sequence
process.p = cms.Path(
#process.selectedHadronsAndPartons*process.genJetFlavourInfos*process.matchGenCHadron*process.matchGenBHadron*
process.egmGsfElectronIDSequence*
process.BJetness*
#process.primaryVertexFilter*
#process.CSCTightHaloFilter*process.eeBadScFilter*process.HBHENoiseFilterResultProducer*process.ApplyBaselineHBHENoiseFilter*
process.TNT
)
