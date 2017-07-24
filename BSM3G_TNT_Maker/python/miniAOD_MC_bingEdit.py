import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
import copy

options = VarParsing.VarParsing('analysis')
#options = VarParsing("python")

# Variables one can control from the multicrab configuration file.
# When connecting a variable you need to tell the module certain information
# about the object.
#                   - Object name.
#                   - Default value.
#                   - Is object a single number or a list.
#                   - Object type.
#                   - Details of object.
#

# ===== Register new variables =====
#options.register('optionTriggerInfo', False,
#VarParsing.VarParsing.multiplicity.singleton,
#VarParsing.VarParsing.varType.bool,
#"Bool for accessing and keeping trigger info in MC")
#options.optionTriggerInfo = True #Set to true for tt+jets

options.register('ofName',
'tempName',
VarParsing.VarParsing.multiplicity.singleton,
VarParsing.VarParsing.varType.string,
"Name for output file."
)

options.register("deterministicSeeds",
True,
VarParsing.VarParsing.multiplicity.singleton,
VarParsing.VarParsing.varType.bool,
"create collections with deterministic seeds" )

options.register("electronRegression",
"GT",
VarParsing.VarParsing.multiplicity.singleton,
VarParsing.VarParsing.varType.string,
"'GT' or an absolute path to a sqlite file for electron energy regression")

options.register("electronSmearing",
"Moriond17_23Jan",
VarParsing.VarParsing.multiplicity.singleton,
VarParsing.VarParsing.varType.string,
"correction type for electron energy smearing" )


# Changed to True for Data configs
options.register("realData",
False,
VarParsing.VarParsing.multiplicity.singleton,
VarParsing.VarParsing.varType.bool,
"input dataset contains real data" )


# ===== Get & parse any command line arguments =====
options.parseArguments()

#####
##   Initial standard configs
#####
process = cms.Process("Demo")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger = cms.Service("MessageLogger", destinations = cms.untracked.vstring('info'), info = cms.untracked.PSet(threshold = cms.untracked.string('DEBUG')), default = cms.untracked.PSet(limit = cms.untracked.int32(-1)), debugModules = cms.untracked.vstring("TNT"))
process.load('Configuration.Geometry.GeometryRecoDB_cff')
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = '80X_mcRun2_asymptotic_2016_TrancheIV_v8'
process.prefer("GlobalTag")
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')


#####
##   Input files
#####
process.source = cms.Source("PoolSource",
  fileNames = cms.untracked.vstring(
    ## ttH synch file
    '/store/mc/RunIISummer16MiniAODv2/ttHToNonbb_M125_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/3C70EB0A-6BBE-E611-B094-0025905A606A.root',
    ##ttjets Synch File
    #'/store/mc/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/0693E0E7-97BE-E611-B32F-0CC47A78A3D8.root'
  ),
  skipEvents = cms.untracked.uint32(0)
)
#>>>> Set limit on number events to process (for testing purposes only) <<<<
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(50) )


#======= Object Collections =======
#Default Collections
electronCollection = cms.InputTag("slimmedElectrons", "", "PAT")
photonCollection   = cms.InputTag("slimmedPhotons", "", "PAT")
muonCollection     = cms.InputTag("slimmedMuons", "", "PAT")
tauCollection      = cms.InputTag("slimmedTaus", "", "PAT")
METCollection      = cms.InputTag("slimmedMETs", "", "PAT")
METPuppiCollection = cms.InputTag("slimmedMETsPuppi")
jetCollection      = cms.InputTag("slimmedJets", "", "PAT")
jetPuppiCollection = cms.InputTag("slimmedJetsPuppi", "", "PAT")
fatjetCollection   = cms.InputTag("slimmedJetsAK8")
topsubjetCollection = cms.InputTag("slimmedJetsCMSTopTagCHSPacked", "SubJets")
#==================================



'''# ===== Implementation of Deterministic seed =====

if options.deterministicSeeds:
    process.load("PhysicsTools.PatUtils.deterministicSeeds_cfi")
    process.deterministicSeeds.produceCollections = cms.bool(True)
    process.deterministicSeeds.produceValueMaps   = cms.bool(False)
    process.deterministicSeeds.electronCollection = electronCollection
    process.deterministicSeeds.muonCollection     = muonCollection
    process.deterministicSeeds.tauCollection      = tauCollection
    process.deterministicSeeds.photonCollection   = photonCollection
    process.deterministicSeeds.jetCollection      = jetCollection
    process.deterministicSeeds.METCollection      = METCollection
    # overwrite output collections
    electronCollection = cms.InputTag("deterministicSeeds", "electronsWithSeed", process.name_())
    muonCollection     = cms.InputTag("deterministicSeeds", "muonsWithSeed", process.name_())
    tauCollection      = cms.InputTag("deterministicSeeds", "tausWithSeed", process.name_())
    photonCollection   = cms.InputTag("deterministicSeeds", "photonsWithSeed", process.name_())
    jetCollection      = cms.InputTag("deterministicSeeds", "jetsWithSeed", process.name_())
    METCollection      = cms.InputTag("deterministicSeeds", "METsWithSeed", process.name_())

# ===== Implementation of Electron Regression =====
if options.electronRegression:
    if options.electronRegression == "GT":
        from EgammaAnalysis.ElectronTools.regressionWeights_cfi import regressionWeights
        process = regressionWeights(process)
    else:
        from EgammaAnalysis.ElectronTools.regressionWeights_local_cfi import GBRDWrapperRcd
        GBRDWrapperRcd.connect = cms.string("sqlite_file:" + options.electronRegression)
        process.regressions = GBRDWrapperRcd
        process.regressions.DumpStat = cms.untracked.bool(False)
        process.es_prefer_regressions = cms.ESPrefer("PoolDBESSource", "regressions")
    process.load("EgammaAnalysis.ElectronTools.regressionApplication_cff")

    # set the electron and photon sources
    process.slimmedElectrons.src = electronCollection
    process.slimmedPhotons.src = photonCollection
    # overwrite output collections
    electronCollection = cms.InputTag("slimmedElectrons", "", process.name_())
    photonCollection = cms.InputTag("slimmedPhotons", "", process.name_())


# ===== Implementation of Electron Smearing =====
if options.electronSmearing and options.electronRegression:
    # the smearing procedure requires a preselection
    process.selectedElectrons = cms.EDFilter("PATElectronSelector",
        src = electronCollection,
        cut = cms.string("pt>5 && abs(superCluster.eta)<2.5")
    )
    electronCollection = cms.InputTag("selectedElectrons", "", process.name_())

    # setup the smearing
    process.load("EgammaAnalysis.ElectronTools.calibratedPatElectronsRun2_cfi")
    from EgammaAnalysis.ElectronTools.calibratedPatElectronsRun2_cfi import files
    process.calibratedPatElectrons.isMC           = cms.bool(not options.realData)
    process.calibratedPatElectrons.correctionFile = cms.string(files[options.electronSmearing])
    process.calibratedPatElectrons.electrons      = electronCollection

    # use our deterministic seeds or a random generator service
    if options.deterministicSeeds:
        process.calibratedPatElectrons.seedUserInt = process.deterministicSeeds.seedUserInt
    else:
        process.load("Configuration.StandardSequences.Services_cff")
        process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
            calibratedPatElectrons = cms.PSet(
                initialSeed = cms.untracked.uint32(81),
                engineName  = cms.untracked.string("TRandom3")
            )
        )

    # overwrite output collections
    electronCollection = cms.InputTag("calibratedPatElectrons", "", process.name_())'''



#####
##   ELECTRON ID SECTION
#####
# Set up everything that is needed to compute electron IDs and
# add the ValueMaps with ID decisions into the event data stream
# Load tools and function definitions
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
switchOnVIDElectronIdProducer(process, DataFormat.MiniAOD)
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Summer16_80X_V1_cff',
                 'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV70_cff',
                 'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring15_25ns_nonTrig_V1_cff',
                 'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring15_25ns_Trig_V1_cff',
                 'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring16_GeneralPurpose_V1_cff',
                 'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring16_HZZ_V1_cff'
                ]
# Add them to the VID producer
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

# update some VID modules to work with potentially changed electron collections
'''process.egmGsfElectronIDs.physicsObjectSrc = electronCollection
process.electronRegressionValueMapProducer.srcMiniAOD = electronCollection
process.electronMVAValueMapProducer.srcMiniAOD = electronCollection'''


#####
##   BTAGGING WITH HIP MITIGATION
#####
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.Geometry.GeometryRecoDB_cff")
from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
updateJetCollection(
  process,
  jetSource = jetCollection,
  jetCorrections = ('AK4PFchs', cms.vstring(['L1FastJet','L2Relative','L3Absolute']), 'None'),
  btagDiscriminators = ['pfCombinedInclusiveSecondaryVertexV2BJetTags', 'pfCombinedMVAV2BJetTags', 'pfJetProbabilityBJetTags', 'pfCombinedCvsLJetTags', 'pfCombinedCvsBJetTags'],
  runIVF=True,
  btagPrefix = 'new' # optional, in case interested in accessing both the old and new discriminator values
)
process.options = cms.untracked.PSet( allowUnscheduled = cms.untracked.bool(True) )




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

#=======================
#===== Output file =====
#=======================
# Must be EDM output files - output files produced by PoolOutputModule in DBS.
# Notice that for publication to be possible, the corresponding output files
# have to be transferred to the permanent storage element.
options.ofName += ".root"
process.TFileService = cms.Service("TFileService",
  fileName = cms.string("OutTree.root")
)

#####
##   Analysis parameters
#####
process.TNT = cms.EDAnalyzer("BSM3G_TNT_Maker",
  #### Running options
  # Choose which trigger you want (do NOT need to put * as it will consider all the versions by default)
  ifevtriggers      = cms.bool(False), # True means you want to require the triggers
  maxtriggerversion = cms.double(10), # please leave it as a double
  evtriggers        = cms.vstring(

    # ttH SL and DL triggers Jan 2017
    'HLT_IsoMu22_v',
    'HLT_IsoTkMu22_v',
    'HLT_IsoMu24_v',
    'HLT_IsoTkMu24_v',
    'HLT_Ele27_eta2p1_WPTight_Gsf_v',
    'HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v',
    'HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v',# Not in ttjets reHLT
    'HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v',
    'HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v',# Not in ttjets reHLT
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
     'HLT_Ele27_WPTight_Gsf_v',
     'HLT_Ele25_eta2p1_WPTight_Gsf_v',
     #Muon
     'HLT_IsoTkMu20_v',
     'HLT_IsoMu22_eta2p1_v',
     'HLT_IsoTkMu22_eta2p1_v',
     #Cross Ele-Mu
     'HLT_DiMu9_Ele9_CaloIdL_TrackIdL_v',
     'HLT_Mu8_DiEle12_CaloIdL_TrackIdL_v',
     'HLT_TripleMu_12_10_5_v',
     'HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v',
     'HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v',
     'HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_v',
     #Other
     'HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v',
     'HLT_Mu50_v',
     'HLT_TkMu50_v',
     'HLT_Ele115_CaloIdVT_GsfTrkIdT_v',
  ),
  # Choose which information you want to use
  fillgeninfo           = cms.bool(True),
  fillgenHFCategoryinfo = cms.bool(True),
  filleventinfo         = cms.bool(True),
  filltriggerinfo       = cms.bool(True), #F, for samples without trigger
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
  analysisFilter        = cms.bool(True),
  # Choose format
  MiniAODv2 = cms.bool(True),
  is_data   = cms.bool(False),
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
  muons               = muonCollection,
  patElectrons        = electronCollection,
  electronVetoIdMap   = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-veto"),
  electronLooseIdMap  = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-loose"),
  electronMediumIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-medium"),
  electronTightIdMap  = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-tight"),
  eleMVATrigIdMap        = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-Trig-V1-wp80"),
  eleMVAnonTrigIdMap     = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-nonTrig-V1-wp80"),
  eleMVATrigwp90IdMap    = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-Trig-V1-wp90"),
  eleMVAnonTrigwp90IdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-nonTrig-V1-wp90"),
  eleMVAGPwp90IdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring16-GeneralPurpose-V1-wp90"),
  eleMVAGPwp80IdMap  = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring16-GeneralPurpose-V1-wp80"),
  eleMVAHZZwpLooseIdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring16-HZZ-V1-wpLoose"),
  eleHEEPIdMap                 = cms.InputTag("egmGsfElectronIDs:heepElectronID-HEEPV70"),
  elemvaValuesMap_nonTrig      = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Values"),
  elemvaCategoriesMap_nonTrig  = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Categories"),
  elemvaValuesMap_Trig         = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15Trig25nsV1Values"),
  elemvaCategoriesMap_Trig     = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15Trig25nsV1Categories"),
  elemvaValuesMap_GP           = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring16GeneralPurposeV1Values"),
  elemvaCategoriesMap_GP       = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring16GeneralPurposeV1Categories"),
  elemvaValuesMap_HZZ          = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring16HZZV1Values"),
  elemvaCategoriesMap_HZZ      = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring16HZZV1Categories"),
  taus                = tauCollection,
  jets                = jetCollection,
  jetsPUPPI           = jetPuppiCollection,
  fatjets             = fatjetCollection,
  topsubjets          = topsubjetCollection,
  mets                = METCollection,
  metsPUPPI           = METPuppiCollection,
  metFilterBits       = cms.InputTag("TriggerResults", "", "PAT"),
  photons             = photonCollection,
  packedPFCandidates  = cms.InputTag("packedPFCandidates"),
  pruned              = cms.InputTag("prunedGenParticles"),
  # =========== JER (Only applied to MC) ============
  jerAK4PFchs     =  cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JER/Spring16_25nsV10_MC_PtResolution_AK4PFchs.txt"),
  jerAK4PFchsSF   =  cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JER/Spring16_25nsV10_MC_SF_AK4PFchs.txt"),
  jerAK4PFPuppi   =  cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JER/Spring16_25nsV10_MC_PtResolution_AK4PFPuppi.txt"),
  jerAK4PFPuppiSF =  cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JER/Spring16_25nsV10_MC_SF_AK4PFPuppi.txt"),
  jerAK8PFchs     =  cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JER/Spring16_25nsV10_MC_PtResolution_AK8PFchs.txt"),
  jerAK8PFchsSF   =  cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JER/Spring16_25nsV10_MC_SF_AK8PFchs.txt"),
  jerAK8PFPuppi   =  cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JER/Spring16_25nsV10_MC_PtResolution_AK8PFPuppi.txt"),
  jerAK8PFPuppiSF =  cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JER/Spring16_25nsV10_MC_SF_AK8PFchs.txt"),
  # ===========  JEC - CORRECTIONS ON FLY ===========
  #=== MC ===
  # L1FastJet
  # L2Relative
  # L3Absolute
  # MC_Uncertainty
  jecPayloadNamesAK4PFchsMC1   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Summer16_23Sep2016V4_MC/Summer16_23Sep2016V4_MC_L1FastJet_AK4PFchs.txt"),
  jecPayloadNamesAK4PFchsMC2   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Summer16_23Sep2016V4_MC/Summer16_23Sep2016V4_MC_L2Relative_AK4PFchs.txt"),
  jecPayloadNamesAK4PFchsMC3   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Summer16_23Sep2016V4_MC/Summer16_23Sep2016V4_MC_L3Absolute_AK4PFchs.txt"),
  #jecPayloadNamesAK4PFchsMCUnc = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Summer16_23Sep2016V4_MC/Summer16_23Sep2016V4_MC_Uncertainty_AK4PFchs.txt"),
  jecPayloadNamesAK4PFchsMCUnc = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Summer16_23Sep2016V4_MC/Summer16_23Sep2016V4_MC_UncertaintySources_AK4PFchs.txt"),
  jecPayloadNamesAK4PFPuppiMC1   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Summer16_23Sep2016V4_MC/Summer16_23Sep2016V4_MC_L1FastJet_AK4PFPuppi.txt"),
  jecPayloadNamesAK4PFPuppiMC2   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Summer16_23Sep2016V4_MC/Summer16_23Sep2016V4_MC_L2Relative_AK4PFPuppi.txt"),
  jecPayloadNamesAK4PFPuppiMC3   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Summer16_23Sep2016V4_MC/Summer16_23Sep2016V4_MC_L3Absolute_AK4PFPuppi.txt"),
  #jecPayloadNamesAK4PFPuppiMCUnc = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Summer16_23Sep2016V4_MC/Summer16_23Sep2016V4_MC_Uncertainty_AK4PFPuppi.txt"),
  jecPayloadNamesAK4PFPuppiMCUnc = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Summer16_23Sep2016V4_MC/Summer16_23Sep2016V4_MC_UncertaintySources_AK4PFPuppi.txt"),
  jecPayloadNamesAK8PFchsMC1   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Summer16_23Sep2016V4_MC/Summer16_23Sep2016V4_MC_L1FastJet_AK8PFchs.txt"),
  jecPayloadNamesAK8PFchsMC2   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Summer16_23Sep2016V4_MC/Summer16_23Sep2016V4_MC_L2Relative_AK8PFchs.txt"),
  jecPayloadNamesAK8PFchsMC3   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Summer16_23Sep2016V4_MC/Summer16_23Sep2016V4_MC_L3Absolute_AK8PFchs.txt"),
  #jecPayloadNamesAK8PFchsMCUnc = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Summer16_23Sep2016V4_MC/Summer16_23Sep2016V4_MC_Uncertainty_AK8PFchs.txt"),
  jecPayloadNamesAK8PFchsMCUnc = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Summer16_23Sep2016V4_MC/Summer16_23Sep2016V4_MC_UncertaintySources_AK8PFchs.txt"),
  #=== DATA ===
  # Note that this is a config for MC.
  # The following Data JEC will not be
  # used but are needed to prevent the
  # BSMFramework code from crashing.
  jecPayloadNamesAK4PFchsDATA1   = cms.FileInPath('BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_L1FastJet_AK4PFchs.txt'),
  jecPayloadNamesAK4PFchsDATA2   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_L2Relative_AK4PFchs.txt"),
  jecPayloadNamesAK4PFchsDATA3   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_L3Absolute_AK4PFchs.txt"),
  jecPayloadNamesAK4PFchsDATA4   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_L2L3Residual_AK4PFchs.txt"),
  #jecPayloadNamesAK4PFchsDATAUnc = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_Uncertainty_AK4PFchs.txt"),
  jecPayloadNamesAK4PFchsDATAUnc = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_UncertaintySources_AK4PFchs.txt"),
  jecPayloadNamesAK8PFchsDATA1   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_L1FastJet_AK4PFPuppi.txt"),
  jecPayloadNamesAK8PFchsDATA2   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_L2Relative_AK4PFPuppi.txt"),
  jecPayloadNamesAK8PFchsDATA3   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_L3Absolute_AK4PFPuppi.txt"),
  jecPayloadNamesAK8PFchsDATA4   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_L2L3Residual_AK4PFPuppi.txt"),
  #jecPayloadNamesAK8PFchsDATAUnc = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_Uncertainty_AK4PFPuppi.txt"),
  jecPayloadNamesAK8PFchsDATAUnc = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_UncertaintySources_AK4PFPuppi.txt"),
  jecPayloadNamesAK4PFPuppiDATA1   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_L1FastJet_AK8PFchs.txt"),
  jecPayloadNamesAK4PFPuppiDATA2   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_L2Relative_AK8PFchs.txt"),
  jecPayloadNamesAK4PFPuppiDATA3   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_L3Absolute_AK8PFchs.txt"),
  jecPayloadNamesAK4PFPuppiDATA4   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_L2L3Residual_AK8PFchs.txt"),
  #jecPayloadNamesAK4PFPuppiDATAUnc = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_Uncertainty_AK8PFchs.txt"),
  jecPayloadNamesAK4PFPuppiDATAUnc = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_UncertaintySources_AK8PFchs.txt"),
  # PILEUP REWEIGHTING
  PUReweightfile      = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/PUReweight/PileUpReweightingMoriond17.root"),
  # BTAG REWEIGHTING
  BTAGReweightfile1   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/BTAGReweight/2017csvSFs/csv_rwt_fit_hf_v2_final_2017_1_10test.root"),
  BTAGReweightfile2   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/BTAGReweight/2017csvSFs/csv_rwt_fit_lf_v2_final_2017_1_10test.root"),
  BTAGReweightfile3   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/BTAGReweight/2017csvSFs/CSVv2_Moriond17_B_H.csv"),
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
  Muon_pt_min         = cms.double(5.), #25.),#5.),
  Muon_eta_max        = cms.double(50), #2.1),#50),
  # Electron cuts
  patElectron_pt_min  = cms.double(5.), #30.),#5.),
  patElectron_eta_max = cms.double(50), #2.1),#50),
  # Tau cuts
  Tau_pt_min          = cms.double(15.), #30),#15.),
  Tau_eta_max         = cms.double(50), #2.1),#50.),
  # Jet cuts
  Jet_pt_min = cms.double(10.),#20.),#10.),
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
  badGlobalMuonTagger       = cms.InputTag("badGlobalMuonTagger"),
  cloneGlobalMuonTagger     = cms.InputTag("cloneGlobalMuonTagger"),

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
process.BJetness.is_data = cms.bool(False)
process.BJetness.vertices = cms.InputTag("offlineSlimmedPrimaryVertices")
process.BJetness.eleMVAnonTrigIdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-nonTrig-V1-wp80")
process.BJetness.electronTightIdMap  = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-tight")
process.BJetness.patElectrons = electronCollection
process.BJetness.muons = muonCollection
process.BJetness.jets = jetCollection
process.BJetness.jecPayloadNamesAK4PFchsMC1 = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Summer16_23Sep2016V4_MC/Summer16_23Sep2016V4_MC_L1FastJet_AK4PFchs.txt")
process.BJetness.jecPayloadNamesAK4PFchsMC2 = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Summer16_23Sep2016V4_MC/Summer16_23Sep2016V4_MC_L2Relative_AK4PFchs.txt")
process.BJetness.jecPayloadNamesAK4PFchsMC3 = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Summer16_23Sep2016V4_MC/Summer16_23Sep2016V4_MC_L3Absolute_AK4PFchs.txt")
process.BJetness.jecPayloadNamesAK4PFchsMCUnc = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Summer16_23Sep2016V4_MC/Summer16_23Sep2016V4_MC_UncertaintySources_AK4PFchs.txt")
process.BJetness.jecPayloadNamesAK4PFchsDATA1 = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_L1FastJet_AK4PFchs.txt")
process.BJetness.jecPayloadNamesAK4PFchsDATA2 = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_L2Relative_AK4PFchs.txt")
process.BJetness.jecPayloadNamesAK4PFchsDATA3 = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_L3Absolute_AK4PFchs.txt")
process.BJetness.jecPayloadNamesAK4PFchsDATA4 = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_L2L3Residual_AK4PFchs.txt")
process.BJetness.jecPayloadNamesAK4PFchsDATAUnc = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_UncertaintySources_AK4PFchs.txt")
process.BJetness.jerAK4PFchs   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JER/Spring16_25nsV10_MC_PtResolution_AK4PFchs.txt")
process.BJetness.jerAK4PFchsSF = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JER/Spring16_25nsV10_MC_SF_AK4PFchs.txt")
#QG likelihood
process.load('BSMFramework.BSM3G_TNT_Maker.QGTagger_cfi')
process.QGTagger.srcJets       = jetCollection
process.QGTagger.jetsLabel     = cms.string('QGL_AK4PFchs')


# Access badGlobalMuonTagger and cloneGlobalMuonTagger
# "_cff" files need to be associated with the process before you can access them/their methods via your process object.
# In order to do so one needs to make a clone of the python script template:

process.load("RecoMET.METFilters.BadPFMuonFilter_cfi")
from RecoMET.METFilters.BadPFMuonFilter_cfi import BadPFMuonFilter
process.BadPFMuonFilter = BadPFMuonFilter.clone(
    muons = muonCollection,
    PFCandidates = cms.InputTag("packedPFCandidates"),
)


process.load("RecoMET.METFilters.BadChargedCandidateFilter_cfi")
from RecoMET.METFilters.BadChargedCandidateFilter_cfi import BadChargedCandidateFilter
process.BadChargedCandidateFilter = BadChargedCandidateFilter.clone(
    muons = muonCollection,
    PFCandidates = cms.InputTag("packedPFCandidates"),
)

process.load("RecoMET.METFilters.badGlobalMuonTaggersMiniAOD_cff")
from RecoMET.METFilters.badGlobalMuonTaggersMiniAOD_cff import badGlobalMuonTaggerMAOD

process.badGlobalMuonTagger = badGlobalMuonTaggerMAOD.clone(
    taggingMode = cms.bool(True),
)

from RecoMET.METFilters.badGlobalMuonTaggersMiniAOD_cff import cloneGlobalMuonTaggerMAOD
process.cloneGlobalMuonTagger = cloneGlobalMuonTaggerMAOD.clone(
    taggingMode = cms.bool(True),
)

#Run analysis sequence
process.p = cms.Path(
#process.regressionApplication
#*process.calibratedPatElectrons
process.selectedHadronsAndPartons
*process.genJetFlavourInfos
*process.matchGenCHadron
*process.matchGenBHadron
*process.egmGsfElectronIDSequence
*process.BJetness
*process.BadPFMuonFilter
*process.BadChargedCandidateFilter
*process.badGlobalMuonTagger
*process.cloneGlobalMuonTagger
*process.TNT
)
#process.e = cms.EndPath(process.out)
