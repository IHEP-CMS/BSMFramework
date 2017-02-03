import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
import copy
options = VarParsing.VarParsing('analysis')
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
options.register('optionTriggerInfo', False,
VarParsing.VarParsing.multiplicity.singleton,
VarParsing.VarParsing.varType.bool,
"Bool for accessing and keeping trigger info in MC")
options.optionTriggerInfo = True #Set to true for tt+jets

options.register('ofName',
'',
VarParsing.VarParsing.multiplicity.singleton,
VarParsing.VarParsing.varType.string,
"Name for output file."
)

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
process.GlobalTag.globaltag = '80X_mcRun2_asymptotic_2016_TrancheIV_v7'
process.prefer("GlobalTag")
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
#####
##   Input files
#####
process.source = cms.Source("PoolSource",
  fileNames = cms.untracked.vstring(
    #ttjets Synch File
    '/store/mc/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/0693E0E7-97BE-E611-B32F-0CC47A78A3D8.root'
  ),
  skipEvents = cms.untracked.uint32(0)
)

#>>>> Set limit on number events to process (for testing purposes only) <<<<
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

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
switchOnVIDElectronIdProducer(process, DataFormat.MiniAOD)
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Summer16_80X_V1_cff',
                 'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV70_cff',
                 'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring15_25ns_nonTrig_V1_cff',
                 'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring15_25ns_Trig_V1_cff'
                ]
# Add them to the VID producer
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

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
  fileName = cms.string(options.ofName)
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
  fillgeninfo           = cms.bool(True),
  fillgenHFCategoryinfo = cms.bool(True),
  filleventinfo         = cms.bool(True),
  filltriggerinfo       = cms.bool(options.optionTriggerInfo), #F, for samples without trigger
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
  muons               = cms.InputTag("slimmedMuons"),
  patElectrons        = cms.InputTag("slimmedElectrons"),
  electronVetoIdMap   = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-veto"),
  electronLooseIdMap  = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-loose"),
  electronMediumIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-medium"),
  electronTightIdMap  = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-tight"),
  eleMVATrigIdMap        = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-Trig-V1-wp80"),
  eleMVAnonTrigIdMap     = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-nonTrig-V1-wp80"),
  eleMVATrigwp90IdMap    = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-Trig-V1-wp90"),
  eleMVAnonTrigwp90IdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-nonTrig-V1-wp90"),
  eleHEEPIdMap                 = cms.InputTag("egmGsfElectronIDs:heepElectronID-HEEPV70"),
  elemvaValuesMap_nonTrig      = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Values"),
  elemvaCategoriesMap_nonTrig  = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Categories"),
  elemvaValuesMap_Trig         = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15Trig25nsV1Values"),
  elemvaCategoriesMap_Trig     = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15Trig25nsV1Categories"),
  taus                = cms.InputTag("slimmedTaus"),
  #jets                = cms.InputTag("selectedPatJetsAK8PFCHS"),
  jets                = cms.InputTag("slimmedJets"), #selectedUpdatedPatJets #slimmedJets"),
  jetsPUPPI           = cms.InputTag("slimmedJetsPuppi"),
  fatjets             = cms.InputTag("slimmedJetsAK8"),
  topsubjets          = cms.InputTag("slimmedJetsCMSTopTagCHSPacked", "SubJets"),
  mets                = cms.InputTag("slimmedMETs"),
  metsPUPPI           = cms.InputTag("slimmedMETsPuppi"),
  metFilterBits       = cms.InputTag("TriggerResults", "", "PAT"),
  photons             = cms.InputTag("slimmedPhotons"),
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
  jecPayloadNamesAK4PFchsMC1   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Spring16_23Sep2016V2_MC/Spring16_23Sep2016V2_MC_L1FastJet_AK4PFchs.txt"),
  jecPayloadNamesAK4PFchsMC2   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Spring16_23Sep2016V2_MC/Spring16_23Sep2016V2_MC_L2Relative_AK4PFchs.txt"),
  jecPayloadNamesAK4PFchsMC3   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Spring16_23Sep2016V2_MC/Spring16_23Sep2016V2_MC_L3Absolute_AK4PFchs.txt"),
  jecPayloadNamesAK4PFchsMCUnc = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Spring16_23Sep2016V2_MC/Spring16_23Sep2016V2_MC_Uncertainty_AK4PFchs.txt"),
  jecPayloadNamesAK4PFPuppiMC1   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Spring16_23Sep2016V2_MC/Spring16_23Sep2016V2_MC_L1FastJet_AK4PFPuppi.txt"),
  jecPayloadNamesAK4PFPuppiMC2   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Spring16_23Sep2016V2_MC/Spring16_23Sep2016V2_MC_L2Relative_AK4PFPuppi.txt"),
  jecPayloadNamesAK4PFPuppiMC3   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Spring16_23Sep2016V2_MC/Spring16_23Sep2016V2_MC_L3Absolute_AK4PFPuppi.txt"),
  jecPayloadNamesAK4PFPuppiMCUnc = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Spring16_23Sep2016V2_MC/Spring16_23Sep2016V2_MC_Uncertainty_AK4PFPuppi.txt"),
  jecPayloadNamesAK8PFchsMC1   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Spring16_23Sep2016V2_MC/Spring16_23Sep2016V2_MC_L1FastJet_AK8PFchs.txt"),
  jecPayloadNamesAK8PFchsMC2   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Spring16_23Sep2016V2_MC/Spring16_23Sep2016V2_MC_L2Relative_AK8PFchs.txt"),
  jecPayloadNamesAK8PFchsMC3   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Spring16_23Sep2016V2_MC/Spring16_23Sep2016V2_MC_L3Absolute_AK8PFchs.txt"),
  jecPayloadNamesAK8PFchsMCUnc = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Spring16_23Sep2016V2_MC/Spring16_23Sep2016V2_MC_Uncertainty_AK8PFchs.txt"),
  #=== DATA ===
  # Note that this is a config for MC.
  # The following Data JEC will not be
  # used but are needed to prevent the
  # BSMFramework code from crashing.
  jecPayloadNamesAK4PFchsDATA1   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Spring16_25nsV6_DATA_L1FastJet_AK4PFchs.txt"),
  jecPayloadNamesAK4PFchsDATA2   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Spring16_25nsV6_DATA_L2Relative_AK4PFchs.txt"),
  jecPayloadNamesAK4PFchsDATA3   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Spring16_25nsV6_DATA_L3Absolute_AK4PFchs.txt"),
  jecPayloadNamesAK4PFchsDATA4   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Spring16_25nsV6_DATA_L2L3Residual_AK4PFchs.txt"),
  jecPayloadNamesAK4PFchsDATAUnc = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Spring16_25nsV6_DATA_Uncertainty_AK4PFchs.txt"),
  jecPayloadNamesAK8PFchsDATA1   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Spring16_25nsV6_DATA_L1FastJet_AK8PFchs.txt"),
  jecPayloadNamesAK8PFchsDATA2   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Spring16_25nsV6_DATA_L2Relative_AK8PFchs.txt"),
  jecPayloadNamesAK8PFchsDATA3   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Spring16_25nsV6_DATA_L3Absolute_AK8PFchs.txt"),
  jecPayloadNamesAK8PFchsDATA4   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Spring16_25nsV6_DATA_L2L3Residual_AK8PFchs.txt"),
  jecPayloadNamesAK8PFchsDATAUnc = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Spring16_25nsV6_DATA_Uncertainty_AK8PFchs.txt"),
  jecPayloadNamesAK4PFPuppiDATA1   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Spring16_25nsV6_DATA_L1FastJet_AK4PFPuppi.txt"),
  jecPayloadNamesAK4PFPuppiDATA2   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Spring16_25nsV6_DATA_L2Relative_AK4PFPuppi.txt"),
  jecPayloadNamesAK4PFPuppiDATA3   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Spring16_25nsV6_DATA_L3Absolute_AK4PFPuppi.txt"),
  jecPayloadNamesAK4PFPuppiDATA4   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Spring16_25nsV6_DATA_L2L3Residual_AK4PFPuppi.txt"),
  jecPayloadNamesAK4PFPuppiDATAUnc = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Spring16_25nsV6_DATA_Uncertainty_AK4PFPuppi.txt"),
  # PILEUP REWEIGHTING
  PUReweightfile      = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/PUReweight/PileUpReweightingMoriond17.root"),
  # BTAG REWEIGHTING
  BTAGReweightfile1   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/BTAGReweight/2017csvSFs/csv_rwt_fit_hf_v2_final_2017_1_10test.root"),
  BTAGReweightfile2   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/BTAGReweight/2017csvSFs/csv_rwt_fit_lf_v2_final_2017_1_10test.root"),
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
process.BJetness.patElectrons = cms.InputTag("slimmedElectrons")
process.BJetness.muons = cms.InputTag("slimmedMuons")
process.BJetness.jets = cms.InputTag("slimmedJets")#selectedUpdatedPatJets #slimmedJets"),
process.BJetness.jecPayloadNamesAK4PFchsMC1 = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Spring16_23Sep2016V2_MC/Spring16_23Sep2016V2_MC_L1FastJet_AK4PFchs.txt")
process.BJetness.jecPayloadNamesAK4PFchsMC2 = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Spring16_23Sep2016V2_MC/Spring16_23Sep2016V2_MC_L2Relative_AK4PFchs.txt")
process.BJetness.jecPayloadNamesAK4PFchsMC3 = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Spring16_23Sep2016V2_MC/Spring16_23Sep2016V2_MC_L3Absolute_AK4PFchs.txt")
process.BJetness.jecPayloadNamesAK4PFchsMCUnc = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Spring16_23Sep2016V2_MC/Spring16_23Sep2016V2_MC_Uncertainty_AK4PFchs.txt")
process.BJetness.jecPayloadNamesAK4PFchsDATA1 = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Spring16_25nsV6_DATA_L1FastJet_AK4PFchs.txt")
process.BJetness.jecPayloadNamesAK4PFchsDATA2 = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Spring16_25nsV6_DATA_L2Relative_AK4PFchs.txt")
process.BJetness.jecPayloadNamesAK4PFchsDATA3 = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Spring16_25nsV6_DATA_L3Absolute_AK4PFchs.txt")
process.BJetness.jecPayloadNamesAK4PFchsDATA4 = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Spring16_25nsV6_DATA_L2L3Residual_AK4PFchs.txt")
process.BJetness.jecPayloadNamesAK4PFchsDATAUnc = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Spring16_25nsV6_DATA_Uncertainty_AK4PFchs.txt")
process.BJetness.jerAK4PFchs   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JER/Spring16_25nsV10_MC_PtResolution_AK4PFchs.txt")
process.BJetness.jerAK4PFchsSF = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JER/Spring16_25nsV10_MC_SF_AK4PFchs.txt")
#QG likelihood
process.load('BSMFramework.BSM3G_TNT_Maker.QGTagger_cfi')
process.QGTagger.srcJets       = cms.InputTag('slimmedJets') #selectedUpdatedPatJets #slimmedJets')
process.QGTagger.jetsLabel     = cms.string('QGL_AK4PFchs')
#Run analysis sequence
process.p = cms.Path(
process.selectedHadronsAndPartons*process.genJetFlavourInfos*process.matchGenCHadron*process.matchGenBHadron*
process.egmGsfElectronIDSequence*
process.BJetness*
process.TNT
)
#process.e = cms.EndPath(process.out)
