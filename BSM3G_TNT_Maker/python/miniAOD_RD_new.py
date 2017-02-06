import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
import copy
from PhysicsTools.PatAlgos.patEventContent_cff import patEventContent

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
options.register('optionJECAK4PFchsDATA1',
'',
VarParsing.VarParsing.multiplicity.singleton,
VarParsing.VarParsing.varType.string,
"Name for AK4PFchsDATA1 JEC file")
#options.optionJECAK4PFchsDATA1 = "BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Spring16_23Sep2016HV2_DATA/Spring16_23Sep2016HV2_DATA_L1FastJet_AK4PFchs.txt"

options.register('optionJECAK4PFchsDATA2',
'',
VarParsing.VarParsing.multiplicity.singleton,
VarParsing.VarParsing.varType.string,
"Name for AK4PFchsDATA2 JEC file"
)
options.register('optionJECAK4PFchsDATA3',
'',
VarParsing.VarParsing.multiplicity.singleton,
VarParsing.VarParsing.varType.string,
"Name for AK4PFchsDATA3 JEC file"
)
options.register('optionJECAK4PFchsDATA4',
'',
VarParsing.VarParsing.multiplicity.singleton,
VarParsing.VarParsing.varType.string,
"Name for AK4PFchsDATA4 JEC file"
)
options.register('optionJECAK4PFchsDATAUnc',
'',
VarParsing.VarParsing.multiplicity.singleton,
VarParsing.VarParsing.varType.string,
"Name for AK4PFchsDATAUnc JEC file"
)


options.register('optionJECAK4PFPuppiDATA1',
'',
VarParsing.VarParsing.multiplicity.singleton,
VarParsing.VarParsing.varType.string,
"Name for AK4PFPuppiDATA1 JEC file"
)
options.register('optionJECAK4PFPuppiDATA2',
'',
VarParsing.VarParsing.multiplicity.singleton,
VarParsing.VarParsing.varType.string,
"Name for AK4PFPuppiDATA2 JEC file"
)
options.register('optionJECAK4PFPuppiDATA3',
'',
VarParsing.VarParsing.multiplicity.singleton,
VarParsing.VarParsing.varType.string,
"Name for AK4PFPuppiDATA3 JEC file"
)
options.register('optionJECAK4PFPuppiDATA4',
'',
VarParsing.VarParsing.multiplicity.singleton,
VarParsing.VarParsing.varType.string,
"Name for AK4PFPuppiDATA4 JEC file"
)
options.register('optionJECAK4PFPuppiDATAUnc',
'',
VarParsing.VarParsing.multiplicity.singleton,
VarParsing.VarParsing.varType.string,
"Name for AK4PFPuppiDATAUnc JEC file"
)



options.register('optionJECAK8PFchsDATA1',
'',
VarParsing.VarParsing.multiplicity.singleton,
VarParsing.VarParsing.varType.string,
"Name for AK8PFchsDATA1 JEC file"
)
options.register('optionJECAK8PFchsDATA2',
'',
VarParsing.VarParsing.multiplicity.singleton,
VarParsing.VarParsing.varType.string,
"Name for AK8PFchsDATA2 JEC file"
)
options.register('optionJECAK8PFchsDATA3',
'',
VarParsing.VarParsing.multiplicity.singleton,
VarParsing.VarParsing.varType.string,
"Name for AK8PFchsDATA3 JEC file"
)
options.register('optionJECAK8PFchsDATA4',
'',
VarParsing.VarParsing.multiplicity.singleton,
VarParsing.VarParsing.varType.string,
"Name for AK8PFchsDATA4 JEC file"
)
options.register('optionJECAK8PFchsDATAUnc',
'',
VarParsing.VarParsing.multiplicity.singleton,
VarParsing.VarParsing.varType.string,
"Name for AK8PFchsDATAUnc JEC file"
)

options.register('ofName',
'',
VarParsing.VarParsing.multiplicity.singleton,
VarParsing.VarParsing.varType.string,
"Name for output file."
)


# ===== Get & parse any command line arguments =====
options.parseArguments()



# =====   Standard config variables =====
process = cms.Process("Demo")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10
process.load('Configuration.Geometry.GeometryRecoDB_cff')
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = '80X_dataRun2_Prompt_v15'
process.prefer("GlobalTag")
process.load('Configuration.StandardSequences.Services_cff')
#####
##   Input files
#####
process.source = cms.Source("PoolSource",
  fileNames = cms.untracked.vstring(

    #Data, e
     '/store/data/Run2016B/SingleElectron/MINIAOD/23Sep2016-v3/00000/00099863-E799-E611-A876-141877343E6D.root',
     '/store/data/Run2016B/SingleElectron/MINIAOD/23Sep2016-v3/00000/0077B9D8-1E9A-E611-963A-0CC47A57CCE8.root',
     '/store/data/Run2016C/SingleElectron/MINIAOD/23Sep2016-v1/50000/001B31C0-248C-E611-B7BD-0025905B860E.root',
     '/store/data/Run2016C/SingleElectron/MINIAOD/23Sep2016-v1/50000/001C375D-098B-E611-9453-0025905B85E8.root',
     '/store/data/Run2016D/SingleElectron/MINIAOD/23Sep2016-v1/70000/04E8F72C-AF89-E611-9D2F-FA163E1D7951.root',
     '/store/data/Run2016D/SingleElectron/MINIAOD/23Sep2016-v1/70000/0654EC37-FA8A-E611-820F-FA163E066046.root',

     #Data, mu
     '/store/data/Run2016B/SingleMuon/MINIAOD/23Sep2016-v3/00000/00AE0629-1F98-E611-921A-008CFA1112CC.root',
     '/store/data/Run2016B/SingleMuon/MINIAOD/23Sep2016-v3/00000/041C2515-1698-E611-B909-901B0E542804.root',
     '/store/data/Run2016C/SingleMuon/MINIAOD/23Sep2016-v1/70000/001F13A2-7E8D-E611-B910-FA163E782438.root',
     '/store/data/Run2016C/SingleMuon/MINIAOD/23Sep2016-v1/70000/00DFF8E4-828C-E611-AF50-0CC47A0AD6AA.root',
     '/store/data/Run2016D/SingleMuon/MINIAOD/23Sep2016-v1/010000/0201B90C-C79B-E611-8763-00266CFFBF80.root',
     '/store/data/Run2016D/SingleMuon/MINIAOD/23Sep2016-v1/010000/0CEAC70F-C79B-E611-BE1E-00266CFFC7CC.root'

  ),
  skipEvents = cms.untracked.uint32(0)
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )
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
  #fileName = cms.string(options.ofName)
  fileName = cms.string("OutTree.root")
)

## Output Module Configuration (expects a path 'p')
#process.out = cms.OutputModule("PoolOutputModule",
    #fileName = cms.untracked.string(options.ofName),
#    fileName = cms.untracked.string("test_output.root"),
    ## save only events passing the full path
#    SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
    ## save PAT output; you need a '*' to unpack the list of commands
    ##'patEventContent'
#    outputCommands = cms.untracked.vstring('drop *', *patEventContent)
#)


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
     'HLT_DoubleEle33_CaloIdL_MW_v',
     'HLT_Mu50_v',
     'HLT_TkMu50_v',
     'HLT_Ele115_CaloIdVT_GsfTrkIdT_v',
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
  analysisFilter        = cms.bool(True),
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
   jecPayloadNamesAK4PFchsDATA1   = cms.FileInPath(options.optionJECAK4PFchsDATA1),
   jecPayloadNamesAK4PFchsDATA2   = cms.FileInPath(options.optionJECAK4PFchsDATA2),
   jecPayloadNamesAK4PFchsDATA3   = cms.FileInPath(options.optionJECAK4PFchsDATA3),
   jecPayloadNamesAK4PFchsDATA4   = cms.FileInPath(options.optionJECAK4PFchsDATA4),
   jecPayloadNamesAK4PFchsDATAUnc = cms.FileInPath(options.optionJECAK4PFchsDATAUnc),
   jecPayloadNamesAK4PFPuppiDATA1   = cms.FileInPath(options.optionJECAK4PFPuppiDATA1),
   jecPayloadNamesAK4PFPuppiDATA2   = cms.FileInPath(options.optionJECAK4PFPuppiDATA2),
   jecPayloadNamesAK4PFPuppiDATA3   = cms.FileInPath(options.optionJECAK4PFPuppiDATA3),
   jecPayloadNamesAK4PFPuppiDATA4   = cms.FileInPath(options.optionJECAK4PFPuppiDATA4),
   jecPayloadNamesAK4PFPuppiDATAUnc = cms.FileInPath(options.optionJECAK4PFPuppiDATAUnc),
   jecPayloadNamesAK8PFchsDATA1   = cms.FileInPath(options.optionJECAK8PFchsDATA1),
   jecPayloadNamesAK8PFchsDATA2   = cms.FileInPath(options.optionJECAK8PFchsDATA2),
   jecPayloadNamesAK8PFchsDATA3   = cms.FileInPath(options.optionJECAK8PFchsDATA3),
   jecPayloadNamesAK8PFchsDATA4   = cms.FileInPath(options.optionJECAK8PFchsDATA4),
   jecPayloadNamesAK8PFchsDATAUnc = cms.FileInPath(options.optionJECAK8PFchsDATAUnc),

  # PILEUP REWEIGHTING
  PUReweightfile      = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/PUReweight/PileUpReweighting2016.root"),
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
process.BJetness.jecPayloadNamesAK4PFchsDATA1 = cms.FileInPath(options.optionJECAK4PFchsDATA1)
process.BJetness.jecPayloadNamesAK4PFchsDATA2 = cms.FileInPath(options.optionJECAK4PFchsDATA2)
process.BJetness.jecPayloadNamesAK4PFchsDATA3 = cms.FileInPath(options.optionJECAK4PFchsDATA3)
process.BJetness.jecPayloadNamesAK4PFchsDATA4 = cms.FileInPath(options.optionJECAK4PFchsDATA4)
process.BJetness.jecPayloadNamesAK4PFchsDATAUnc = cms.FileInPath(options.optionJECAK4PFchsDATAUnc)
process.BJetness.jerAK4PFchs   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JER/Spring16_25nsV10_MC_PtResolution_AK4PFchs.txt")
process.BJetness.jerAK4PFchsSF = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JER/Spring16_25nsV10_MC_SF_AK4PFchs.txt")
#QG likelihood
process.load('BSMFramework.BSM3G_TNT_Maker.QGTagger_cfi')
process.QGTagger.srcJets       = cms.InputTag('slimmedJets')
process.QGTagger.jetsLabel     = cms.string('QGL_AK4PFchs')
#Run analysis sequence
process.p = cms.Path(
process.egmGsfElectronIDSequence*
process.BJetness*
process.TNT
)
#process.outpath = cms.EndPath(process.out)
#del process.out
