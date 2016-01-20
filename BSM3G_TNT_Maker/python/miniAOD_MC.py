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
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
process.GlobalTag.globaltag = '74X_mcRun2_asymptotic_v4'
process.prefer("GlobalTag")
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
#####
##   Input files
#####
process.source = cms.Source("PoolSource",
  fileNames = cms.untracked.vstring(
    #TPrime b -> tZb (M=1.0TeV)
    #'/store/mc/RunIISpring15DR74/TprimeBToTZ_M-1000_LH_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/60000/1036477B-3A3A-E511-B6B8-002590593920.root',
    #TTHbb
    #'/store/mc/RunIISpring15DR74/ttHTobb_M125_13TeV_powheg_pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/00000/141B9915-1F08-E511-B9FF-001E675A6AB3.root',
#'/store/mc/RunIISpring15DR74/ttHTobb_M125_13TeV_powheg_pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/00000/1440BF72-A308-E511-9812-90B11C1453E1.root',
    #TT
    #'/store/mc/RunIISpring15DR74/TT_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v2/00000/0AB045B5-BB0C-E511-81FD-0025905A60B8.root'
    #DYJetsToLL
    #'/store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v3/10000/009D49A5-7314-E511-84EF-0025905A605E.root'
    #SUSYGluGluToHToTauTau_M-160
    #'/store/mc/RunIISpring15DR74/SUSYGluGluToHToTauTau_M-160_TuneCUETP8M1_13TeV-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/10000/2A3929AE-5303-E511-9EFE-0025905A48C0.root'
    #Zprime
    #v1
    #'/store/mc/RunIISpring15DR74/ZprimeToTauTau_M_1500_TuneCUETP8M1_tauola_13TeV_pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/30000/20773840-892F-E511-A40C-002590AC4BF8.root',
    #v2
    #'root://eoscms.cern.ch//eos/cms/store/mc/RunIISpring15MiniAODv2/ZprimeToTauTau_M_2000_TuneCUETP8M1_tauola_13TeV_pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/80000/9E637F69-DC74-E511-A619-00259074AE98.root'
    #HN
    #'/store/mc/RunIISpring15MiniAODv2/ExtendedWeakIsospinModel_mumujj_L15000_M3500_CalcHEP/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/20000/0CE3AC18-1B80-E511-BD72-44A842B2D631.root'
    #TT v2
    #'/store/mc/RunIISpring15MiniAODv2/TT_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/40000/00087FEB-236E-E511-9ACB-003048FF86CA.root'
    #TTHLep
    #v2
    '/store/mc/RunIISpring15MiniAODv2/ttHToNonbb_M125_13TeV_powheg_pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/10000/02FE2DB6-D06D-E511-8BC7-0025905C431C.root'
    #v1
    #'/store/mc/RunIISpring15DR74/ttHJetToNonbb_M125_13TeV_amcatnloFXFX_madspin_pythia8_mWCutfix/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/30000/088378DB-3D24-E511-8B0E-20CF3027A589.root'
  ),
  skipEvents = cms.untracked.uint32(0)
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )
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
#####
##   For tt+X
#####
# Setting input particle collections to be used by the tools
genJetCollection              = 'ak4GenJetsCustom'
genParticleCollection         = 'prunedGenParticles'
genJetInputParticleCollection = 'packedGenParticles'
# Supplies PDG ID to real name resolution of MC particles
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
# Producing own jets for testing purposes
from RecoJets.JetProducers.ak4GenJets_cfi import ak4GenJets
process.ak4GenJetsCustom = ak4GenJets.clone(
  src = genJetInputParticleCollection,
  rParam = cms.double(0.4),
  jetAlgorithm = cms.string("AntiKt")
)
# Ghost particle collection used for Hadron-Jet association 
# MUST use proper input particle collection
from PhysicsTools.JetMCAlgos.HadronAndPartonSelector_cfi import selectedHadronsAndPartons
process.selectedHadronsAndPartons = selectedHadronsAndPartons.clone(
  particles = genParticleCollection
)
# Input particle collection for matching to gen jets (partons + leptons) 
# MUST use use proper input jet collection: the jets to which hadrons should be associated
# rParam and jetAlgorithm MUST match those used for jets to be associated with hadrons
# More details on the tool: https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideBTagMCTools#New_jet_flavour_definition
from PhysicsTools.JetMCAlgos.sequences.GenHFHadronMatching_cff import genJetFlavourPlusLeptonInfos
process.genJetFlavourPlusLeptonInfos = genJetFlavourPlusLeptonInfos.clone(
  jets = genJetCollection,
  rParam = cms.double(0.4),
  jetAlgorithm = cms.string("AntiKt")
)
# Plugin for analysing B hadrons
# MUST use the same particle collection as in selectedHadronsAndPartons
from PhysicsTools.JetMCAlgos.sequences.GenHFHadronMatching_cff import matchGenBHadron
process.matchGenBHadron = matchGenBHadron.clone(
  genParticles = genParticleCollection
)
# Plugin for analysing C hadrons
# MUST use the same particle collection as in selectedHadronsAndPartons
from PhysicsTools.JetMCAlgos.sequences.GenHFHadronMatching_cff import matchGenCHadron
process.matchGenCHadron = matchGenCHadron.clone(
  genParticles = genParticleCollection
)
#####
##   Output file
#####
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
    #Common
     #Electron
     'HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v',
     #Muon
     'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v',
     'HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v', 
     'HLT_IsoMu20_v',
     'HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v',
     'HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v',
    #TTHbb
     #Electron
     'HLT_Ele105_CaloIdVT_GsfTrkIdT_v',
     'HLT_Ele27_eta2p1_WP75_Gsf_v',
     'HLT_Ele27_WP85_Gsf_v',
     'HLT_Ele27_eta2p1_WPLoose_Gsf_v',
     #Muon
     'HLT_Mu45_eta2p1',
     'HLT_Mu50_v',
     'HLT_IsoMu17_eta2p1_v',
     'HLT_IsoMu24_eta2p1_v',
     'HLT_IsoMu18_v',
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
  ),
  # Choose which information you want to use
  fillgeninfo           = cms.bool(True),
  fillgenHFCategoryinfo = cms.bool(True),
  filleventinfo         = cms.bool(True),
  filltriggerinfo       = cms.bool(True),
  fillPVinfo            = cms.bool(True),
  fillmuoninfo          = cms.bool(True),
  fillelectronpatinfo   = cms.bool(True),
  filltauinfo           = cms.bool(True),
  filljetinfo           = cms.bool(True),
  filltthjetinfo        = cms.bool(False), #F
  fillBoostedJetinfo    = cms.bool(True),
  fillTopSubJetinfo     = cms.bool(True),
  fillBJetnessinfo      = cms.bool(True),
  fillBJetnessFVinfo    = cms.bool(True),
  fillBTagReweight      = cms.bool(True),
  fillPileupReweight    = cms.bool(True),
  fillMETinfo           = cms.bool(True),
  fillphotoninfo        = cms.bool(False), #F   
  # Choose format 
  MiniAODv2 = cms.bool(True),
  is_data   = cms.bool(False),
  debug_    = cms.bool(False),
  super_TNT = cms.bool(False),
  AJVar     = cms.bool(False),
  tthlepVar = cms.bool(True),
  PuppiVar  = cms.bool(False),
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
  eleMVATrigIdMap     = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-Trig-V1-wp80"),
  eleHEEPIdMap        = cms.InputTag("egmGsfElectronIDs:heepElectronID-HEEPV60"),
  elemvaValuesMap_nonTrig      = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Values"),
  elemvaCategoriesMap_nonTrig  = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Categories"),
  elemvaValuesMap_Trig         = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15Trig25nsV1Values"),
  elemvaCategoriesMap_Trig     = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15Trig25nsV1Categories"),
  taus                = cms.InputTag("slimmedTaus"),
  #jets                = cms.InputTag("selectedPatJetsAK8PFCHS"),
  jets                = cms.InputTag("slimmedJets"),
  jetsPUPPI           = cms.InputTag("slimmedJetsPuppi"),
  fatjets             = cms.InputTag("slimmedJetsAK8"),
  topsubjets          = cms.InputTag("slimmedJetsCMSTopTagCHSPacked", "SubJets"),
  mets                = cms.InputTag("slimmedMETs"),
  metsPUPPI           = cms.InputTag("slimmedMETsPuppi"),
  photons             = cms.InputTag("slimmedPhotons"),
  packedPFCandidates  = cms.InputTag("packedPFCandidates"), 
  # JEC - CORRECTIONS ON FLY
  jecPayloadNamesAK4PFchsMC1   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Summer15_25nsV6_MC_L1FastJet_AK4PFchs.txt"),
  jecPayloadNamesAK4PFchsMC2   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Summer15_25nsV6_MC_L2Relative_AK4PFchs.txt"),
  jecPayloadNamesAK4PFchsMC3   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Summer15_25nsV6_MC_L3Absolute_AK4PFchs.txt"),
  jecPayloadNamesAK4PFchsMCUnc = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Summer15_25nsV6_MC_Uncertainty_AK4PFchs.txt"),
  jecPayloadNamesAK4PFchsDATA1   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer15_25nsV6_DATA_L1FastJet_AK4PFchs.txt"),
  jecPayloadNamesAK4PFchsDATA2   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer15_25nsV6_DATA_L2Relative_AK4PFchs.txt"),
  jecPayloadNamesAK4PFchsDATA3   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer15_25nsV6_DATA_L3Absolute_AK4PFchs.txt"),
  jecPayloadNamesAK4PFchsDATAUnc = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer15_25nsV6_DATA_Uncertainty_AK4PFchs.txt"),
  jecPayloadNamesAK4PFPuppiMC1   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Summer15_25nsV6_MC_L1FastJet_AK4PFPuppi.txt"),
  jecPayloadNamesAK4PFPuppiMC2   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Summer15_25nsV6_MC_L2Relative_AK4PFPuppi.txt"),
  jecPayloadNamesAK4PFPuppiMC3   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Summer15_25nsV6_MC_L3Absolute_AK4PFPuppi.txt"),
  jecPayloadNamesAK4PFPuppiMCUnc = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Summer15_25nsV6_MC_Uncertainty_AK4PFPuppi.txt"),
  jecPayloadNamesAK4PFPuppiDATA1   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer15_25nsV6_DATA_L1FastJet_AK4PFPuppi.txt"),
  jecPayloadNamesAK4PFPuppiDATA2   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer15_25nsV6_DATA_L2Relative_AK4PFPuppi.txt"),
  jecPayloadNamesAK4PFPuppiDATA3   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer15_25nsV6_DATA_L3Absolute_AK4PFPuppi.txt"),
  jecPayloadNamesAK4PFPuppiDATAUnc = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer15_25nsV6_DATA_Uncertainty_AK4PFPuppi.txt"),
  jecPayloadNamesAK8PFchsMC1   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Summer15_25nsV6_MC_L1FastJet_AK8PFchs.txt"),
  jecPayloadNamesAK8PFchsMC2   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Summer15_25nsV6_MC_L2Relative_AK8PFchs.txt"),
  jecPayloadNamesAK8PFchsMC3   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Summer15_25nsV6_MC_L3Absolute_AK8PFchs.txt"),
  jecPayloadNamesAK8PFchsMCUnc = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Summer15_25nsV6_MC_Uncertainty_AK8PFchs.txt"),
  jecPayloadNamesAK8PFchsDATA1   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer15_25nsV6_DATA_L1FastJet_AK8PFchs.txt"),
  jecPayloadNamesAK8PFchsDATA2   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer15_25nsV6_DATA_L2Relative_AK8PFchs.txt"),
  jecPayloadNamesAK8PFchsDATA3   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer15_25nsV6_DATA_L3Absolute_AK8PFchs.txt"),
  jecPayloadNamesAK8PFchsDATAUnc = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer15_25nsV6_DATA_Uncertainty_AK8PFchs.txt"),
  # PILEUP REWEIGHTING
  PUReweightfile      = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/PUReweight/MyDataPileupHistogram_true.root"),
  # BTAG REWEIGHTING
  BTAGReweightfile1   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/BTAGReweight/csv_rwt_hf_IT_FlatSF_2015_07_27.root"),
  BTAGReweightfile2   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/BTAGReweight/csv_rwt_lf_IT_FlatSF_2015_07_27.root"),
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
  patElectron_pt_min  = cms.double(7.),
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
  genJetPtMin               = cms.double(20),
  genJetAbsEtaMax           = cms.double(2.4),
  genJets                   = cms.InputTag("ak4GenJetsCustom"),
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
#QG likelihood
process.load('BSMFramework.BSM3G_TNT_Maker.QGTagger_cfi')
process.QGTagger.srcJets          = cms.InputTag('slimmedJets')
process.QGTagger.jetsLabel        = cms.string('QGL_AK4PFchs')
#process.p = cms.Path(process.printGenParticleList)
process.p = cms.Path(
process.selectedHadronsAndPartons*process.ak4GenJetsCustom*process.genJetFlavourPlusLeptonInfos*process.matchGenCHadron*process.selectedHadronsAndPartons*process.genJetFlavourPlusLeptonInfos*process.matchGenBHadron*
process.egmGsfElectronIDSequence*
process.QGTagger*
#process.primaryVertexFilter* 
#process.CSCTightHaloFilter*process.eeBadScFilter*process.HBHENoiseFilterResultProducer*process.ApplyBaselineHBHENoiseFilter*
process.TNT
)
