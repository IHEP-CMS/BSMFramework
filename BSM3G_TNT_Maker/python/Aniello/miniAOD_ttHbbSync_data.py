import FWCore.ParameterSet.Config as cms
#####
##   Initial standard configs
#####
process = cms.Process("Demo")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 500
#process.load("Configuration.StandardSequences.Geometry_cff")
##process.load('Configuration.Geometry.GeometryIdeal_cff')
process.load('Configuration.Geometry.GeometryRecoDB_cff')
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
#process.load("Configuration.StandardSequences.MagneticField_38T_PostLS1_cff")
#process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
#process.GlobalTag.globaltag = 'MCRUN2_74_V9::All'
process.GlobalTag.globaltag = '74X_dataRun2_v5'
process.prefer("GlobalTag")
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
#####
##   Input files
#####
process.source = cms.Source("PoolSource",
  fileNames = cms.untracked.vstring(
    #data
    'root://xrootd-cms.infn.it///store/user/hmildner/el_skim_3loosejets.root',
    'root://xrootd-cms.infn.it///store/user/hmildner/mu_skim_3loosejets.root',
    'root://xrootd-cms.infn.it///store/user/hmildner/muel_skim_2loosejets.root'
  ),
  skipEvents = cms.untracked.uint32(0)
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
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
# Do not forget to add the egmGsfElectronIDSequence to the path,
# as in the example below!
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
  fileName = cms.string("OutTree_data.root")
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
  fillgeninfo           = cms.bool(False),
  fillgenHFCategoryinfo = cms.bool(False),
  filleventinfo         = cms.bool(True),
  filltriggerinfo       = cms.bool(True),
  fillPVinfo            = cms.bool(False),
  fillmuoninfo          = cms.bool(True),
  fillelectronpatinfo   = cms.bool(True),
  filltauinfo           = cms.bool(False),
  filljetinfo           = cms.bool(True),
  fillBoostedJetinfo    = cms.bool(False),
  fillTopSubJetinfo     = cms.bool(False),
  fillBJetnessinfo      = cms.bool(False),
  fillBTagReweight      = cms.bool(True),
  fillPileupReweight    = cms.bool(False),
  fillMETinfo           = cms.bool(True),
  fillphotoninfo        = cms.bool(False), 
  filltthjetinfo        = cms.bool(False),  
  # Choose format
  MiniAODv2 = cms.bool(True),
  is_data   = cms.bool(True),
  debug_    = cms.bool(False),
  super_TNT = cms.bool(False),
  AJVar     = cms.bool(False),
  tthlepVar = cms.bool(False),
  PuppiVar  = cms.bool(False),
  ##### Input tags 
  triggerResults      = cms.InputTag('TriggerResults', '', 'HLT' ),
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
  #JEC - CORRECTIONS ON FLY
  jecPayloadNamesAK4PFchsMC1   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/MCRUN2_74_V9_L1FastJet_AK4PFchs.txt"),
  jecPayloadNamesAK4PFchsMC2   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/MCRUN2_74_V9_L2Relative_AK4PFchs.txt"),
  jecPayloadNamesAK4PFchsMC3   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/MCRUN2_74_V9_L3Absolute_AK4PFchs.txt"),
  jecPayloadNamesAK4PFchsMCUnc = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/MCRUN2_74_V9_Uncertainty_AK4PFchs.txt"),
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
  #PILEUP REWEIGHTING
  PUReweightfile      = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/PUReweight/MyDataPileupHistogram_true.root"),
  #BTAG REWEIGHTING
  BTAGReweightfile1   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/BTAGReweight/csv_rwt_fit_hf_2015_11_20.root"),
  BTAGReweightfile2   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/BTAGReweight/csv_rwt_fit_lf_2015_11_20.root"),
  #### Object selection
  # Primary vertex cuts
  Pvtx_ndof_min   = cms.double(4.),
  Pvtx_vtx_max    = cms.double(24.),
  Pvtx_vtxdxy_max = cms.double(24.),
  # Muon cuts
  Muon_pt_min         = cms.double(5.),
  Muon_eta_max        = cms.double(3.0),
  vtx_ndof_min        = cms.int32(4),
  vtx_rho_max         = cms.int32(2),
  vtx_position_z_max  = cms.double(24.),
  # Electron cuts
  patElectron_pt_min  = cms.double(5.),
  patElectron_eta_max = cms.double(3.0),
  # Tau cuts
  Tau_pt_min              = cms.double(10.),
  Tau_eta_max             = cms.double(3.0),
  Tau_vtx_ndof_min        = cms.int32(4),
  Tau_vtx_rho_max         = cms.int32(2),
  Tau_vtx_position_z_max  = cms.double(24.),
  # Jet cuts
  Jet_pt_min = cms.double(5.),
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
#process.p = cms.Path(process.printGenParticleList)
process.p = cms.Path(
#process.hltFilter*
#process.selectedHadronsAndPartons*process.ak4GenJetsCustom*process.genJetFlavourPlusLeptonInfos*process.matchGenCHadron*process.selectedHadronsAndPartons*process.genJetFlavourPlusLeptonInfos*process.matchGenBHadron*
process.egmGsfElectronIDSequence*
#process.primaryVertexFilter* 
#process.CSCTightHaloFilter*process.eeBadScFilter*process.HBHENoiseFilterResultProducer*process.ApplyBaselineHBHENoiseFilter*
process.TNT)
