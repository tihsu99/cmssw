import os

import FWCore.ParameterSet.Config as cms
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
from Configuration.AlCa.GlobalTag import GlobalTag

from Configuration.StandardSequences.Eras import eras
process = cms.Process("ElectronMVANtuplizer", eras.Phase2C17I13M9)


process.load("Configuration.Geometry.GeometryExtended2026D88_cff")
process.load("Configuration.Geometry.GeometryExtended2026D88Reco_cff")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

process.GlobalTag = GlobalTag(process.GlobalTag, "auto:phase2_realistic_T21", "")

import FWCore.ParameterSet.VarParsing as VarParsing
options = VarParsing.VarParsing("analysis")

options.register("sourceFile",
    "", # Default value
    VarParsing.VarParsing.multiplicity.singleton, # singleton or list
    VarParsing.VarParsing.varType.string, # string, int, or float
    "File containing list of input files" # Description
)

options.register("mvaVariablesFile",
    "RecoEgamma/ElectronIdentification/data/ElectronIDVariables.txt", # Default value
    VarParsing.VarParsing.multiplicity.singleton, # singleton or list
    VarParsing.VarParsing.varType.string, # string, int, or float
    "MVA variables file" # Description
)

options.register("electronLabel",
    "gedGsfElectrons", # Default value
    VarParsing.VarParsing.multiplicity.singleton, # singleton or list
    VarParsing.VarParsing.varType.string, # string, int, or float
    "Electron collection name (e.g. gedGsfElectrons, ecalDrivenGsfElectronsHGC)" # Description
)

options.register("outDir",
    "", # Default value
    VarParsing.VarParsing.multiplicity.singleton, # singleton or list
    VarParsing.VarParsing.varType.string, # string, int, or float
    "Output directory" # Description
)

options.register("outFileNumber",
    -1, # Default value
    VarParsing.VarParsing.multiplicity.singleton, # singleton or list
    VarParsing.VarParsing.varType.int, # string, int, or float
    "File number (will be added to the filename if >= 0)" # Description
)

options.register("useAOD",
    0, # Default value
    VarParsing.VarParsing.multiplicity.singleton, # singleton or list
    VarParsing.VarParsing.varType.int, # string, int, or float
    "Use AOD" # Description
)

options.register("debugFile",
    0, # Default value
    VarParsing.VarParsing.multiplicity.singleton, # singleton or list
    VarParsing.VarParsing.varType.int, # string, int, or float
    "Create debug file" # Description
)

options.parseArguments()

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )


# Input
if (len(options.sourceFile)) :
    
    sourceFile = options.sourceFile

fNames = []

if (len(options.inputFiles)) :
    
    fNames = options.inputFiles

else :
    
    with open(sourceFile) as f:
        
        fNames = f.readlines()

fNames = [f for f in fNames if f[0] != "#"]

for iFile, fName in enumerate(fNames) :
    
    if (
        "file:" not in fName and
        "root:" not in fName
    ) :
        
        fNames[iFile] = "file:%s" %(fName)

sourceFileNames = cms.untracked.vstring(fNames)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        sourceFileNames
    ),
    #inputCommands=cms.untracked.vstring(
    #    "drop l1tTkPrimaryVertexs_L1TkPrimaryVertex__HLT",
    #)
)

# Output
outFileSuffix = ""

if (options.outFileNumber >= 0) :
    
    outFileSuffix = "%s_%d" %(outFileSuffix, options.outFileNumber)

outFile = "ntupleTree%s.root" %(outFileSuffix)

if (len(options.outDir)) :
    
    os.system("mkdir -p %s" %(options.outDir))
    outFile = "%s/%s" %(options.outDir, outFile)

process.TFileService = cms.Service("TFileService", fileName = cms.string(outFile))


useAOD = bool(options.useAOD)

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
# turn on VID producer, indicate data format  to be
# DataFormat.AOD or DataFormat.MiniAOD, as appropriate

electronTag = cms.InputTag(options.electronLabel)

if useAOD == True :
    dataFormat = DataFormat.AOD
    input_tags = dict(
        src = electronTag,
        vertices = cms.InputTag("offlinePrimaryVertices"),
        vertices4d = cms.InputTag("offlinePrimaryVertices4D"),
        pileup = cms.InputTag("addPileupInfo"),
        genParticles = cms.InputTag("genParticles"),
        ebReducedRecHitCollection = cms.InputTag("reducedEcalRecHitsEB"),
        eeReducedRecHitCollection = cms.InputTag("reducedEcalRecHitsEE"),
    )
else :
    
    dataFormat = DataFormat.MiniAOD
    input_tags = dict(
        src = electronTag,
    )


#hltEgammaCandidatesUnseeded = cms.EDProducer("EgammaHLTRecoEcalCandidateProducers",
#    recoEcalCandidateCollection = cms.string(''),
#    scHybridBarrelProducer = cms.InputTag("particleFlowSuperClusterECAL", "particleFlowBasicClusterECALBarrel"),
#    scIslandEndcapProducer = cms.InputTag("particleFlowSuperClusterHGCal")
#)



# https://cmssdt.cern.ch/lxr/source/RecoEgamma/EgammaHLTProducers/plugins/EgammaHLTElectronTrackIsolationProducers.cc
# https://gitlab.cern.ch/sharper/EgHLTPhase2/-/blob/master/EDProducers/hltEgammaEleGsfTrackIsoUnseeded_cfi.py
#process.hltEgammaEleGsfTrackIso = cms.EDProducer("EgammaHLTElectronTrackIsolationProducers",
#    beamSpotProducer = cms.InputTag("offlineBeamSpot"),
#    egTrkIsoConeSize = cms.double(0.3),
#    egTrkIsoPtMin = cms.double(1.0),
#    egTrkIsoRSpan = cms.double(999999.0),
#    egTrkIsoStripBarrel = cms.double(0.01),
#    egTrkIsoStripEndcap = cms.double(0.01),
#    egTrkIsoVetoConeSizeBarrel = cms.double(0.01),
#    egTrkIsoVetoConeSizeEndcap = cms.double(0.01),
#    egTrkIsoZSpan = cms.double(0.15),
#    electronProducer = cms.InputTag("slimmedElectrons"),
#    #recoEcalCandidateProducer = cms.InputTag("hltEgammaCandidatesUnseeded"),
#    #trackProducer = cms.InputTag("generalTracks"),
#    trackProducer = cms.InputTag("trackExtenderWithMTD"),
#    useGsfTrack = cms.bool(True),
#    #useSCRefs = cms.bool(True)
#    useSCRefs = cms.bool(False)
#)

# https://cmssdt.cern.ch/lxr/source/RecoEgamma/EgammaIsolationAlgos/python/electronTrackIsolationScone_cfi.py#0003
electronTrackIsolation = cms.EDProducer("EgammaElectronTkIsolationProducerNew",
    absolut = cms.bool(True),
    trackProducer = cms.InputTag("generalTracks"),
    #trackProducer = cms.InputTag("trackExtenderWithMTD"),
    mtdt0 = cms.InputTag("tofPID:t0"),
    mtdSigmat0 = cms.InputTag("tofPID:sigmat0"),
    mtdTrkQualMVA = cms.InputTag("mtdTrackQualityMVA:mtdQualMVA"),
    intRadiusBarrel = cms.double(0.01),
    intRadiusEndcap = cms.double(0.01),
    stripBarrel = cms.double(0.01),
    stripEndcap = cms.double(0.01),
    electronProducer = electronTag,
    extRadius = cms.double(0.3),
    ptMin = cms.double(1.0),
    maxVtxDist = cms.double(0.15),
    BeamspotProducer = cms.InputTag("offlineBeamSpot"),
    maxVtxDistXY     = cms.double(9999.0),
    vertexProducer = cms.InputTag("offlineSlimmedPrimaryVertices4D"),
    dtRef = cms.int32(0),
    dtType = cms.int32(0),
    dtMax = cms.double(9999.0),
    trkMtdMvaMin = cms.double(0.5),
)

#process.electronTrackIsolationMod1 = process.electronTrackIsolation.clone(
#    dtRef = cms.int32(0),
#    dtType = cms.int32(0),
#    dtMax = cms.double(3.0),
#)

#l_dtMax = [1, 2, 3, 5, 7, 10, 15, 9999]
l_dtMax = [0.01, 9999]
l_isoProdLabel = []
extraTask = cms.Task()

for iVal, dtMax in enumerate(l_dtMax) :
    
    dtMaxStr = "DtMax%s" %(str(dtMax).replace(".", "p"))
    
    prodLabel = "eleTrkIso%s" %(dtMaxStr)
    
    setattr(
        process,
        prodLabel,
        electronTrackIsolation.clone(
            dtRef = cms.int32(1),
            dtType = cms.int32(0),
            dtMax = cms.double(dtMax),
        )
    )
    
    l_isoProdLabel.append(prodLabel)
    extraTask.add(getattr(process, prodLabel))


process.ntuplizer = cms.EDAnalyzer("ElectronMVANtuplizer",
    variableDefinition   = cms.string(options.mvaVariablesFile),
    ptThreshold = cms.double(5.0),
    #
    doEnergyMatrix = cms.bool(False), # disabled by default due to large size
    energyMatrixSize = cms.int32(2), # corresponding to 5x5\
    #
    eleMVAValMapLabels = cms.vstring(
        l_isoProdLabel
        #"electronTrackIsolation",
        #"electronTrackIsolationMod1",
    ),
    eleMVAValMaps = cms.vstring(
        l_isoProdLabel
        #"electronTrackIsolation",
        #"electronTrackIsolationMod1",
    ),
    
    #
    **input_tags
)

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.RecoSim_cff')
#proces.lo
# RecoParticleFlow/PFClusterProducer/python/particleFlowRecHitHGC_cfi

#process.load('CalibTracker.SiStripLorentzAngle.SiStripLorentzAngle_cfi')
process.load('CalibTracker.SiStripESProducers.SiStripLorentzAngleDepESProducer_cfi')

#from RecoEgamma.EgammaElectronProducers.ecalDrivenGsfElectronCores_cfi import ecalDrivenGsfElectronCores
#from RecoEgamma.EgammaElectronProducers.ecalDrivenGsfElectronCoresHGC_cff import ecalDrivenGsfElectronCoresHGC
#from RecoEgamma.EgammaElectronProducers.gsfElectrons_cfi import *
#
#process.ecalDrivenGsfElectronsHGCmod = ecalDrivenGsfElectronsHGC.clone()
#process.ecalDrivenGsfElectronCoresHGCmod = ecalDrivenGsfElectronCoresHGC.clone()
#gsfEcalDrivenElectronTaskHGCmod = cms.Task(process.ecalDrivenGsfElectronCoresHGCmod, process.ecalDrivenGsfElectronsHGCmod)


#from RecoParticleFlow.PFClusterProducer.particleFlowRecHitHGC_cfi import *
#from RecoEcal.EgammaClusterProducers.particleFlowSuperClusteringSequence_cff import *
#from RecoEgamma.EgammaElectronProducers.ecalDrivenElectronSeeds_cff import *
#from RecoEgamma.EgammaElectronProducers.ecalDrivenElectronSeeds_cfi import *
#from RecoEgamma.EgammaElectronProducers.ecalDrivenGsfElectronCoresHGC_cff import *
#from RecoEgamma.EgammaElectronProducers.gsfElectrons_cfi import *
#from RecoParticleFlow.PFClusterProducer.particleFlowClusterHGC_cfi import *
#from RecoParticleFlow.PFTracking.mergedElectronSeeds_cfi import *
#from RecoTracker.IterativeTracking.ElectronSeeds_cff import *
#from TrackingTools.GsfTracking.CkfElectronCandidateMaker_cff import *
#from TrackingTools.GsfTracking.GsfElectronGsfFit_cff import *
#
#process.particleFlowClusterHGCalmod         = particleFlowClusterHGCal.clone()
#process.particleFlowSuperClusterHGCalmod    = particleFlowSuperClusterHGCal.clone()
#process.ecalDrivenElectronSeedsmod          = ecalDrivenElectronSeeds.clone()
#process.electronMergedSeedsmod              = electronMergedSeeds.clone()
#process.electronCkfTrackCandidatesmod       = electronCkfTrackCandidates.clone()
#process.electronGsfTracksmod                = electronGsfTracks.clone()
#process.ecalDrivenGsfElectronCoresHGCmod    = ecalDrivenGsfElectronCoresHGC.clone()
#process.ecalDrivenGsfElectronsHGCmod        = ecalDrivenGsfElectronsHGC.clone()
#
#
#process.particleFlowClusterHGCalmod.initialClusteringStep.filterByTracksterPID = cms.bool(True)
#process.particleFlowClusterHGCalmod.initialClusteringStep.filterByTracksterIteration = cms.bool(False)
##process.particleFlowClusterHGCalmod.recHitsSource = cms.InputTag("particleFlowRecHitHGC", "Cleaned", "reRECO")
#process.particleFlowClusterHGCalmod.recHitsSource = cms.InputTag("particleFlowRecHitHGC", "", "ElectronMVANtuplizer")
##process.particleFlowSuperClusterHGCalmod.ESAssociation = cms.InputTag("")
#process.particleFlowSuperClusterHGCalmod.PFClusters = cms.InputTag("particleFlowClusterHGCalmod")
#process.ecalDrivenElectronSeedsmod.endcapSuperClusters = cms.InputTag("particleFlowSuperClusterHGCalmod")
#process.electronMergedSeedsmod.EcalBasedSeeds = cms.InputTag("ecalDrivenElectronSeedsmod")
#process.electronCkfTrackCandidatesmod.src = "electronMergedSeedsmod"
#process.electronGsfTracksmod.src = "electronCkfTrackCandidatesmod"
#process.ecalDrivenGsfElectronCoresHGCmod.gsfTracks = "electronGsfTracksmod"
#process.ecalDrivenGsfElectronsHGCmod.gsfElectronCoresTag = "ecalDrivenGsfElectronCoresHGCmod"
#
#
#process.ecalDrivenGsfElectronsmod_step = cms.Sequence(
#    process.particleFlowRecHitHGC *
#    process.particleFlowClusterHGCalmod *
#    process.particleFlowSuperClusterHGCalmod *
#    process.ecalDrivenElectronSeedsmod *
#    process.electronMergedSeedsmod *
#    process.electronCkfTrackCandidatesmod *
#    process.electronGsfTracksmod *
#    process.ecalDrivenGsfElectronCoresHGCmod *
#    process.ecalDrivenGsfElectronsHGCmod
#)

#gsfEcalDrivenElectronTaskHGCmod = cms.Task(
#    process.siStripLorentzAngleDepESProducer,
#    process.particleFlowClusterHGCalmod,
#    process.particleFlowSuperClusterHGCalmod,
#    process.ecalDrivenElectronSeedsmod,
#    process.electronMergedSeedsmod,
#    process.electronCkfTrackCandidatesmod,
#    process.electronGsfTracksmod,
#    process.ecalDrivenGsfElectronCoresHGCmod,
#    process.ecalDrivenGsfElectronsHGCmod
#)
#
#process.p.associate(gsfEcalDrivenElectronTaskHGCmod)


process.p = cms.Path(
    #process.egmGsfElectronIDSequence *
    #process.hltEgammaEleGsfTrackIso *
    
    #process.electronTrackIsolation *
    #process.electronTrackIsolationMod1 *
    
    #process.ecalDrivenGsfElectronsmod_step *
    process.ntuplizer
)

process.p.associate(extraTask)


process.schedule = cms.Schedule(process.p)



# Debug
if (options.debugFile) :
    
    process.out = cms.OutputModule("PoolOutputModule",
        fileName = cms.untracked.string("debug.root")
    )
    
    process.output_step = cms.EndPath(process.out)
    process.schedule.extend([process.output_step])
