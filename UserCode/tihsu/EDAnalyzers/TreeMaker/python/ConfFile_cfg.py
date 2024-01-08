import os

import FWCore.ParameterSet.Config as cms

processName = "Demo"

from Configuration.StandardSequences.Eras import eras

d_conditionOpt = {}

d_conditionOpt["Phase2HLTTDRSummer20ReRECOMiniAOD"] = {
    "GT": "auto:phase2_realistic_T15",
    "Geom": "GeometryExtended2026D49",
    "Era": eras.Phase2C9,
}

d_conditionOpt["PhaseIISpring22DRMiniAOD"] = {
    "GT": "auto:phase2_realistic_T21",
    #"GT": "123X_mcRun4_realistic_v11",
    "Geom": "GeometryExtended2026D88",
    "Era": eras.Phase2C17I13M9,
}

############################## Parse arguments ##############################

import FWCore.ParameterSet.VarParsing as VarParsing
options = VarParsing.VarParsing("analysis")

options.register("sourceFile",
    "", # Default value
    VarParsing.VarParsing.multiplicity.singleton, # singleton or list
    VarParsing.VarParsing.varType.string, # string, int, or float
    "File containing list of input files" # Description
)

options.register("outputDir",
    "", # Default value
    VarParsing.VarParsing.multiplicity.singleton, # singleton or list
    VarParsing.VarParsing.varType.string, # string, int, or float
    "Output directory" # Description
)

options.register("outFileBaseName",
    "ntupleTree", # Default value
    VarParsing.VarParsing.multiplicity.singleton, # singleton or list
    VarParsing.VarParsing.varType.string, # string, int, or float
    "Output file base name. Number and extenstion will be added automatically." # Description
)

options.register("outFileNumber",
    -1, # Default value
    VarParsing.VarParsing.multiplicity.singleton, # singleton or list
    VarParsing.VarParsing.varType.int, # string, int, or float
    "File number (will be added to the filename if >= 0)" # Description
)

options.register("eventRange",
    [], # Default value
    VarParsing.VarParsing.multiplicity.list, # singleton or list
    VarParsing.VarParsing.varType.string, # string, int, or float
    "Syntax: Run1:Event1-Run2:Event2 Run3:Event3-Run4:Event4(includes both)" # Description
)

options.register("debugFile",
    0, # Default value
    VarParsing.VarParsing.multiplicity.singleton, # singleton or list
    VarParsing.VarParsing.varType.int, # string, int, or float
    "Create debug file" # Description
)

options.register("onRaw",
    0, # Default value
    VarParsing.VarParsing.multiplicity.singleton, # singleton or list
    VarParsing.VarParsing.varType.int, # string, int, or float
    "Running on RAW" # Description
)

options.register("conditionOpt",
    None, # Default value
    VarParsing.VarParsing.multiplicity.singleton, # singleton or list
    VarParsing.VarParsing.varType.string, # string, int, or float
    "Condition option. Choices: %s" %(", ".join(list(d_conditionOpt.keys()))) # Description
)

options.register("storeSimHit",
    0, # Default value
    VarParsing.VarParsing.multiplicity.singleton, # singleton or list
    VarParsing.VarParsing.varType.int, # string, int, or float
    "Store sim-hits" # Description
)

options.register("storeRecHit",
    0, # Default value
    VarParsing.VarParsing.multiplicity.singleton, # singleton or list
    VarParsing.VarParsing.varType.int, # string, int, or float
    "Store rec-hits" # Description
)

options.register("storeHGCALlayerClus",
    0, # Default value
    VarParsing.VarParsing.multiplicity.singleton, # singleton or list
    VarParsing.VarParsing.varType.int, # string, int, or float
    "Store HGCal layer clusters" # Description
)

options.register("storeSuperClusTICLclus",
    0, # Default value
    VarParsing.VarParsing.multiplicity.singleton, # singleton or list
    VarParsing.VarParsing.varType.int, # string, int, or float
    "Store info about TICL-electron SC, SC seed, and TICL-cluster matches" # Description
)

options.register("TICLeleGenMatchDR",
    99999, # Default value
    VarParsing.VarParsing.multiplicity.singleton, # singleton or list
    VarParsing.VarParsing.varType.float, # string, int, or float
    "DeltaR to use for TICL-electron gen-matching (will store only the gen-matched ones)" # Description
)

options.register("TICLphoGenMatchDR",
    99999, # Default value
    VarParsing.VarParsing.multiplicity.singleton, # singleton or list
    VarParsing.VarParsing.varType.float, # string, int, or float
    "DeltaR to use for TICL-photon gen-matching (will store only the gen-matched ones)" # Description
)

options.register("isGunSample",
    0, # Default value
    VarParsing.VarParsing.multiplicity.singleton, # singleton or list
    VarParsing.VarParsing.varType.int, # string, int, or float
    "Is it a particle gun sample" # Description
)

options.register("genEleFilter",
    0, # Default value
    VarParsing.VarParsing.multiplicity.singleton, # singleton or list
    VarParsing.VarParsing.varType.int, # string, int, or float
    "Apply gen-electron filter" # Description
)

options.register("genPhoFilter",
    0, # Default value
    VarParsing.VarParsing.multiplicity.singleton, # singleton or list
    VarParsing.VarParsing.varType.int, # string, int, or float
    "Apply gen-photon filter" # Description
)

options.register("genPartonFilter",
    0, # Default value
    VarParsing.VarParsing.multiplicity.singleton, # singleton or list
    VarParsing.VarParsing.varType.int, # string, int, or float
    "Apply gen-parton filter" # Description
)

options.register("trace",
    0, # Default value
    VarParsing.VarParsing.multiplicity.singleton, # singleton or list
    VarParsing.VarParsing.varType.int, # string, int, or float
    "Trace modules" # Description
)

options.register("memoryCheck",
    0, # Default value
    VarParsing.VarParsing.multiplicity.singleton, # singleton or list
    VarParsing.VarParsing.varType.int, # string, int, or float
    "Check memory usage" # Description
)

options.register("printTime",
    0, # Default value
    VarParsing.VarParsing.multiplicity.singleton, # singleton or list
    VarParsing.VarParsing.varType.int, # string, int, or float
    "Print timing information" # Description
)

options.register("depGraph",
    0, # Default value
    VarParsing.VarParsing.multiplicity.singleton, # singleton or list
    VarParsing.VarParsing.varType.int, # string, int, or float
    "Produce dependency graph only" # Description
)

options.parseArguments()

assert(options.conditionOpt in d_conditionOpt.keys())

process = cms.Process(processName, d_conditionOpt[options.conditionOpt]["Era"])

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
##process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
##process.load('SimGeneral.MixingModule.mixNoPU_cfi')
##process.load('Configuration.Geometry.GeometryExtended2023D41Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.RecoSim_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
##process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, d_conditionOpt[options.conditionOpt]["GT"], "")
process.load('Configuration.Geometry.%sReco_cff' %(d_conditionOpt[options.conditionOpt]["Geom"]))
process.load('Configuration.Geometry.%s_cff' %(d_conditionOpt[options.conditionOpt]["Geom"]))

process.load("Geometry.HGCalGeometry.HGCalGeometryESProducer_cfi")


#maxEvents = -1
#options.maxEvents = 15
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(options.maxEvents))


#sourceFile = "sourceFiles/SingleElectronFlatPtGun_pT-15_rajdeep/SingleElectronFlatPtGun_pT-15_rajdeep.txt"
#sourceFile = "sourceFiles/SingleElectronFlatPtGun_pT-15_eta-1p5-3p0_sobhatta-crab_SingleElectronFlatPtGun_pT-15_eta-1p5-3p0_GEN-SIM-RECO-ffc2278112c688bef3890fc698a39794_USER/SingleElectronFlatPtGun_pT-15_eta-1p5-3p0_sobhatta-crab_SingleElectronFlatPtGun_pT-15_eta-1p5-3p0_GEN-SIM-RECO-ffc2278112c688bef3890fc698a39794_USER.txt"
#sourceFile = "sourceFiles/SingleElectronFlatPtGun_pT-15_eta-1p5-3p0_withTICLfractions/SingleElectronFlatPtGun_pT-15_eta-1p5-3p0_withTICLfractions.txt"

#sourceFile = "sourceFiles/SinglePi0FlatPtGun_pT-15_eta-1p5-3p0_GEN-SIM-RECO/SinglePi0FlatPtGun_pT-15_eta-1p5-3p0_GEN-SIM-RECO.txt"
#sourceFile = "sourceFiles/SingleZprimeToEEflatPtGun_m-5_pT-35_eta-1p5-3p0_GEN-SIM-RECO/SingleZprimeToEEflatPtGun_m-5_pT-35_eta-1p5-3p0_GEN-SIM-RECO.txt"

#sourceFile = "sourceFiles/RelValElectronGunPt2To100_CMSSW_10_6_0_pre4-106X_upgrade2023_realistic_v2_2023D41noPU-v1_GEN-SIM-RECO/RelValElectronGunPt2To100_CMSSW_10_6_0_pre4-106X_upgrade2023_realistic_v2_2023D41noPU-v1_GEN-SIM-RECO.txt"
#sourceFile = "sourceFiles/RelValElectronGunPt2To100_CMSSW_10_6_0_patch2-106X_upgrade2023_realistic_v3_2023D41noPU-v1_GEN-SIM-RECO/RelValElectronGunPt2To100_CMSSW_10_6_0_patch2-106X_upgrade2023_realistic_v3_2023D41noPU-v1_GEN-SIM-RECO.txt"

sourceFile = "sourceFiles/SingleElectronFlatPtGun_fpantale_pT-0-200_eta-1p5-3p0_GEN-SIM-RECO/SingleElectronFlatPtGun_fpantale_pT-0-200_eta-1p5-3p0_GEN-SIM-RECO_mod.txt"

if (len(options.sourceFile)) :
    
    sourceFile = options.sourceFile


fNames = []

if (len(options.inputFiles)) :
    
    fNames = options.inputFiles

else :
    
    with open(sourceFile) as f:
        
        fNames = f.readlines()

#fNames = ["file:/afs/cern.ch/work/s/sobhatta/private/HGCal_ele-reco/CMSSW_11_0_0_pre4/src/output_GEN-SIM-RECO_numEvent20.root"]
#fNames = ["file:/afs/cern.ch/work/s/sobhatta/private/HGCal_ele-reco/CMSSW_11_0_0_pre4/src/output_GEN-SIM-RECO_numEvent100.root"]
#fNames = ["file:/afs/cern.ch/work/s/sobhatta/private/HGCal_ele-reco/CMSSW_11_0_0_pre4/src/output_GEN-SIM-RECO_numEvent2000.root"]
#fNames = ["file:/eos/cms/store/group/phys_egamma/fpantale/output_GEN-SIM-RECO.root"]


for iFile, fName in enumerate(fNames) :
    
    if (
        "file:" not in fName and
        "root:" not in fName
    ) :
        
        fNames[iFile] = "file:%s" %(fName)


outFileSuffix = ""


if (options.onRaw) :
    
    outFileSuffix = "%s_onRaw" %(outFileSuffix)


if (options.outFileNumber >= 0) :
    
    outFileSuffix = "%s_%d" %(outFileSuffix, options.outFileNumber)


outFile = "%s%s.root" %(options.outFileBaseName, outFileSuffix)

if (len(options.outputDir)) :
    
    os.system("mkdir -p %s" %(options.outputDir))
    
    outFile = "%s/%s" %(options.outputDir, outFile)


sourceFileNames = cms.untracked.vstring(fNames)
#print sourceFileNames

process.source = cms.Source("PoolSource",
    fileNames = sourceFileNames,
    
    # Run1:Event1 to Run2:Event2
    #eventsToProcess = cms.untracked.VEventRange("1:78722-1:78722"),
    
    #duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
)


if (len(options.eventRange)) :
    
    process.source.eventsToProcess = cms.untracked.VEventRange(options.eventRange)


if (options.depGraph) :
    
    process.DependencyGraph = cms.Service("DependencyGraph")
    process.source = cms.Source("EmptySource")
    process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(0))


inputProcessName = ""

if (options.onRaw) :
    
    process.reconstruction_mod = process.reconstruction.copy()
    inputProcessName = processName


label_gsfEleFromTICL = cms.InputTag("ecalDrivenGsfElectronsHGC", "", inputProcessName)
label_phoFromTICL = cms.InputTag("photonsHGC", "", inputProcessName)

label_TICLtrackster = cms.InputTag("ticlTrackstersMerge", "", inputProcessName)
label_TICLmultiCluster = cms.InputTag("particleFlowClusterHGCal", "", inputProcessName)

label_hgcalLayerClus = cms.InputTag("hgcalLayerClusters",  "", inputProcessName),

# TICL-ele variables
from MyTools.EDProducers.producers_cfi import *

l_var_TICLele = []
l_mapProdVars_TICLele = []
TICLeleVar_task = cms.Task()

l_var_TICLpho = []
l_mapProdVars_TICLpho = []
TICLphoVar_task = cms.Task()

#l_rad = [1.8, 2.0, 2.8, 3.5]
l_rad = [2.8]
l_hitEnCut = [0.0]#, 0.02, 0.05, 0.07, 0.1, 0.5, 1.0, 1.5, 2.0, 3.0]

#l_isoCone = [0.12, 0.15, 0.2, 0.25]
l_isoCone = [0.15]
l_clusEnCut = [0.0]#, 0.02, 0.05, 0.07, 0.1, 0.5, 1.0, 1.5, 2.0, 3.0]


for iCone, cone in enumerate(l_isoCone) :
    
    for iEnCut, enCut in enumerate(l_clusEnCut) :
        
        radStr = "R%s" %(str(cone).replace(".", "p"))
        enCutStr = "En%s" %(str(enCut).replace(".", "p"))
        
        
        # Ele H/E
        prodLabel = "TICLeleHoverEProducer%s%s" %(radStr, enCutStr)
        
        setattr(
            process,
            prodLabel,
            HGCalElectronHoverEProducer.clone(
                electrons = label_gsfEleFromTICL,
                #layerClusters = label_hgcalLayerClus,
                coneDR = cone,
                minClusE = enCut,
                minPt = 10.0,
                #debug = cms.bool(True),
            )
        )
        
        TICLeleVar_task.add(getattr(process, prodLabel))
        
        l_mapProdVars_TICLele.append(cms.InputTag(prodLabel, getattr(process, prodLabel).instanceName.value(), processName))
        l_var_TICLele.append("%s_%s" %(prodLabel, getattr(process, prodLabel).instanceName.value()))
        
        
        # Pho H/E
        prodLabel = "TICLphoHoverEProducer%s%s" %(radStr, enCutStr)
        
        setattr(
            process,
            prodLabel,
            HGCalPhotonHoverEProducer.clone(
                photons = label_phoFromTICL,
                #layerClusters = label_hgcalLayerClus,
                coneDR = cone,
                minClusE = enCut,
                minPt = 10.0,
                #debug = cms.bool(True),
            )
        )
        
        TICLphoVar_task.add(getattr(process, prodLabel))
        
        l_mapProdVars_TICLpho.append(cms.InputTag(prodLabel, getattr(process, prodLabel).instanceName.value(), processName))
        l_var_TICLpho.append("%s_%s" %(prodLabel, getattr(process, prodLabel).instanceName.value()))
    
    
    ## Track-iso
    #prodLabel = "TICLeleTrackIsoProducer%s" %(radStr)
    #
    #setattr(
    #    process,
    #    prodLabel,
    #    HGCalElectronTrackIsoProducer.clone(
    #        electrons = label_gsfEleFromTICL,
    #        isoConeDR = cone,
    #        minTrackPt = 0.0,
    #        maxTrackEleDz = 99999.0,
    #        #debug = cms.bool(True),
    #    )
    #)
    #
    #TICLeleVar_task.add(getattr(process, prodLabel))
    #
    #l_mapProdVars_TICLele.append(
    #    cms.InputTag(prodLabel, getattr(process, prodLabel).instanceName.value(), processName)
    #)
    #
    #
    #prodLabel = "TICLeleTrackIsoProducer%sTrkDz0p15TrkPt1" %(radStr)
    #
    #setattr(
    #    process,
    #    prodLabel,
    #    HGCalElectronTrackIsoProducer.clone(
    #        electrons = label_gsfEleFromTICL,
    #        isoConeDR = cone,
    #        #debug = cms.bool(True),
    #    )
    #)
    #
    #TICLeleVar_task.add(getattr(process, prodLabel))
    #
    #l_mapProdVars_TICLele.append(
    #    cms.InputTag(prodLabel, getattr(process, prodLabel).instanceName.value(), processName)
    #)


for iRad, rad in enumerate(l_rad) :
    
    for iEnCut, enCut in enumerate(l_hitEnCut) :
        
        radStr = "R%s" %(str(rad).replace(".", "p"))
        enCutStr = "En%s" %(str(enCut).replace(".", "p"))
        
        
        # Ele Rvar
        prodLabel = "TICLeleRvarProducer%s%s" %(radStr, enCutStr)
        
        setattr(
            process,
            prodLabel,
            HGCalElectronRvarProducer.clone(
                electrons = label_gsfEleFromTICL,
                cylinderR = rad,
                minHitE = enCut,
                minPt = 10.0,
                #debug = cms.bool(True),
            )
        )
        
        TICLeleVar_task.add(getattr(process, prodLabel))
        
        l_mapProdVars_TICLele.append(cms.InputTag(prodLabel, getattr(process, prodLabel).instanceName.value(), processName))
        l_var_TICLele.append("%s_%s" %(prodLabel, getattr(process, prodLabel).instanceName.value()))
        
        
        # Pho Rvar
        prodLabel = "TICLphoRvarProducer%s%s" %(radStr, enCutStr)
        
        setattr(
            process,
            prodLabel,
            HGCalPhotonRvarProducer.clone(
                photons = label_phoFromTICL,
                cylinderR = rad,
                minHitE = enCut,
                minPt = 10.0,
                #debug = cms.bool(True),
            )
        )
        
        TICLphoVar_task.add(getattr(process, prodLabel))
        
        l_mapProdVars_TICLpho.append(cms.InputTag(prodLabel, getattr(process, prodLabel).instanceName.value(), processName))
        l_var_TICLpho.append("%s_%s" %(prodLabel, getattr(process, prodLabel).instanceName.value()))
        
        
        # Ele PCA
        prodLabel = "TICLelePCAProducer%s%s" %(radStr, enCutStr)
        
        setattr(
            process,
            prodLabel,
            HGCalElectronPCAProducer.clone(
                electrons = label_gsfEleFromTICL,
                cylinderR = rad,
                minHitE = enCut,
                minPt = 10.0,
                #debug = cms.bool(True),
            )
        )
        
        TICLeleVar_task.add(getattr(process, prodLabel))
        
        l_mapProdVars_TICLele.append(cms.InputTag(prodLabel, getattr(process, prodLabel).instanceName.value()+"Sigma2UU", processName))
        l_mapProdVars_TICLele.append(cms.InputTag(prodLabel, getattr(process, prodLabel).instanceName.value()+"Sigma2VV", processName))
        l_mapProdVars_TICLele.append(cms.InputTag(prodLabel, getattr(process, prodLabel).instanceName.value()+"Sigma2WW", processName))
        
        l_var_TICLele.append("%s_%s" %(prodLabel, getattr(process, prodLabel).instanceName.value()+"Sigma2UU"))
        l_var_TICLele.append("%s_%s" %(prodLabel, getattr(process, prodLabel).instanceName.value()+"Sigma2VV"))
        l_var_TICLele.append("%s_%s" %(prodLabel, getattr(process, prodLabel).instanceName.value()+"Sigma2WW"))
        
        
        # Pho PCA
        prodLabel = "TICLphoPCAProducer%s%s" %(radStr, enCutStr)
        
        setattr(
            process,
            prodLabel,
            HGCalPhotonPCAProducer.clone(
                photons = label_phoFromTICL,
                cylinderR = rad,
                minHitE = enCut,
                minPt = 10.0,
                #debug = cms.bool(True),
            )
        )
        
        TICLphoVar_task.add(getattr(process, prodLabel))
        
        l_mapProdVars_TICLpho.append(cms.InputTag(prodLabel, getattr(process, prodLabel).instanceName.value()+"Sigma2UU", processName))
        l_mapProdVars_TICLpho.append(cms.InputTag(prodLabel, getattr(process, prodLabel).instanceName.value()+"Sigma2VV", processName))
        l_mapProdVars_TICLpho.append(cms.InputTag(prodLabel, getattr(process, prodLabel).instanceName.value()+"Sigma2WW", processName))
        
        l_var_TICLpho.append("%s_%s" %(prodLabel, getattr(process, prodLabel).instanceName.value()+"Sigma2UU"))
        l_var_TICLpho.append("%s_%s" %(prodLabel, getattr(process, prodLabel).instanceName.value()+"Sigma2VV"))
        l_var_TICLpho.append("%s_%s" %(prodLabel, getattr(process, prodLabel).instanceName.value()+"Sigma2WW"))


# Ele cluster isolation
prodLabel = "TICLeleClusIsoProducer"

setattr(
    process,
    prodLabel,
    HGCalElectronClusIsoProducer.clone(
        electrons = label_gsfEleFromTICL,
        minPt = 10.0,
        #debug = cms.bool(True),
    )
)

TICLeleVar_task.add(getattr(process, prodLabel))
l_mapProdVars_TICLele.append(cms.InputTag(prodLabel, getattr(process, prodLabel).instanceName.value(), processName))
l_var_TICLele.append("%s_%s" %(prodLabel, getattr(process, prodLabel).instanceName.value()))


# Pho cluster isolation
prodLabel = "TICLphoClusIsoProducer"

setattr(
    process,
    prodLabel,
    HGCalPhotonClusIsoProducer.clone(
        photons = label_phoFromTICL,
        minPt = 10.0,
        #debug = cms.bool(True),
    )
)

TICLphoVar_task.add(getattr(process, prodLabel))
l_mapProdVars_TICLpho.append(cms.InputTag(prodLabel, getattr(process, prodLabel).instanceName.value(), processName))
l_var_TICLpho.append("%s_%s" %(prodLabel, getattr(process, prodLabel).instanceName.value()))


# TICL ele
process.HGCalElectronVarMap = mapProducer.clone(
    collections = cms.VInputTag(l_mapProdVars_TICLele),
    
    debug = cms.bool(False),
)

# TICL pho
process.HGCalPhotonVarMap = mapProducer.clone(
    collections = cms.VInputTag(l_mapProdVars_TICLpho),
    
    debug = cms.bool(False),
)


process.HGCalVar_seq = cms.Sequence(
    process.HGCalElectronVarMap *
    process.HGCalPhotonVarMap
)

process.HGCalVar_seq.associate(TICLeleVar_task)
process.HGCalVar_seq.associate(TICLphoVar_task)


process.treeMaker = cms.EDAnalyzer(
    "TreeMaker",
    
    ############################## My stuff ##############################
    debug = cms.bool(False),
    
    isGunSample = cms.bool(bool(options.isGunSample)),
    
    storeSimHit = cms.bool(bool(options.storeSimHit)),
    storeRecHit = cms.bool(bool(options.storeRecHit)),
    
    TICLeleGenMatchDR = cms.double(options.TICLeleGenMatchDR),
    TICLphoGenMatchDR = cms.double(options.TICLphoGenMatchDR),
    
    
    ############################## GEN ##############################
    
    label_generator = cms.InputTag("generator"),
    label_genParticle = cms.InputTag("genParticles"),
    
    
    ############################## RECO ##############################
    
    label_pileup = cms.InputTag("addPileupInfo"),
    label_rho = cms.InputTag("fixedGridRhoFastjetAll"),
    
    label_EcalEBRecHit = cms.InputTag("ecalRecHit", "EcalRecHitsEB"),
    
    label_HGCEESimHit = cms.InputTag("g4SimHits", "HGCHitsEE"),
    label_HGCHEFSimHit = cms.InputTag("g4SimHits", "HGCHitsHEfront"), # HE Silicon
    label_HGCEEBSimHit = cms.InputTag("g4SimHits", "HGCHitsHEback"), # HE Scintillator
    
    label_HGCEERecHit = cms.InputTag("HGCalRecHit", "HGCEERecHits"),
    label_HGCHEFRecHit = cms.InputTag("HGCalRecHit", "HGCHEFRecHits"),
    label_HGCHEBRecHit = cms.InputTag("HGCalRecHit", "HGCHEBRecHits"),
    
    label_HGCALlayerCluster = cms.InputTag("hgcalLayerClusters"),
    label_HGCALlayerClusterTime = cms.InputTag("hgcalLayerClusters", "timeLayerCluster"),
    
    label_TICLtrackster = label_TICLtrackster,
    label_TICLmultiCluster = label_TICLmultiCluster,
    
    label_PFRecHitHGC = cms.InputTag("particleFlowRecHitHGC"),
    
    label_caloParticle = cms.InputTag("mix", "MergedCaloTruth"),
    
    label_gsfEleFromMultiClus = cms.InputTag(""),
    label_phoFromMultiClus = cms.InputTag(""),
    
    label_gsfEleFromTICL = label_gsfEleFromTICL,
    label_gsfEleFromTICLvarList = cms.vstring(l_var_TICLele),
    label_gsfEleFromTICLvarMap = cms.InputTag("HGCalElectronVarMap", process.HGCalElectronVarMap.instanceName.value(), processName),
    
    label_generalTrack = cms.InputTag("generalTracks"),
    
    #label_phoFromMultiClus = cms.InputTag("photonsHGC", "", "RECO"),
    
    label_phoFromTICL = label_phoFromTICL,
    label_phoFromTICLvarList = cms.vstring(l_var_TICLpho),
    label_phoFromTICLvarMap = cms.InputTag("HGCalPhotonVarMap", process.HGCalPhotonVarMap.instanceName.value(), processName),
)


########## Filters ##########

from EDFilters.MyFilters.GenParticleFilter_cfi import *

# Gen-ele filter
process.GenParticleFilter_ele = GenParticleFilter.clone()
process.GenParticleFilter_ele.atLeastN = cms.int32(1)
process.GenParticleFilter_ele.pdgIds = cms.vint32(11)
process.GenParticleFilter_ele.minPt = cms.double(10)
process.GenParticleFilter_ele.minEta = cms.double(1.479)
process.GenParticleFilter_ele.maxEta = cms.double(3.1)
process.GenParticleFilter_ele.isGunSample = cms.bool(bool(options.isGunSample))
#process.GenParticleFilter_ele.debug = cms.bool(True)

process.filter_seq_genEle = cms.Sequence()

if (options.genEleFilter) :

    process.filter_seq_genEle = cms.Sequence(process.GenParticleFilter_ele)


# Gen-pho filter
process.GenParticleFilter_pho = GenParticleFilter.clone()
process.GenParticleFilter_pho.atLeastN = cms.int32(1)
process.GenParticleFilter_pho.pdgIds = cms.vint32(22)
process.GenParticleFilter_pho.minPt = cms.double(10)
process.GenParticleFilter_pho.minEta = cms.double(1.479)
process.GenParticleFilter_pho.maxEta = cms.double(3.1)
process.GenParticleFilter_pho.isGunSample = cms.bool(bool(options.isGunSample))
#process.GenParticleFilter_pho.debug = cms.bool(True)

process.filter_seq_genPho = cms.Sequence()

if (options.genPhoFilter) :

    process.filter_seq_genPho = cms.Sequence(process.GenParticleFilter_pho)


# Gen-parton filter
process.GenParticleFilter_part = GenParticleFilter.clone()
process.GenParticleFilter_part.atLeastN = cms.int32(1)
process.GenParticleFilter_part.pdgIds = cms.vint32(1, 2, 3, 4, 5, 21)
process.GenParticleFilter_part.minPt = cms.double(10)
process.GenParticleFilter_part.minEta = cms.double(1.479)
process.GenParticleFilter_part.maxEta = cms.double(3.1)
process.GenParticleFilter_part.isGunSample = cms.bool(bool(options.isGunSample))
#process.GenParticleFilter_part.debug = cms.bool(True)

process.filter_seq_genPart = cms.Sequence()

if (options.genPartonFilter) :

    process.filter_seq_genPart = cms.Sequence(process.GenParticleFilter_part)


print("Deleting existing output file.")
os.system("rm %s" %(outFile))


# Output file name modification
if (outFile.find("/eos/cms") ==  0) :
    
    outFile = outFile.replace("/eos/cms", "root://eoscms.cern.ch//eos/cms")


# Output
process.TFileService = cms.Service(
    "TFileService",
    fileName = cms.string(outFile)
)


process.schedule = cms.Schedule()


# Aging
from SLHCUpgradeSimulations.Configuration.aging import customise_aging_1000
customise_aging_1000(process)


process.reco_seq = cms.Sequence()

if (options.onRaw) :
    
    process.reco_seq = cms.Sequence(
        process.RawToDigi *
        #process.L1Reco *
        process.reconstruction_mod# *
        #process.recosim_step
    )
    
    process.reco_seq.associate(process.recosim)


###### PixelCPE issue
#process.TrackProducer.TTRHBuilder = "WithTrackAngle"
#process.PixelCPEGenericESProducer.UseErrorsFromTemplates = False
#process.PixelCPEGenericESProducer.LoadTemplatesFromDB = False
#process.PixelCPEGenericESProducer.TruncatePixelCharge = False
#process.PixelCPEGenericESProducer.IrradiationBiasCorrection = False
#process.PixelCPEGenericESProducer.DoCosmics = False
#process.PixelCPEGenericESProducer.Upgrade = cms.bool(True) 
######


process.p = cms.Path(
    
    process.filter_seq_genEle *
    process.filter_seq_genPho *
    process.filter_seq_genPart *
    
    process.reco_seq *
    
    process.HGCalVar_seq *
    process.treeMaker
)

process.schedule.insert(0, process.p)

print("\n")
print("*"*50)
print("process.schedule:", process.schedule)
print("*"*50)
#print "process.schedule.__dict__:", process.schedule.__dict__
#print "*"*50
print("\n")


# Tracer
if (options.trace) :
    
    process.Tracer = cms.Service("Tracer")


if (options.memoryCheck) :
    
    process.SimpleMemoryCheck = cms.Service(
        "SimpleMemoryCheck",
        moduleMemorySummary = cms.untracked.bool(True),
    )


#Timing
if (options.printTime) :

    process.Timing = cms.Service("Timing",
        summaryOnly = cms.untracked.bool(False),
        useJobReport = cms.untracked.bool(True)
    )


# Debug
if (options.debugFile) :
    
    process.out = cms.OutputModule("PoolOutputModule",
        fileName = cms.untracked.string("debug.root")
    )
    
    process.output_step = cms.EndPath(process.out)
    process.schedule.extend([process.output_step])


process.MessageLogger = cms.Service(
    "MessageLogger",
    destinations = cms.untracked.vstring(
        "cerr",
    ),
    cerr = cms.untracked.PSet(
        #threshold  = cms.untracked.string("ERROR"),
        DEBUG = cms.untracked.PSet(limit = cms.untracked.int32(0)),
        WARNING = cms.untracked.PSet(limit = cms.untracked.int32(0)),
        ERROR = cms.untracked.PSet(limit = cms.untracked.int32(0)),
    )
)


#from FWCore.ParameterSet.Utilities import convertToUnscheduled
#process = convertToUnscheduled(process)


# Add early deletion of temporary data products to reduce peak memory need
#from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
#process = customiseEarlyDelete(process)
# End adding early deletion
