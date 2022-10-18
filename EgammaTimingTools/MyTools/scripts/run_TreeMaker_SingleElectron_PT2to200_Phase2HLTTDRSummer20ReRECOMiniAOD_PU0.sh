#!/bin/bash


cmsRun EDAnalyzers/TreeMaker/python/ConfFile_cfg.py \
sourceFile=\
sourceFiles/SingleElectron_PT2to200_Phase2HLTTDRSummer20ReRECOMiniAOD-NoPU_111X_mcRun4_realistic_T15_v1-v1_FEVT/SingleElectron_PT2to200_Phase2HLTTDRSummer20ReRECOMiniAOD-NoPU_111X_mcRun4_realistic_T15_v1-v1_FEVT.txt \
conditionOpt=Phase2HLTTDRSummer20ReRECOMiniAOD \
isGunSample=1 \
onRaw=1 \
debugFile=0 \
maxEvents=5000 \
outFileBaseName="ntupleTree_SingleElectron_PT2to200_Phase2HLTTDRSummer20ReRECOMiniAOD-NoPU"

