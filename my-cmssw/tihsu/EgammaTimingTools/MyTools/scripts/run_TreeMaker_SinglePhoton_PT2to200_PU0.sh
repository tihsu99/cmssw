#!/bin/bash


cmsRun EDAnalyzers/TreeMaker/python/ConfFile_cfg.py \
    sourceFile=\
sourceFiles/SinglePhoton_PT2to200_Phase2HLTTDRSummer20ReRECOMiniAOD-NoPU_111X_mcRun4_realistic_T15_v1-v1_FEVT/SinglePhoton_PT2to200_Phase2HLTTDRSummer20ReRECOMiniAOD-NoPU_111X_mcRun4_realistic_T15_v1-v1_FEVT.txt \
    isGunSample=1 \
    onRaw=1 \
    debugFile=0 \
    maxEvents=500 \

