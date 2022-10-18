#!/bin/bash

cmsRun EDAnalyzers/TreeMaker/python/ConfFile_cfg.py \
sourceFile=\
sourceFiles/SinglePhoton_Pt-2To200-gun_PhaseIISpring22DRMiniAOD-PU200_123X_mcRun4_realistic_v11-v1_GEN-SIM-DIGI-RAW-MINIAOD/SinglePhoton_Pt-2To200-gun_PhaseIISpring22DRMiniAOD-PU200_123X_mcRun4_realistic_v11-v1_GEN-SIM-DIGI-RAW-MINIAOD.txt \
conditionOpt=PhaseIISpring22DRMiniAOD \
isGunSample=1 \
onRaw=1 \
debugFile=0 \
maxEvents=5 \
outFileBaseName="ntupleTree_SinglePhoton_PT2to200_PhaseIISpring22DRMiniAOD-PU200"
