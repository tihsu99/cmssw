#!/usr/bin/env python3

from __future__ import print_function

import os


l_sampleName = [
    "SinglePhoton_PT2to200_Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_111X_mcRun4_realistic_T15_v1_ext2-v1_FEVT",
]


l_sourceFile = ["sourceFiles/{sampleName}/{sampleName}.txt".format(sampleName = sampleName) for sampleName in l_sampleName]


cmd = ("python3 EgammaTimingTools/MyTools/scripts/run_condor.py "
    "--processNames {processNames} "
    "--inputFileLists {inputFileList} "
    "--cmsRunFile EDAnalyzers/TreeMaker/python/ConfFile_cfg.py "
    "--outputDir /eos/cms/store/group/phys_egamma/ec/`whoami` "
    "--nUnitPerJob 1 "
    "--nInputFileMax 1 " # Will use N files from each sample
    "--cmsRunOptions \"isGunSample=1 onRaw=1\" " # Options for the cmsRun config
    "--test " # Will create the job files, but will not submit
).format(
    processNames = " ".join(l_sampleName),
    inputFileList = " ".join(l_sourceFile),
)

print(cmd)


os.system(cmd)
