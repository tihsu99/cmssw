#!/usr/bin/env python3

from __future__ import print_function

import os


l_sampleName = [
    "sampleName1",
    "sampleName2"
]


l_sourceFile = ["sourceFiles/{sampleName}.txt".format(sampleName = sampleName) for sampleName in l_sampleName]


cmd = ("python3 EgammaTimingTools/MyTools/scripts/run_condor.py "
    "--processNames {processNames} "
    "--inputFileLists {inputFileList} "
    "--cmsRunFile <cmsRun config file> "
    "--outputDir <output directory> "
    "--nUnitPerJob 1 "
    #"--nInputFileMax <N> " # Will use N files from each sample
    #"--cmsRunOptions \"isGunSample=1 onRaw=1\" " # Options for the cmsRun config
    #"--test " # Will create the job files, but will not submit
).format(
    processNames = " ".join(l_sampleName),
    inputFileList = " ".join(l_sourceFile),
)

print(cmd)


os.system(cmd)
