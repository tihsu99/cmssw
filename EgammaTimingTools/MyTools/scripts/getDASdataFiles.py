import argparse
import os


# Argument parser
parser = argparse.ArgumentParser(formatter_class = argparse.RawTextHelpFormatter)


parser.add_argument(
    "--getCount",
    help = "Get the number of events",
    default = False,
    action = "store_true",
)

parser.add_argument(
    "--getFiles",
    help = "Get the list of files",
    default = False,
    action = "store_true",
)

parser.add_argument(
    "--replace",
    help = "Will replace \"str1\" with \"str2\": \"str1\" \"str2\" (e.g. \"/store\" \"root://dcache-cms-xrootd.desy.de://pnfs/desy.de/cms/tier2/store\")",
    type = str,
    nargs = 2,
    required = False,
    default = ["/store", "root://cms-xrd-global.cern.ch//store"],
)

parser.add_argument(
    "--createRucioRule",
    help = "Create Rucio rules to copy datasets",
    default = False,
    action = "store_true",
)

parser.add_argument(
    "--tier2site",
    help = "T2 site name for Rucio rules (e.g. T2_DE_DESY)",
    type = str,
    required = False,
    default = None,
)

parser.add_argument(
    "--sampleNames",
    help = "Sample names (if provided, will not use the internal list)",
    type = str,
    required = False,
    nargs = "*",
    default = None,
)



# Parse arguments
args = parser.parse_args()
d_args = vars(args)


outDir = "sourceFiles"

prefix = args.replace[1]
toReplace = args.replace[0]


l_sampleName = [
    
    "/DYToLL_M-50_TuneCP5_14TeV-pythia8/Phase2HLTTDRWinter20DIGI-PU200_pilot2_110X_mcRun4_realistic_v3-v2/GEN-SIM-DIGI-RAW",
    
    "/QCD_Pt-15to3000_TuneCP5_Flat_14TeV-pythia8/Phase2HLTTDRWinter20DIGI-PU200_castor_110X_mcRun4_realistic_v3-v2/GEN-SIM-DIGI-RAW",
]


if (args.sampleNames is not None) :
    
    l_sampleName = args.sampleNames


# Reconfirm rucio rule creation request with user input
rucioConfirmed = False

if (args.createRucioRule and args.tier2site is not None) :
    
    print("Samples:")
    print("\n".join(l_sampleName))
    print("")
    
    print("REALLY create Rucio rules?")
    inputStr = str(raw_input("Enter CONFIRM to confirm: ")).strip()
    
    rucioConfirmed = (inputStr == "CONFIRM")
    
    if (not rucioConfirmed) :
        
        print("Rucio rule creation not confirmed. Exiting...")
        exit()


for iSample, sampleName in enumerate(l_sampleName) :
    
    print "\n"
    print "*"*50
    print "Sample %d/%d: %s" %(iSample+1, len(l_sampleName), sampleName)
    print "*"*50
    
    
    if (args.getCount) :
        
        command = "dasgoclient -query=\"file dataset=%s | sum(file.nevents)\"" %(sampleName)
        os.system(command)
    
    
    if (args.getFiles) :
        
        sampleName_mod = sampleName[1:].replace("/", "_")
        
        outDir_mod = "%s/%s" %(outDir, sampleName_mod)
        
        command = "mkdir -p %s" %(outDir_mod)
        print "Command:", command
        print ""
        os.system(command)
        
        outFile = "%s/%s.txt" %(outDir_mod, sampleName_mod)
        
        command = "dasgoclient -query=\"file dataset=%s\" > %s" %(sampleName, outFile)
        print "Command:", command
        print ""
        os.system(command)
        
        fileContent = ""
        
        print "Replacing \"%s\" with \"%s\" in file." %(toReplace, prefix)
        print ""
        
        print "Number of lines:"
        os.system("wc -l %s" %(outFile))
        print ""
        
        with open(outFile, "r") as f :
            
            fileContent = f.read()
        
        fileContent = fileContent.replace(toReplace, prefix)
        
        with open(outFile, "w") as f :
            
            f.write(fileContent)
    
    # https://twiki.cern.ch/twiki/bin/view/CMS/Rucio
    # https://twiki.cern.ch/twiki/bin/view/CMSPublic/RucioUserDocsQuotas
    # https://twiki.cern.ch/twiki/bin/view/CMSPublic/RucioUserDocsRules
    if (rucioConfirmed) :
        
        #command = "rucio add-rule cms:%s 1 %s" %(sampleName, args.tier2site)
        command = "rucio add-rule --ask-approval cms:%s 1 %s" %(sampleName, args.tier2site)
        print "Command:", command
        print ""
        os.system(command)
