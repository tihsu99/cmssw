import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
# Enable LogInfo
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
     limit = cms.untracked.int32(-1)
 )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'root://cms-xrd-global.cern.ch///store/data/Run2024C/ScoutingPFMonitor/MINIAOD/PromptReco-v1/000/379/420/00000/5aaef7fd-ab63-497c-ac23-ca1f98d99ab8.root'
    )
)


from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
# turn on VID producer, indicate data format  to be
# DataFormat.AOD or DataFormat.MiniAOD, as appropriate 

dataFormat = DataFormat.MiniAOD
switchOnVIDElectronIdProducer(process, dataFormat)

# define which IDs we want to produce
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Winter22_122X_V1_cff']

#add them to the VID producer
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

process.load("DQMOffline.Scouting.ScoutingMonitoring_cfi")
process.load("DQMOffline.Scouting.ScoutingElectronTagProbeAnalyzer_cfi")
process.load("DQMOffline.Scouting.PatElectronTagProbeAnalyzer_cfi")
process.DQMStore = cms.Service("DQMStore")

process.load("DQMServices.FileIO.DQMFileSaverOnline_cfi")
process.dqmSaver.tag = 'SCOUTMONIT'

process.content = cms.EDAnalyzer("EventContentAnalyzer")
process.p = cms.Path(process.egmGsfElectronIDSequence + process.scoutingMonitoring + process.scoutingMonitoringTagProbe + process.scoutingMonitoringPatElectronTagProbe + process.dqmSaver)
