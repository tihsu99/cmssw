import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
# Enable LogInfo
# process.MessageLogger.cerr.INFO = cms.untracked.PSet(
#     limit = cms.untracked.int32(-1)
# )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'root://cms-xrd-global.cern.ch///store/data/Run2024C/ScoutingPFMonitor/MINIAOD/PromptReco-v1/000/379/420/00000/5aaef7fd-ab63-497c-ac23-ca1f98d99ab8.root'
    )
)

process.load("DQMOffline.Scouting.ScoutingMonitoring_cfi")
process.DQMStore = cms.Service("DQMStore")

process.load("DQMServices.FileIO.DQMFileSaverOnline_cfi")
process.dqmSaver.tag = 'SCOUTMONIT'

process.p = cms.Path(process.scoutingMonitoring + process.dqmSaver)
