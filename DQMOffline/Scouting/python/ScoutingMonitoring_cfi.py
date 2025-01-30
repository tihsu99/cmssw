import FWCore.ParameterSet.Config as cms
from DQMServices.Core.DQMEDAnalyzer import DQMEDAnalyzer

ScoutingMonitoringAnalysis = DQMEDAnalyzer('ScoutingMonitoring',

    OutputInternalPath = cms.string('ScoutingMonitoring'),
    TriggerResultTag   = cms.InputTag("TriggerResults", "", "HLT"),
    ElectronCollection = cms.InputTag('slimmedElectrons'),
    ScoutingElectronCollection = cms.InputTag("hltScoutingEgammaPacker"),
    eleIdMapTight = cms.InputTag('egmGsfElectronIDsForScoutingDQM:cutBasedElectronID-RunIIIWinter22-V1-loose')
)


scoutingMonitoring = cms.Sequence(ScoutingMonitoringAnalysis)
