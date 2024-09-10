import FWCore.ParameterSet.Config as cms
from DQMServices.Core.DQMEDAnalyzer import DQMEDAnalyzer

ScoutingMonitoringAnalysis = DQMEDAnalyzer('ScoutingMonitoring',

    OutputInternalPath = cms.string('ScoutingMonitoring'),
    ElectronCollection = cms.InputTag('slimmedElectrons'),
    ScoutingElectronCollection = cms.InputTag("hltScoutingEgammaPacker::HLT"),
    eleIdMapTight = cms.InputTag('egmGsfElectronIDs:cutBasedElectronID-RunIIIWinter22-V1-loose')
)

scoutingMonitoring = cms.Sequence(ScoutingMonitoringAnalysis)
