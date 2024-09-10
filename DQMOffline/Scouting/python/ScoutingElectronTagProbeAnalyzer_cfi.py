import FWCore.ParameterSet.Config as cms
from DQMServices.Core.DQMEDAnalyzer import DQMEDAnalyzer

ScoutingElectronTagProbeAnalysis = DQMEDAnalyzer('ScoutingElectronTagProbeAnalyzer',

    OutputInternalPath = cms.string('ScoutingMonitoring'),
    ElectronCollection = cms.InputTag('slimmedElectrons'),
    ScoutingElectronCollection = cms.InputTag('hltScoutingEgammaPacker')

)

scoutingMonitoringTagProbe = cms.Sequence(ScoutingElectronTagProbeAnalysis)
