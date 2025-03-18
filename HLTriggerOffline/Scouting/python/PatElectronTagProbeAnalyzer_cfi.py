import FWCore.ParameterSet.Config as cms
from DQMServices.Core.DQMEDAnalyzer import DQMEDAnalyzer

DoubleEGL1 = [
  "L1_DoubleEG_LooseIso16_LooseIso12_er1p5",
  "L1_DoubleEG_LooseIso18_LooseIso12_er1p5",
  "L1_DoubleEG_LooseIso20_LooseIso12_er1p5",
  "L1_DoubleEG_LooseIso22_LooseIso12_er1p5",
  "L1_DoubleEG11_er1p2_dR_Max0p6"
]



PatElectronTagProbeAnalysis = DQMEDAnalyzer('PatElectronTagProbeAnalyzer',

    OutputInternalPath = cms.string('/HLT/ScoutingOffline/EGamma/TnP/Tag_PatElectron'),

    triggerSelection = cms.vstring(["DST_PFScouting_ZeroBias_v*", "DST_PFScouting_DoubleEG_v*", "DST_PFScouting_JetHT_v*"]), #Denominator
    triggerConfiguration = cms.PSet(
        hltResults          = cms.InputTag('TriggerResults', '', 'HLT'),
        l1tResults          = cms.InputTag('','',''),
        l1tIgnoreMaskAndPrescale = cms.bool(False),
        throw               = cms.bool(True),
        usePathStatus       = cms.bool(False)
    ),
    AlgInputTag        = cms.InputTag("gtStage2Digis"),
    l1tAlgBlkInputTag  = cms.InputTag("gtStage2Digis"),
    l1tExtBlkInputTag  = cms.InputTag("gtStage2Digis"),
    L1Seeds            = cms.vstring(DoubleEGL1),
    ReadPrescalesFromFile = cms.bool(False),

    TriggerResultTag   = cms.InputTag("TriggerResults", "", "HLT"),
    TriggerObjects     = cms.InputTag("slimmedPatTrigger"),
    ElectronCollection = cms.InputTag('slimmedElectrons'),
    ScoutingElectronCollection = cms.InputTag('hltScoutingEgammaPacker'),
    eleIdMapTight = cms.InputTag('egmGsfElectronIDsForScoutingDQM:cutBasedElectronID-RunIIIWinter22-V1-tight')

)

scoutingMonitoringPatElectronTagProbe = cms.Sequence(PatElectronTagProbeAnalysis)
