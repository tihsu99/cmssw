// -*- C++ -*-
//
// Package:    DQMOffline/Scouting
// Class:      ScoutingMonitoring
//
/**\class ScoutingMonitoring ScoutingMonitoring.cc
 DQMOffline/Scouting/plugins/ScoutingMonitoring.cc

 Description: ScoutingMonitoring is developed to enable us to monitor the
              comparison between pat::Object and Run3Scouting<Object>.

 Implementation:
     * Current runs on top of MINIAOD dataformat of the ScoutingMonitoring
 dataset.
     * Implemented only for electrons as of now.
*/
//
// Original Author:  Abanti Ranadhir Sahasransu
//         Created:  Sun, 18 Aug 2024 13:02:11 GMT
//
//
#ifndef DQMOffline_Scouting_ScoutingMonitoring_h
#define DQMOffline_Scouting_ScoutingMonitoring_h

#include <string>
#include <vector>

// user include files
#include "DQMServices/Core/interface/DQMGlobalEDAnalyzer.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/Scouting/interface/Run3ScoutingElectron.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEventWithRefs.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerRefsCollections.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"

//
// class declaration
//

struct kThreeMomentumHistos {
  dqm::reco::MonitorElement* h1Pt;
  dqm::reco::MonitorElement* h1Eta;
  dqm::reco::MonitorElement* h1Phi;
};

struct kInvmHistos {
  dqm::reco::MonitorElement* h1N;
  kThreeMomentumHistos electrons;
  kThreeMomentumHistos electron1;
  kThreeMomentumHistos electron2;
  dqm::reco::MonitorElement* h1InvMass12;
  dqm::reco::MonitorElement* h1InvMassID;
  dqm::reco::MonitorElement* h1InvMassIDEBEB;
  dqm::reco::MonitorElement* h1InvMassIDEBEE;
  dqm::reco::MonitorElement* h1InvMassIDEEEE;
  dqm::reco::MonitorElement* h1InvMassID_passDoubleEG_DST;
  dqm::reco::MonitorElement* h1InvMassIDEBEB_passDoubleEG_DST;
  dqm::reco::MonitorElement* h1InvMassIDEBEE_passDoubleEG_DST;
  dqm::reco::MonitorElement* h1InvMassIDEEEE_passDoubleEG_DST;
  dqm::reco::MonitorElement* h1InvMassID_passSinglePhoton_DST;
  dqm::reco::MonitorElement* h1InvMassIDEBEB_passSinglePhoton_DST;
  dqm::reco::MonitorElement* h1InvMassIDEBEE_passSinglePhoton_DST;
  dqm::reco::MonitorElement* h1InvMassIDEEEE_passSinglePhoton_DST;
};

struct kHistogramsScoutingMonitoring {
  kInvmHistos patElectron;
  kInvmHistos sctElectron;
};

class ScoutingMonitoring
    : public DQMGlobalEDAnalyzer<kHistogramsScoutingMonitoring> {
 public:
  explicit ScoutingMonitoring(const edm::ParameterSet&);
  ~ScoutingMonitoring() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
 private:
  void bookHistograms(DQMStore::IBooker&, edm::Run const&,
                      edm::EventSetup const&,
                      kHistogramsScoutingMonitoring&) const override;

  void dqmAnalyze(edm::Event const&, edm::EventSetup const&,
                  kHistogramsScoutingMonitoring const&) const override;
  bool scoutingElectronID(const Run3ScoutingElectron el) const;
  bool scoutingElectronGsfTrackID(const Run3ScoutingElectron el, size_t) const;
  bool scoutingElectronGsfTrackIdx(const Run3ScoutingElectron el, size_t &) const;
  bool hasPatternInHLTPath(const edm::TriggerNames& triggerNames, const std::string& pattern) const;

  // ------------ member data ------------
  std::string outputInternalPath_;
  edm::EDGetToken triggerResultsToken_;
  edm::EDGetTokenT<edm::View<pat::Electron> > electronCollection_;
  edm::EDGetTokenT<std::vector<Run3ScoutingElectron> >
      scoutingElectronCollection_;
  edm::EDGetTokenT<edm::ValueMap<bool> > eleIdMapTightToken_;
};

#endif
