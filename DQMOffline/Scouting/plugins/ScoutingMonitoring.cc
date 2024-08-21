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
#include "ScoutingMonitoring.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
ScoutingMonitoring::ScoutingMonitoring(const edm::ParameterSet& iConfig)
    : outputInternalPath_(
          iConfig.getParameter<std::string>("OutputInternalPath")),
      electronCollection_(consumes<std::vector<pat::Electron> >(
          iConfig.getParameter<edm::InputTag>("ElectronCollection"))),
      scoutingElectronCollection_(consumes<std::vector<Run3ScoutingElectron> >(
          iConfig.getParameter<edm::InputTag>("ScoutingElectronCollection"))) {
  // now do what ever initialization is needed
}

ScoutingMonitoring::~ScoutingMonitoring() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
}

//
// member functions
//

// ------------ method called for each event  ------------

void ScoutingMonitoring::dqmAnalyze(
    edm::Event const& iEvent, edm::EventSetup const& iSetup,
    kHistogramsScoutingMonitoring const& histos) const {
  edm::Handle<std::vector<pat::Electron> > patEls;
  iEvent.getByToken(electronCollection_, patEls);
  if (patEls.failedToGet()) {
    edm::LogWarning("ScoutingMonitoring")
        << "pat::Electron collection not found.";
    return;
  }

  edm::Handle<std::vector<Run3ScoutingElectron> > sctEls;
  iEvent.getByToken(scoutingElectronCollection_, sctEls);
  if (sctEls.failedToGet()) {
    edm::LogWarning("ScoutingMonitoring")
        << "Run3ScoutingElectron collection not found.";
    return;
  }

  edm::LogInfo("ScoutingMonitoring")
      << "Process pat::Electrons: " << patEls->size();
  edm::LogInfo("ScoutingMonitoring")
      << "Process Run3ScoutingElectrons: " << sctEls->size();

  // Loop to verify the sorting of pat::Electron collection - REMOVE IN ONLINE
  // DQM
  for (size_t i = 1; i < patEls->size(); ++i) {
    if (patEls->at(i - 1).pt() < patEls->at(i).pt()) {
      edm::LogWarning("ScoutingMonitoring")
          << "pat::Electron collection not sorted by PT in descending order"
          << " will result in random histo filling. \n"
          << "pat::Electron[" << i << "].pt() = " << patEls->at(i).pt() << "\n"
          << "pat::Electron[" << i + 1 << "].pt() = " << patEls->at(i + 1).pt();
    }
  }

  // Fill pat::Electron histograms
  histos.patElectron.h1N->Fill(patEls->size());
  for (const auto& el : *patEls) {
    histos.patElectron.electrons.h1Pt->Fill(el.pt());
    histos.patElectron.electrons.h1Eta->Fill(el.eta());
    histos.patElectron.electrons.h1Phi->Fill(el.phi());
  }

  // Expect pat electrons to be sorted by pt
  if (patEls->size() >= 1) {
    histos.patElectron.electron1.h1Pt->Fill(patEls->at(0).pt());
    histos.patElectron.electron1.h1Eta->Fill(patEls->at(0).eta());
    histos.patElectron.electron1.h1Phi->Fill(patEls->at(0).phi());
  }

  if (patEls->size() >= 2) {
    histos.patElectron.electron2.h1Pt->Fill(patEls->at(1).pt());
    histos.patElectron.electron2.h1Eta->Fill(patEls->at(1).eta());
    histos.patElectron.electron2.h1Phi->Fill(patEls->at(1).phi());
    histos.patElectron.h1InvMass12->Fill(
        (patEls->at(0).p4() + patEls->at(1).p4()).mass());
  }

  // Fill the Run3ScoutingElectron histograms. No sorting assumed.
  histos.sctElectron.h1N->Fill(sctEls->size());
  // unsigned int leadSctElIndx = 0, subleadSctElIndx = -1;
  for (const auto& el : *sctEls) {
    histos.sctElectron.electrons.h1Pt->Fill(el.pt());
    histos.sctElectron.electrons.h1Eta->Fill(el.eta());
    histos.sctElectron.electrons.h1Phi->Fill(el.phi());

    // Find the leading and subleading Run3ScoutingElectron
    // if (el.pt() > sctEls->at(leadSctElIndx).pt()) {
    //   subleadSctElIndx = leadSctElIndx;
    //   leadSctElIndx = &el - &sctEls->at(0);
    // } else if (el.pt() > sctEls->at(subleadSctElIndx).pt()) {
    //   subleadSctElIndx = &el - &sctEls->at(0);
    // }
  }
}

void ScoutingMonitoring::bookHistograms(
    DQMStore::IBooker& ibook, edm::Run const& run,
    edm::EventSetup const& iSetup,
    kHistogramsScoutingMonitoring& histos) const {
  ibook.setCurrentFolder(outputInternalPath_);
  histos.patElectron.h1N =
      ibook.book1D("all_patElectron_electrons_N", "all_patElectron_electrons_N",
                   20, 0., 20.);
  histos.patElectron.electrons.h1Pt =
      ibook.book1D("all_patElectron_electrons_Pt",
                   "all_patElectron_electrons_Pt", 5000, 0., 500.);
  histos.patElectron.electrons.h1Eta =
      ibook.book1D("all_patElectron_electrons_Eta",
                   "all_patElectron_electrons_Eta", 1000, -5., 5.);
  histos.patElectron.electrons.h1Phi =
      ibook.book1D("all_patElectron_electrons_Phi",
                   "all_patElectron_electrons_Phi", 660, -3.3, 3.3);
  histos.patElectron.electron1.h1Pt = ibook.book1D(
      "patElectron_leading_Pt", "patElectron_leading_Pt", 5000, 0., 500.);
  histos.patElectron.electron1.h1Eta = ibook.book1D(
      "patElectron_leading_Eta", "patElectron_leading_Eta", 1000, -5., 5.);
  histos.patElectron.electron1.h1Phi = ibook.book1D(
      "patElectron_leading_Phi", "patElectron_leading_Phi", 660, -3.3, 3.3);
  histos.patElectron.electron2.h1Pt = ibook.book1D(
      "patElectron_subleading_Pt", "patElectron_subleading_Pt", 5000, 0., 500.);
  histos.patElectron.electron2.h1Eta =
      ibook.book1D("patElectron_subleading_Eta", "patElectron_subleading_Eta",
                   1000, -5., 5.);
  histos.patElectron.electron2.h1Phi =
      ibook.book1D("patElectron_subleading_Phi", "patElectron_subleading_Phi",
                   660, -3.3, 3.3);
  histos.patElectron.h1InvMass12 = ibook.book1D(
      "patElectron_E1E2_invMass", "patElectron_E1E2_invMass", 2000, 0., 200.);

  histos.sctElectron.h1N =
      ibook.book1D("all_sctElectron_electrons_N", "all_sctElectron_electrons_N",
                   20, 0., 20.);
  histos.sctElectron.electrons.h1Pt =
      ibook.book1D("all_sctElectron_electrons_Pt",
                   "all_sctElectron_electrons_Pt", 5000, 0., 500.);
  histos.sctElectron.electrons.h1Eta =
      ibook.book1D("all_sctElectron_electrons_Eta",
                   "all_sctElectron_electrons_Eta", 1000, -5., 5.);
  histos.sctElectron.electrons.h1Phi =
      ibook.book1D("all_sctElectron_electrons_Phi",
                   "all_sctElectron_electrons_Phi", 660, -3.3, 3.3);
  histos.sctElectron.electron1.h1Pt = ibook.book1D(
      "sctElectron_leading_Pt", "sctElectron_leading_Pt", 5000, 0., 500.);
  histos.sctElectron.electron1.h1Eta = ibook.book1D(
      "sctElectron_leading_Eta", "sctElectron_leading_Eta", 1000, -5., 5.);
  histos.sctElectron.electron1.h1Phi = ibook.book1D(
      "sctElectron_leading_Phi", "sctElectron_leading_Phi", 660, -3.3, 3.3);
  histos.sctElectron.electron2.h1Pt = ibook.book1D(
      "sctElectron_subleading_Pt", "sctElectron_subleading_Pt", 5000, 0., 500.);
  histos.sctElectron.electron2.h1Eta =
      ibook.book1D("sctElectron_subleading_Eta", "sctElectron_subleading_Eta",
                   1000, -5., 5.);
  histos.sctElectron.electron2.h1Phi =
      ibook.book1D("sctElectron_subleading_Phi", "sctElectron_subleading_Phi",
                   660, -3.3, 3.3);
  histos.sctElectron.h1InvMass12 = ibook.book1D(
      "sctElectron_E1E2_invMass", "sctElectron_E1E2_invMass", 2000, 0., 200.);
}

// ------------ method fills 'descriptions' with the allowed parameters for the
// module  ------------
void ScoutingMonitoring::fillDescriptions(
    edm::ConfigurationDescriptions& descriptions) {
  // The following says we do not know what parameters are allowed so do no
  // validation Please change this to state exactly what you do use, even if it
  // is no parameters
  edm::ParameterSetDescription desc;
  desc.add<std::string>("OutputInternalPath", "MY_FOLDER");
  desc.add<edm::InputTag>("ElectronCollection",
                          edm::InputTag("slimmedElectrons"));
  desc.add<edm::InputTag>("ScoutingElectronCollection",
                          edm::InputTag("Run3ScoutingElectrons"));
  descriptions.add("ScoutingMonitoring", desc);
}

// define this as a plug-in
DEFINE_FWK_MODULE(ScoutingMonitoring);
