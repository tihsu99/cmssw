#ifndef DQMOffline_Scouting_PatElectronTagProbeAnalyzer_h
#define DQMOffline_Scouting_PatElectronTagProbeAnalyzer_h

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

/////////////////////////
//  Class declaration  //
/////////////////////////

struct kProbeKinematicHistos{
  dqm::reco::MonitorElement* hPt_Barrel;
  dqm::reco::MonitorElement* hPt_Endcap;
  dqm::reco::MonitorElement* hEta;
  dqm::reco::MonitorElement* hEtavPhi;
  dqm::reco::MonitorElement* hPhi;
  dqm::reco::MonitorElement* hHoverE_Barrel;
  dqm::reco::MonitorElement* hHoverE_Endcap;
  dqm::reco::MonitorElement* hOoEMOoP_Barrel;
  dqm::reco::MonitorElement* hOoEMOoP_Endcap;
  dqm::reco::MonitorElement* hdPhiIn_Barrel;
  dqm::reco::MonitorElement* hdPhiIn_Endcap;
  dqm::reco::MonitorElement* hdEtaIn_Barrel;
  dqm::reco::MonitorElement* hdEtaIn_Endcap;
  dqm::reco::MonitorElement* hSigmaIetaIeta_Barrel;
  dqm::reco::MonitorElement* hSigmaIetaIeta_Endcap;
  dqm::reco::MonitorElement* hMissingHits_Barrel;
  dqm::reco::MonitorElement* hMissingHits_Endcap;
  dqm::reco::MonitorElement* hTrackfbrem_Barrel;
  dqm::reco::MonitorElement* hTrackfbrem_Endcap;
  dqm::reco::MonitorElement* hRelEcalIsolation_Barrel;
  dqm::reco::MonitorElement* hRelEcalIsolation_Endcap;
  dqm::reco::MonitorElement* hRelHcalIsolation_Barrel;
  dqm::reco::MonitorElement* hRelHcalIsolation_Endcap;
  dqm::reco::MonitorElement* hRelTrackIsolation_Barrel;
  dqm::reco::MonitorElement* hRelTrackIsolation_Endcap;
  dqm::reco::MonitorElement* hInvMass;
};

struct kTagProbeHistos {
  kProbeKinematicHistos resonanceZ;
  kProbeKinematicHistos resonanceJ;
  kProbeKinematicHistos resonanceY;
  kProbeKinematicHistos resonanceAll;
};

class PatElectronTagProbeAnalyzer: public DQMGlobalEDAnalyzer<kTagProbeHistos> {
      public: 
        explicit PatElectronTagProbeAnalyzer(const edm::ParameterSet& conf);
        ~PatElectronTagProbeAnalyzer() override;
        static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

      private:
        void dqmAnalyze(const edm::Event & e, const edm::EventSetup & c, kTagProbeHistos const&) const override;
        void bookHistograms(DQMStore::IBooker &, edm::Run const &, edm::EventSetup const &, kTagProbeHistos &) const override;
        void bookHistograms_resonance(DQMStore::IBooker &, edm::Run const &, edm::EventSetup const &, kProbeKinematicHistos &, const std::string &) const;
        void fillHistograms_resonance(const kProbeKinematicHistos histos, const pat::Electron el, const float inv_mass) const;
        bool scoutingElectronID(const Run3ScoutingElectron el) const;

        // --------------------- member data  ----------------------
        std::string outputInternalPath_;
        edm::EDGetTokenT<edm::View<pat::Electron>> electronCollection_;
        edm::EDGetTokenT<std::vector<Run3ScoutingElectron>> scoutingElectronCollection_;
        edm::EDGetTokenT<edm::ValueMap<bool> > eleIdMapTightToken_;

    };

#endif
