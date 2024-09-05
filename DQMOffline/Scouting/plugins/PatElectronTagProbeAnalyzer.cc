#include "PatElectronTagProbeAnalyzer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

PatElectronTagProbeAnalyzer::PatElectronTagProbeAnalyzer(const edm::ParameterSet& iConfig)
  : outputInternalPath_(iConfig.getParameter<std::string>("OutputInternalPath")),
    electronCollection_(consumes<edm::View<pat::Electron>>(iConfig.getParameter<edm::InputTag>("ElectronCollection"))),
    scoutingElectronCollection_(consumes<std::vector<Run3ScoutingElectron>>(iConfig.getParameter<edm::InputTag>("ScoutingElectronCollection"))),
    eleIdMapTightToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleIdMapTight"))){
}

PatElectronTagProbeAnalyzer::~PatElectronTagProbeAnalyzer(){
}

void PatElectronTagProbeAnalyzer::dqmAnalyze(edm::Event const& iEvent,
                                                  edm::EventSetup const& iSetup,
                                                  kTagProbeHistos const& histos) const {

  edm::Handle<edm::View<pat::Electron>> patEls;
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

  edm::Handle<edm::ValueMap<bool> > tight_ele_id_decisions;
  iEvent.getByToken(eleIdMapTightToken_ ,tight_ele_id_decisions);


  edm::LogInfo("ScoutingMonitoring")
      << "Process pat::Electrons: " << patEls->size();
  edm::LogInfo("ScoutingMonitoring")
      << "Process Run3ScoutingElectrons: " << sctEls->size();

 // for (const auto& pat_el : *patEls){
  for (size_t i = 0; i < patEls->size(); ++i){
    const auto pat_el = patEls->ptrAt(i);
    if (!((*tight_ele_id_decisions)[pat_el])) continue;
   
    math::PtEtaPhiMLorentzVector tag_pat_el(pat_el->pt(), pat_el->eta(), pat_el->phi(), pat_el->mass());
    for (size_t j =0; j < patEls->size(); ++j){
//    for (const auto& pat_el_second : *patEls){
    const auto pat_el_second = patEls->ptrAt(j);
      if(i == j) continue;
      math::PtEtaPhiMLorentzVector probe_pat_el(pat_el_second->pt(), pat_el_second->eta(), pat_el_second->phi(), pat_el_second->mass());
      float invMass = (tag_pat_el + probe_pat_el).mass();
      edm::LogInfo("ScoutingMonitoring")
          << "Inv Mass: " << invMass;
      if((80 < invMass) && (invMass < 100)){  fillHistograms_resonance(histos.resonanceZ,  *pat_el_second, invMass);
                                              fillHistograms_resonance(histos.resonanceAll,*pat_el_second, invMass);}
      if((2.8 < invMass) && (invMass < 3.8)){ fillHistograms_resonance(histos.resonanceJ,  *pat_el_second, invMass); // J/Psi mass: 3.3 +/- 0.2 GeV
                                              fillHistograms_resonance(histos.resonanceAll,*pat_el_second, invMass);}
      if((9.0 < invMass) && (invMass < 12.6)){ fillHistograms_resonance(histos.resonanceY,  *pat_el_second, invMass); // Y mass: 9.8 +/- 0.4 GeV & 10.6 +/- 1 GeV
                                               fillHistograms_resonance(histos.resonanceAll,*pat_el_second, invMass);}
      
    }
  }
}

bool PatElectronTagProbeAnalyzer::scoutingElectronID(const Run3ScoutingElectron el) const{

  bool isEB = (fabs(el.eta()) < 1.5);
  if (isEB){
    if(el.sigmaIetaIeta() > 0.015) return false;
    if(el.hOverE() > 0.2) return false;
    if(fabs(el.dEtaIn()) > 0.008) return false;
    if(fabs(el.dPhiIn()) > 0.06)  return false;
    if(el.ecalIso()/el.rawEnergy() > 0.25) return false;
    return true;

  }
  else{
    if(el.sigmaIetaIeta() > 0.045) return false;
    if(el.hOverE() > 0.2) return false;
    if(fabs(el.dEtaIn()) > 0.012) return false;
    if(fabs(el.dPhiIn()) > 0.06) return false;
    if(el.ecalIso()/el.rawEnergy() > 0.1) return false;
    return true;
  }
}

void PatElectronTagProbeAnalyzer::fillHistograms_resonance(const kProbeKinematicHistos histos, 
                                                           const pat::Electron el,
                                                           const float inv_mass) const{
  histos.hEta->Fill(el.eta());
  histos.hPhi->Fill(el.phi());
  histos.hInvMass->Fill(inv_mass);

  if(el.isEB()){
    histos.hPt_Barrel->Fill(el.pt());
    histos.hHoverE_Barrel->Fill(el.hadronicOverEm());
    histos.hOoEMOoP_Barrel->Fill((1.0 / el.ecalEnergy() - el.eSuperClusterOverP() / el.ecalEnergy()));
    histos.hdPhiIn_Barrel->Fill(fabs(el.deltaPhiSuperClusterTrackAtVtx()));
    histos.hdEtaIn_Barrel->Fill(fabs(el.deltaEtaSuperClusterTrackAtVtx()));
    histos.hSigmaIetaIeta_Barrel->Fill(el.sigmaIetaIeta());
    histos.hMissingHits_Barrel->Fill(el.gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS));
    histos.hTrackfbrem_Barrel->Fill(el.fbrem());
    histos.hRelEcalIsolation_Barrel->Fill(el.ecalIso() / el.pt());
    histos.hRelHcalIsolation_Barrel->Fill(el.hcalIso() / el.pt());
    histos.hRelTrackIsolation_Barrel->Fill(el.trackIso() / el.pt());
  }
  else{
    histos.hPt_Endcap->Fill(el.pt());
    histos.hHoverE_Endcap->Fill(el.hadronicOverEm());
    histos.hOoEMOoP_Endcap->Fill((1.0 / el.ecalEnergy() - el.eSuperClusterOverP() / el.ecalEnergy()));
    histos.hdPhiIn_Endcap->Fill(fabs(el.deltaPhiSuperClusterTrackAtVtx()));
    histos.hdEtaIn_Endcap->Fill(fabs(el.deltaEtaSuperClusterTrackAtVtx()));
    histos.hSigmaIetaIeta_Endcap->Fill(el.sigmaIetaIeta());
    histos.hMissingHits_Endcap->Fill(el.gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS));
    histos.hTrackfbrem_Endcap->Fill(el.fbrem());
    histos.hRelEcalIsolation_Endcap->Fill(el.ecalIso() / el.pt());
    histos.hRelHcalIsolation_Endcap->Fill(el.hcalIso() / el.pt());
    histos.hRelTrackIsolation_Endcap->Fill(el.trackIso() / el.pt());
  }
}

void PatElectronTagProbeAnalyzer::bookHistograms(DQMStore::IBooker& ibook,
                                                       edm::Run const& run,
                                                       edm::EventSetup const& iSetup,
                                                       kTagProbeHistos& histos) const{
    ibook.setCurrentFolder(outputInternalPath_);
    bookHistograms_resonance(ibook, run, iSetup, histos.resonanceZ, "resonanceZ");
    bookHistograms_resonance(ibook, run, iSetup, histos.resonanceJ, "resonanceJ");
    bookHistograms_resonance(ibook, run, iSetup, histos.resonanceY, "resonanceY");
    bookHistograms_resonance(ibook, run, iSetup, histos.resonanceAll, "resonanceAll");
 }

void  PatElectronTagProbeAnalyzer::bookHistograms_resonance(DQMStore::IBooker& ibook,
                                                                 edm::Run const& run,
                                                                 edm::EventSetup const& iSetup,
                                                                 kProbeKinematicHistos& histos,
                                                                 const std::string& name) const{
     ibook.setCurrentFolder(outputInternalPath_);
     histos.hPt_Barrel = 
       ibook.book1D(name + "_Probe_patElectron_Pt_Barrel",
                    name + "_Probe_patElectron_Pt_Barrel", 500, 0, 500);
     histos.hPt_Endcap = 
       ibook.book1D(name + "_Probe_patElectron_Pt_Endcap",
                    name + "_Probe_patElectron_Pt_Endcap", 500, 0, 500);

     histos.hEta = 
       ibook.book1D(name + "_Probe_patElectron_Eta",
                    name + "_Probe_patElectron_Eta", 100, -5.0, 5.0);
     histos.hPhi = 
       ibook.book1D(name + "_Probe_patElectron_Phi",
                    name + "_Probe_patElectron_Phi", 100, -3.3, 3.3);

     histos.hHoverE_Barrel = 
       ibook.book1D(name + "_Probe_patElectron_HoverE_Barrel",
                    name + "_Probe_patElectron_HoverE_Barrel", 500, 0, 0.5);
     histos.hHoverE_Endcap = 
       ibook.book1D(name + "_Probe_patElectron_HoverE_Endcap",
                    name + "_Probe_patElectron_HoverE_Endcap", 500, 0, 0.5);

     histos.hOoEMOoP_Barrel = 
       ibook.book1D(name + "_Probe_patElectron_OoEMOoP_Barrel",
                    name + "_Probe_patElectron_OoEMOoP_Barrel", 500, -0.2, 0.2);
     histos.hOoEMOoP_Endcap = 
       ibook.book1D(name + "_Probe_patElectron_OoEMOoP_Endcap",
                    name + "_Probe_patElectron_OoEMOoP_Endcap", 500, -0.2, 0.2);

     histos.hdPhiIn_Barrel = 
       ibook.book1D(name + "_Probe_patElectron_dPhiIn_Barrel",
                    name + "_Probe_patElectron_dPhiIn_Barrel", 100, 0, 0.1);
     histos.hdPhiIn_Endcap = 
       ibook.book1D(name + "_Probe_patElectron_dPhiIn_Endcap",
                    name + "_Probe_patElectron_dPhiIn_Endcap", 100, 0, 0.1);

     histos.hdEtaIn_Barrel = 
       ibook.book1D(name + "_Probe_patElectron_dEtaIn_Barrel",
                    name + "_Probe_patElectron_dEtaIn_Barrel", 100, 0, 0.1);
     histos.hdEtaIn_Endcap = 
       ibook.book1D(name + "_Probe_patElectron_dEtaIn_Endcap",
                    name + "_Probe_patElectron_dEtaIn_Endcap", 100, 0, 0.1);

     histos.hSigmaIetaIeta_Barrel = 
       ibook.book1D(name + "_Probe_patElectron_SigmaIetaIeta_Barrel",
                    name + "_Probe_patElectron_SigmaIetaIeta_Barrel", 500, 0, 0.05);
     histos.hSigmaIetaIeta_Endcap = 
       ibook.book1D(name + "_Probe_patElectron_SigmaIetaIeta_Endcap",
                    name + "_Probe_patElectron_SigmaIetaIeta_Endcap", 500, 0, 0.05);


     histos.hMissingHits_Barrel = 
       ibook.book1D(name + "_Probe_patElectron_MissingHits_Barrel",
                    name + "_Probe_patElectron_MissingHits_Barrel", 21, -0.5, 20.5);
     histos.hMissingHits_Endcap = 
       ibook.book1D(name + "_Probe_patElectron_MissingHits_Endcap",
                    name + "_Probe_patElectron_MissingHits_Endcap", 21, -0.5, 20.5);

     histos.hTrackfbrem_Barrel = 
       ibook.book1D(name + "_Probe_patElectron_Trackfbrem_Barrel",
                    name + "_Probe_patElectron_Trackfbrem_Barrel", 100, 0, 1.0);
     histos.hTrackfbrem_Endcap = 
       ibook.book1D(name + "_Probe_patElectron_Trackfbrem_Endcap",
                    name + "_Probe_patElectron_Trackfbrem_Endcap", 100, 0, 1.0);

     histos.hRelEcalIsolation_Barrel = 
       ibook.book1D(name + "_Probe_patElectron_RelEcalIsolation_Barrel",
                    name + "_Probe_patElectron_RelEcalIsolation_Barrel", 100, 0, 1.0);
     histos.hRelEcalIsolation_Endcap = 
       ibook.book1D(name + "_Probe_patElectron_RelEcalIsolation_Endcap",
                    name + "_Probe_patElectron_RelEcalIsolation_Endcap", 100, 0, 1.0);

     histos.hRelHcalIsolation_Barrel = 
       ibook.book1D(name + "_Probe_patElectron_RelHcalIsolation_Barrel",
                    name + "_Probe_patElectron_RelHcalIsolation_Barrel", 100, 0, 1.0);
     histos.hRelHcalIsolation_Endcap = 
       ibook.book1D(name + "_Probe_patElectron_RelHcalIsolation_Endcap",
                    name + "_Probe_patElectron_RelHcalIsolation_Endcap", 100, 0, 1.0);


     histos.hRelTrackIsolation_Barrel = 
       ibook.book1D(name + "_Probe_patElectron_RelTrackIsolation_Barrel",
                    name + "_Probe_patElectron_RelTrackIsolation_Barrel", 100, 0, 1.0);
     histos.hRelTrackIsolation_Endcap = 
       ibook.book1D(name + "_Probe_patElectron_RelTrackIsolation_Endcap",
                    name + "_Probe_patElectron_RelTrackIsolation_Endcap", 100, 0, 1.0);


     histos.hInvMass = 
       ibook.book1D(name + "_patElectron_Invariant_Mass",
                    name + "_patElectron_Invariant_Mass", 800, 0, 200);
}


// ------------ method fills 'descriptions' with the allowed parameters for the
// module  ------------
void PatElectronTagProbeAnalyzer::fillDescriptions(
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
  desc.add<edm::InputTag>("eleIdMapTight",
                          edm::InputTag("egmGsfElectronIDs:cutBasedElectronID-RunIIIWinter22-V1-tight"));
  descriptions.add("PatElectronTagProbeAnalyzer", desc);
}

DEFINE_FWK_MODULE(PatElectronTagProbeAnalyzer);
