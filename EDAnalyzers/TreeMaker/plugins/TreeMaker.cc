// -*- C++ -*-
//
// Package:    EDAnalyzers/TreeMaker
// Class:      TreeMaker
//
/**\class TreeMaker TreeMaker.cc EDAnalyzers/TreeMaker/plugins/TreeMaker.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  
//         Created:  Sat, 11 May 2019 13:14:55 GMT
//
//


// system include files
# include <memory>

// user include files



# include "CommonTools/UtilAlgos/interface/TFileService.h"
# include "DataFormats/CaloRecHit/interface/CaloCluster.h"
# include "DataFormats/CaloTowers/interface/CaloTowerDefs.h"
# include "DataFormats/Common/interface/MapOfVectors.h"
# include "DataFormats/Common/interface/ValueMap.h"
# include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
# include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
# include "DataFormats/EgammaCandidates/interface/Photon.h"
# include "DataFormats/ForwardDetId/interface/HGCEEDetId.h"
# include "DataFormats/ForwardDetId/interface/HGCalDetId.h"
# include "DataFormats/FWLite/interface/ESHandle.h"
# include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
# include "DataFormats/HGCalReco/interface/Trackster.h"
# include "DataFormats/HGCRecHit/interface/HGCRecHit.h"
# include "DataFormats/HepMCCandidate/interface/GenParticle.h"
# include "DataFormats/JetReco/interface/PFJet.h"
# include "DataFormats/Math/interface/LorentzVector.h"
# include "DataFormats/Math/interface/deltaPhi.h"
# include "DataFormats/Math/interface/deltaR.h"
# include "DataFormats/ParticleFlowReco/interface/HGCalMultiCluster.h"
# include "DataFormats/ParticleFlowReco/interface/PFRecHit.h"
# include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
# include "DataFormats/TrackReco/interface/Track.h"
# include "DataFormats/TrackReco/interface/TrackFwd.h"
# include "DataFormats/VertexReco/interface/Vertex.h"
# include "FWCore/Framework/interface/Event.h"
# include "FWCore/Framework/interface/ESHandle.h"
# include "FWCore/Framework/interface/Frameworkfwd.h"
# include "FWCore/Framework/interface/MakerMacros.h"
# include "FWCore/Framework/interface/one/EDAnalyzer.h"
# include "FWCore/ParameterSet/interface/ParameterSet.h"
# include "FWCore/ServiceRegistry/interface/Service.h"
# include "FWCore/Utilities/interface/InputTag.h"
# include "Geometry/CaloTopology/interface/HGCalTopology.h"
# include "Geometry/Records/interface/CaloGeometryRecord.h"
# include "Geometry/Records/interface/IdealGeometryRecord.h"
# include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"
# include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"
# include "SimDataFormats/CaloAnalysis/interface/SimCluster.h"
# include "SimDataFormats/CaloHit/interface/PCaloHit.h"
# include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
# include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

# include "EDAnalyzers/TreeMaker/interface/Common.h"
# include "EDAnalyzers/TreeMaker/interface/Constants.h"
# include "EDAnalyzers/TreeMaker/interface/TreeOutputInfo.h"

# include <CLHEP/Matrix/Matrix.h>
# include <CLHEP/Vector/ThreeVector.h>
# include <CLHEP/Vector/ThreeVector.h>

# include <Compression.h>
# include <TH1F.h>
# include <TH2F.h>
# include <TMatrixD.h>
# include <TTree.h> 
# include <TVector2.h> 
# include <TVectorD.h> 

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.



double HGCal_minEta = 1.479;
double HGCal_maxEta = 3.1;

double el_minPt = 10; //15;
double el_maxPt = 99999; //30;

double ph_minPt = 10; //15;
double ph_maxPt = 99999; //30;

double _largeVal = 999999999;


class TreeMaker : public edm::one::EDAnalyzer<edm::one::SharedResources>
{
    public:
    
    explicit TreeMaker(const edm::ParameterSet&);
    ~TreeMaker();
    
    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
    
    
    private:
    
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;
    
    hgcal::RecHitTools recHitTools;
    
    int minLayer;
    int maxLayer;
    
    
    TreeOutputInfo::TreeOutput *treeOutput;
    
    
    // My stuff //
    bool debug;
    bool isGunSample;
    
    bool storeSimHit;
    bool storeRecHit;
    
    double TICLeleGenMatchDR;
    double TICLphoGenMatchDR;
    
    std::vector <std::string> v_gsfEleFromTICLvar;
    std::vector <std::string> v_phoFromTICLvar;
    
    
    // GenEventInfoProduct //
    edm::EDGetTokenT <GenEventInfoProduct> tok_generator;
    
    
    // Gen particles //
    edm::EDGetTokenT <std::vector <reco::GenParticle> > tok_genParticle;
    
    
    // Pileup //
    edm::EDGetTokenT <std::vector <PileupSummaryInfo> > tok_pileup;
    
    
    // Rho //
    edm::EDGetTokenT <double> tok_rho;
    
    
    // SimHits //
    edm::EDGetTokenT <std::vector <PCaloHit> > tok_HGCEESimHit;
    
    
    // RecHits //
    edm::EDGetTokenT <edm::SortedCollection <EcalRecHit, edm::StrictWeakOrdering <EcalRecHit> > > tok_EcalEBRechit;
    
    edm::EDGetTokenT <std::vector <reco::PFRecHit> > tok_PFRecHitHGC;
    edm::EDGetTokenT <edm::SortedCollection <HGCRecHit, edm::StrictWeakOrdering <HGCRecHit> > > tok_HGCEERecHit;
    edm::EDGetTokenT <edm::SortedCollection <HGCRecHit, edm::StrictWeakOrdering <HGCRecHit> > > tok_HGCHEFRecHit;
    edm::EDGetTokenT <edm::SortedCollection <HGCRecHit, edm::StrictWeakOrdering <HGCRecHit> > > tok_HGCHEBRecHit;
    
    
    // HGCAL layer clusters //
    edm::EDGetTokenT <reco::CaloClusterCollection> tok_HGCALlayerCluster;
    edm::EDGetTokenT <edm::ValueMap <std::pair<float,float> > > tok_HGCALlayerClusterTime;
    
    
    // TICL //
    edm::EDGetTokenT <std::vector <ticl::Trackster> > tok_TICLtrackster;
    edm::EDGetTokenT <std::vector <reco::PFCluster> > tok_TICLmultiCluster;
    //edm::EDGetTokenT <std::vector <reco::HGCalMultiCluster> > tok_TICLmultiClusterMIP;
    
    
    // Calo particles //
    edm::EDGetTokenT <std::vector <CaloParticle> > tok_caloParticle;
    
    
    // Gsf electrons from multiclusters //
    edm::EDGetTokenT <std::vector <reco::GsfElectron> > tok_gsfEleFromMultiClus;
    
    
    // Gsf electrons from TICL //
    edm::EDGetTokenT <std::vector <reco::GsfElectron> > tok_gsfEleFromTICL;
    edm::EDGetTokenT <edm::MapOfVectors <std::string, double> > tok_gsfEleFromTICLvarMap;
    
    
    // Photons from multiclusters //
    edm::EDGetTokenT <std::vector <reco::Photon> > tok_phoFromMultiClus;
    
    
    // Photons from TICL //
    edm::EDGetTokenT <std::vector <reco::Photon> > tok_phoFromTICL;
    edm::EDGetTokenT <edm::MapOfVectors <std::string, double> > tok_phoFromTICLvarMap;
    
    
    // General tracks //
    edm::EDGetTokenT <std::vector <reco::Track> > tok_generalTrack;
    
    
    edm::ESGetToken<CaloGeometry, CaloGeometryRecord> tok_geom;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
TreeMaker::TreeMaker(const edm::ParameterSet& iConfig) :
    tok_geom{esConsumes<CaloGeometry, CaloGeometryRecord>()}
{
    usesResource("TFileService");
    edm::Service<TFileService> fs;
    
    // Compression
    //fs->file().SetCompressionAlgorithm(ROOT::kLZMA);
    //fs->file().SetCompressionLevel(8);
    
    
    minLayer = +9999;
    maxLayer = -9999;
    
    
    //now do what ever initialization is needed
    
    treeOutput = new TreeOutputInfo::TreeOutput("tree", fs);
    
    // My stuff //
    debug = iConfig.getParameter <bool>("debug");
    isGunSample = iConfig.getParameter <bool>("isGunSample");
    
    storeSimHit = iConfig.getParameter <bool>("storeSimHit");
    storeRecHit = iConfig.getParameter <bool>("storeRecHit");
    
    TICLeleGenMatchDR = iConfig.getParameter <double>("TICLeleGenMatchDR");
    TICLphoGenMatchDR = iConfig.getParameter <double>("TICLphoGenMatchDR");
    
    v_gsfEleFromTICLvar = iConfig.getParameter <std::vector <std::string> >("label_gsfEleFromTICLvarList");
    v_phoFromTICLvar = iConfig.getParameter <std::vector <std::string> >("label_phoFromTICLvarList");
    
    
    treeOutput->addToCustomVarMap("gsfEleFromTICL", v_gsfEleFromTICLvar);
    treeOutput->addToCustomVarMap("phoFromTICL", v_phoFromTICLvar);
    
    
    // GenEventInfoProduct //
    tok_generator = consumes <GenEventInfoProduct>(iConfig.getParameter <edm::InputTag>("label_generator"));

    
    // Gen particles //
    tok_genParticle = consumes <std::vector <reco::GenParticle> >(iConfig.getParameter <edm::InputTag>("label_genParticle"));
    
    
    // Pileup //
    tok_pileup = consumes <std::vector <PileupSummaryInfo> >(iConfig.getParameter <edm::InputTag>("label_pileup"));
    
    
    // Rho //
    tok_rho = consumes <double>(iConfig.getParameter <edm::InputTag>("label_rho"));
    
    
    // SimHits //
    tok_HGCEESimHit = consumes <std::vector <PCaloHit> >(iConfig.getParameter <edm::InputTag>("label_HGCEESimHit"));
    
    
    // RecHits //
    tok_EcalEBRechit = consumes <edm::SortedCollection <EcalRecHit, edm::StrictWeakOrdering <EcalRecHit> > >(iConfig.getParameter <edm::InputTag>("label_EcalEBRecHit"));
    
    tok_PFRecHitHGC = consumes <std::vector <reco::PFRecHit> >(iConfig.getParameter <edm::InputTag>("label_PFRecHitHGC"));
    
    tok_HGCEERecHit = consumes <edm::SortedCollection <HGCRecHit, edm::StrictWeakOrdering <HGCRecHit> > >(iConfig.getParameter <edm::InputTag>("label_HGCEERecHit"));
    tok_HGCHEFRecHit = consumes <edm::SortedCollection <HGCRecHit, edm::StrictWeakOrdering <HGCRecHit> > >(iConfig.getParameter <edm::InputTag>("label_HGCHEFRecHit"));
    tok_HGCHEBRecHit = consumes <edm::SortedCollection <HGCRecHit, edm::StrictWeakOrdering <HGCRecHit> > >(iConfig.getParameter <edm::InputTag>("label_HGCHEBRecHit"));
    
    
    // HGCAL layer clusters //
    tok_HGCALlayerCluster = consumes <reco::CaloClusterCollection>(iConfig.getParameter <edm::InputTag>("label_HGCALlayerCluster"));
    tok_HGCALlayerClusterTime = consumes <edm::ValueMap<std::pair<float,float> > >(iConfig.getParameter <edm::InputTag>("label_HGCALlayerClusterTime"));
    
    
    // TICL //
    tok_TICLtrackster = consumes <std::vector <ticl::Trackster> >(iConfig.getParameter <edm::InputTag>("label_TICLtrackster"));
    tok_TICLmultiCluster = consumes <std::vector <reco::PFCluster> >(iConfig.getParameter <edm::InputTag>("label_TICLmultiCluster"));
    //tok_TICLmultiClusterMIP = consumes <std::vector <reco::HGCalMultiCluster> >(iConfig.getParameter <edm::InputTag>("label_TICLmultiClusterMIP"));
    
    
    // Calo particles //
    tok_caloParticle = consumes <std::vector <CaloParticle> >(iConfig.getParameter <edm::InputTag>("label_caloParticle"));
    
    
    // Gsf electrons from multiclusters //
    tok_gsfEleFromMultiClus = consumes <std::vector <reco::GsfElectron> >(iConfig.getParameter <edm::InputTag>("label_gsfEleFromMultiClus"));
    
    
    // Gsf electrons from TICL //
    tok_gsfEleFromTICL = consumes <std::vector <reco::GsfElectron> >(iConfig.getParameter <edm::InputTag>("label_gsfEleFromTICL"));
    tok_gsfEleFromTICLvarMap = consumes <edm::MapOfVectors <std::string, double> >(iConfig.getParameter <edm::InputTag>("label_gsfEleFromTICLvarMap"));
    
    
    // Photons from multiclusters //
    tok_phoFromMultiClus = consumes <std::vector <reco::Photon> >(iConfig.getParameter <edm::InputTag>("label_phoFromMultiClus"));
    
    
    // Photons from TICL //
    tok_phoFromTICL = consumes <std::vector <reco::Photon> >(iConfig.getParameter <edm::InputTag>("label_phoFromTICL"));
    tok_phoFromTICLvarMap = consumes <edm::MapOfVectors <std::string, double> >(iConfig.getParameter <edm::InputTag>("label_phoFromTICLvarMap"));
    
    
    // General tracks //
    tok_generalTrack = consumes <std::vector <reco::Track> >(iConfig.getParameter <edm::InputTag>("label_generalTrack"));
}


TreeMaker::~TreeMaker()
{
    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)
    
    delete treeOutput;
}


//
// member functions
//


// ------------ method called for each event  ------------
void TreeMaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;
    
    long long eventNumber = iEvent.id().event();
    //printf("Event %llu \n", eventNumber);
    
    
    treeOutput->clear();
    
    //recHitTools.getEventSetup(iSetup);
    
    //edm::ESHandle<CaloGeometry> geom;
    //iSetup.get<CaloGeometryRecord>().get(geom);
    //recHitTools.setGeometry(*(geom.product()));

    const CaloGeometry *geom = &iSetup.getData(tok_geom);

    recHitTools.setGeometry(*geom);
    
    
    //////////////////// Run info ////////////////////
    treeOutput->runNumber = iEvent.id().run();
    treeOutput->eventNumber = iEvent.id().event();
    treeOutput->luminosityNumber = iEvent.id().luminosityBlock();
    treeOutput->bunchCrossingNumber = iEvent.bunchCrossing();
    
    
    // HGCal Topology
    /*edm::ESHandle <HGCalTopology> handle_topo_HGCalEE;
    iSetup.get<IdealGeometryRecord>().get("HGCalEESensitive", handle_topo_HGCalEE);
    
    if(!handle_topo_HGCalEE.isValid())
    {
        printf("Error: Invalid HGCalEE topology. \n");
        
        exit(EXIT_FAILURE);
    }
    
    const auto& topo_HGCalEE = *handle_topo_HGCalEE;
    */
    
    //////////////////// GenEventInfoProduct ////////////////////
    edm::Handle <GenEventInfoProduct> generatorHandle;
    iEvent.getByToken(tok_generator, generatorHandle);
    GenEventInfoProduct generator = *generatorHandle;
    
    printf("[%llu] Gen. evt. wt. %0.4g \n", eventNumber, generator.weight());
    treeOutput->genEventWeight = generator.weight();
    
    
    //////////////////// Gen particle ////////////////////
    edm::Handle <std::vector <reco::GenParticle> > v_genParticle;
    iEvent.getByToken(tok_genParticle, v_genParticle);
    
    std::vector <CLHEP::HepLorentzVector> v_genEl_4mom;
    
    std::vector <CLHEP::HepLorentzVector> v_genPh_4mom;
    
    
    for(int iPart = 0; iPart < (int) v_genParticle->size(); iPart++)
    {
        reco::GenParticle part = v_genParticle->at(iPart);
        
        int pdgId = part.pdgId();
        int status = part.status();
        
        // Gen ele
        //if(abs(pdgId) == 11 && status == 1)
        //if(abs(pdgId) == 11 && (part.isHardProcess() || status == 1))
        if(
            abs(pdgId) == 11 && (
                (isGunSample && status == 1) ||
                (!isGunSample && part.isHardProcess())
            )
        )
        {
            //printf("[%llu] Gen electron found: E %0.2f, pT %0.2f, eta %+0.2f \n", eventNumber, part.energy(), part.pt(), part.eta());
            
            printf(
                "[%llu] "
                "Gen-ele found: E %0.2f, pT %0.2f, eta %+0.2f, pz %+0.2f, "
                "\n",
                eventNumber,
                part.energy(), part.pt(), part.eta(), part.pz()
            );
            
            if(fabs(part.eta()) > HGCal_minEta && fabs(part.eta()) < HGCal_maxEta && part.pt() > el_minPt && part.pt() < el_maxPt)
            {
                CLHEP::HepLorentzVector genEl_4mom;
                
                genEl_4mom.setT(part.energy());
                genEl_4mom.setX(part.px());
                genEl_4mom.setY(part.py());
                genEl_4mom.setZ(part.pz());
                
                v_genEl_4mom.push_back(genEl_4mom);
                
                treeOutput->v_genEl_E.push_back(genEl_4mom.e());
                treeOutput->v_genEl_px.push_back(genEl_4mom.px());
                treeOutput->v_genEl_py.push_back(genEl_4mom.py());
                treeOutput->v_genEl_pz.push_back(genEl_4mom.pz());
                treeOutput->v_genEl_pT.push_back(genEl_4mom.perp());
                treeOutput->v_genEl_eta.push_back(genEl_4mom.eta());
                treeOutput->v_genEl_phi.push_back(genEl_4mom.phi());
                
                treeOutput->v_genEl_multiClus_totE.push_back(0);
                
                //if(part.eta() > 0)
                //{
                //    treeOutput->v_genEl_HGCalEEP_EsortedIndex.push_back(treeOutput->genEl_n);
                //}
                //
                //else
                //{
                //    treeOutput->v_genEl_HGCalEEM_EsortedIndex.push_back(treeOutput->genEl_n);
                //}
                
                treeOutput->genEl_n++;
            }
        }
        
        else if(
            abs(pdgId) == 22 && (
                (isGunSample && status == 1) ||
                (!isGunSample && part.isHardProcess())
            )
        )
        {
            printf(
                "[%llu] "
                "Gen-pho found: E %0.2f, pT %0.2f, eta %+0.2f, pz %+0.2f, "
                "\n",
                eventNumber,
                part.energy(), part.pt(), part.eta(), part.pz()
            );
            
            if(fabs(part.eta()) > HGCal_minEta && fabs(part.eta()) < HGCal_maxEta && part.pt() > ph_minPt && part.pt() < ph_maxPt)
            {
                CLHEP::HepLorentzVector genPh_4mom;
                
                genPh_4mom.setT(part.energy());
                genPh_4mom.setX(part.px());
                genPh_4mom.setY(part.py());
                genPh_4mom.setZ(part.pz());
                
                v_genPh_4mom.push_back(genPh_4mom);
                
                treeOutput->v_genPh_E.push_back(genPh_4mom.e());
                treeOutput->v_genPh_px.push_back(genPh_4mom.px());
                treeOutput->v_genPh_py.push_back(genPh_4mom.py());
                treeOutput->v_genPh_pz.push_back(genPh_4mom.pz());
                treeOutput->v_genPh_pT.push_back(genPh_4mom.perp());
                treeOutput->v_genPh_eta.push_back(genPh_4mom.eta());
                treeOutput->v_genPh_phi.push_back(genPh_4mom.phi());
                
                treeOutput->genPh_n++;
            }
        }
        
        else if(abs(pdgId) == 22 && Common::isPromptPhoton(part))
        {
            printf(
                "[%llu] "
                "Prompt gen-pho found: E %0.2f, pT %0.2f, eta %+0.2f, pz %+0.2f, "
                "\n",
                eventNumber,
                part.energy(), part.pt(), part.eta(), part.pz()
            );
            
            if(fabs(part.eta()) > HGCal_minEta && fabs(part.eta()) < HGCal_maxEta && part.pt() > ph_minPt && part.pt() < ph_maxPt)
            {
                CLHEP::HepLorentzVector genPh_4mom;
                
                genPh_4mom.setT(part.energy());
                genPh_4mom.setX(part.px());
                genPh_4mom.setY(part.py());
                genPh_4mom.setZ(part.pz());
                
                v_genPh_4mom.push_back(genPh_4mom);
                
                treeOutput->v_genPh_E.push_back(genPh_4mom.e());
                treeOutput->v_genPh_px.push_back(genPh_4mom.px());
                treeOutput->v_genPh_py.push_back(genPh_4mom.py());
                treeOutput->v_genPh_pz.push_back(genPh_4mom.pz());
                treeOutput->v_genPh_pT.push_back(genPh_4mom.perp());
                treeOutput->v_genPh_eta.push_back(genPh_4mom.eta());
                treeOutput->v_genPh_phi.push_back(genPh_4mom.phi());
                
                treeOutput->genPh_n++;
            }
        }
    }
    
    
    //if(!treeOutput->genEl_n)
    //{
    //    return;
    //}
    
    
    // Sort the electrons
    //std::sort(
    //    treeOutput->v_genEl_HGCalEEP_EsortedIndex.begin(), treeOutput->v_genEl_HGCalEEP_EsortedIndex.end(),
    //    [&](int iEle1, int iEle2)
    //    {
    //        return (treeOutput->v_genEl_E[iEle1] > treeOutput->v_genEl_E[iEle2]);
    //    }
    //);
    //
    //std::sort(
    //    treeOutput->v_genEl_HGCalEEM_EsortedIndex.begin(), treeOutput->v_genEl_HGCalEEM_EsortedIndex.end(),
    //    [&](int iEle1, int iEle2)
    //    {
    //        return (treeOutput->v_genEl_E[iEle1] > treeOutput->v_genEl_E[iEle2]);
    //    }
    //);
    
    
    // Sort the photons
    //std::sort(
    //    treeOutput->v_genPh_HGCalEEP_EsortedIndex.begin(), treeOutput->v_genPh_HGCalEEP_EsortedIndex.end(),
    //    [&](int iEle1, int iEle2)
    //    {
    //        return (treeOutput->v_genPh_E[iEle1] > treeOutput->v_genPh_E[iEle2]);
    //    }
    //);
    //
    //std::sort(
    //    treeOutput->v_genPh_HGCalEEM_EsortedIndex.begin(), treeOutput->v_genPh_HGCalEEM_EsortedIndex.end(),
    //    [&](int iEle1, int iEle2)
    //    {
    //        return (treeOutput->v_genPh_E[iEle1] > treeOutput->v_genPh_E[iEle2]);
    //    }
    //);
    
    
    // Pileup
    edm::Handle <std::vector <PileupSummaryInfo> > pileUps_reco;
    iEvent.getByToken(tok_pileup, pileUps_reco);
    treeOutput->pileup_n = Common::getPileup(pileUps_reco);
    
    
    // Rho
    edm::Handle <double> handle_rho;
    iEvent.getByToken(tok_rho, handle_rho);
    double rho = *handle_rho;
    
    treeOutput->rho = rho;
    
    
    // SimHit dictionary
    edm::Handle <std::vector <PCaloHit> > v_HGCEESimHit;
    iEvent.getByToken(tok_HGCEESimHit, v_HGCEESimHit);
    
    std::map <DetId, const PCaloHit*> m_simHit;
    
    int nSimHit = v_HGCEESimHit->size();
    
    for(int iSimHit = 0; iSimHit < nSimHit; iSimHit++)
    {
        const PCaloHit *simHit = &(v_HGCEESimHit->at(iSimHit));
        
        DetId detId(simHit->id());
        
        m_simHit[detId] = simHit;
    }
    
    
    // RecHit dictionary
    edm::Handle <edm::SortedCollection <EcalRecHit, edm::StrictWeakOrdering <EcalRecHit> > > v_EcalEBRecHit;
    iEvent.getByToken(tok_EcalEBRechit, v_EcalEBRecHit);
    
    edm::Handle <edm::SortedCollection <HGCRecHit, edm::StrictWeakOrdering <HGCRecHit> > > v_HGCEERecHit;
    iEvent.getByToken(tok_HGCEERecHit, v_HGCEERecHit);
    
    edm::Handle <edm::SortedCollection <HGCRecHit, edm::StrictWeakOrdering <HGCRecHit> > > v_HGCHEFRecHit;
    iEvent.getByToken(tok_HGCHEFRecHit, v_HGCHEFRecHit);
    
    edm::Handle <edm::SortedCollection <HGCRecHit, edm::StrictWeakOrdering <HGCRecHit> > > v_HGCHEBRecHit;
    iEvent.getByToken(tok_HGCHEBRecHit, v_HGCHEBRecHit);
    
    edm::Handle <std::vector <reco::PFRecHit> > v_PFRecHitHGC;
    iEvent.getByToken(tok_PFRecHitHGC, v_PFRecHitHGC);
    
    
    std::map <DetId, const HGCRecHit*> m_HGCRecHit;
    std::map <DetId, const HGCRecHit*> m_HGCEERecHit;
    std::unordered_map <DetId, const reco::PFRecHit*> m_PFRecHitHGC;
    
    int nHGCEERecHit = v_HGCEERecHit->size();
    
    for(int iRecHit = 0; iRecHit < nHGCEERecHit; iRecHit++)
    {
        const HGCRecHit *recHit = &(*v_HGCEERecHit)[iRecHit];
        
        m_HGCRecHit[recHit->id()] = recHit;
        m_HGCEERecHit[recHit->id()] = recHit;
    }
    
    
    //
    int nHGCHEFRecHit = v_HGCHEFRecHit->size();
    
    for(int iRecHit = 0; iRecHit < nHGCHEFRecHit; iRecHit++)
    {
        const HGCRecHit *recHit = &(*v_HGCHEFRecHit)[iRecHit];
        
        m_HGCRecHit[recHit->id()] = recHit;
    }
    
    
    //
    int nHGCHEBRecHit = v_HGCHEBRecHit->size();
    
    for(int iRecHit = 0; iRecHit < nHGCHEBRecHit; iRecHit++)
    {
        const HGCRecHit *recHit = &(*v_HGCHEBRecHit)[iRecHit];
        
        m_HGCRecHit[recHit->id()] = recHit;
    }
    
    
    //
    int nPFRecHitHGC = v_PFRecHitHGC->size();
    
    for(int iRecHit = 0; iRecHit < nPFRecHitHGC; iRecHit++)
    {
        const reco::PFRecHit *recHit = &(*v_PFRecHitHGC)[iRecHit];
        
        m_PFRecHitHGC[recHit->detId()] = recHit;
    }
    
    
    // General tracks
    edm::Handle <std::vector <reco::Track> > v_generalTrack;
    iEvent.getByToken(tok_generalTrack, v_generalTrack);
    
    
    /////////////////// HGCAL layer clusters ////////////////////
    edm::Handle <reco::CaloClusterCollection> v_HGCALlayerCluster;
    iEvent.getByToken(tok_HGCALlayerCluster, v_HGCALlayerCluster);
    
    edm::Handle <edm::ValueMap<std::pair<float,float> > > vm_HGCALlayerClusterTime;
    iEvent.getByToken(tok_HGCALlayerClusterTime, vm_HGCALlayerClusterTime);
    
    //int lcIdx = -1;
    //
    //for(const auto lc : (*v_HGCALlayerCluster))
    //{
    //    lcIdx++;
    //    
    //    edm::Ref <reco::CaloClusterCollection> lcRef(v_HGCALlayerCluster, lcIdx);
    //    std::pair<float,float> timeAndErr = (*vm_HGCALlayerClusterTime)[lcRef];
    //    
    //    printf("HGCAL layer cluster %d: E %f, size %d, time %f %f \n", lcIdx, lc.energy(), (int) lc.size(), timeAndErr.first, timeAndErr.second);
    //}
    
    
    // SimHits
    if(storeSimHit)
    {
        for(int iSimHit = 0; iSimHit < nSimHit; iSimHit++)
        {
            auto simHit = v_HGCEESimHit->at(iSimHit);
            
            DetId detId(simHit.id());
            
            int layer = recHitTools.getLayer(simHit.id()) - 1; // Start from 0
            int zside = recHitTools.zside(simHit.id());
            
            
            auto position = recHitTools.getPosition(simHit.id());
            
            if(zside > 0)
            {
                treeOutput->v_simHit_HGCalEEPlayer_totE.at(layer) += simHit.energy();
            }
            
            else
            {
                treeOutput->v_simHit_HGCalEEMlayer_totE.at(layer) += simHit.energy();
            }
            
            //bool isCaloParticleMatched = (m_caloPart_simClustHit.find(detId) != m_caloPart_simClustHit.end());
            
            
            treeOutput->v_simHit_E.push_back(simHit.energy());
            
            treeOutput->v_simHit_x.push_back(position.x());
            treeOutput->v_simHit_y.push_back(position.y());
            treeOutput->v_simHit_z.push_back(position.z());
            
            treeOutput->v_simHit_eta.push_back(position.eta());
            treeOutput->v_simHit_phi.push_back(position.phi());
            
            treeOutput->v_simHit_ET.push_back(recHitTools.getPt(position, simHit.energy()));
            
            treeOutput->v_simHit_layer.push_back(layer+1);
            treeOutput->v_simHit_zside.push_back(zside);
            //treeOutput->v_simHit_isCaloParticleMatched.push_back(isCaloParticleMatched);
            
            
            treeOutput->simHit_n++;
        }
    }
    
    
    /////////////////////////////////////////////////////////////////
    //////////////////// TICL electrons /////////////////////////////
    /////////////////////////////////////////////////////////////////
    edm::Handle <std::vector <reco::GsfElectron> > v_gsfEleFromTICL;
    iEvent.getByToken(tok_gsfEleFromTICL, v_gsfEleFromTICL);
    
    edm::Handle <edm::MapOfVectors <std::string, double> > m_gsfEleFromTICLvarMap;
    iEvent.getByToken(tok_gsfEleFromTICLvarMap, m_gsfEleFromTICLvarMap);
    
    int nEleFromTICL = v_gsfEleFromTICL->size();
    
    std::map <reco::SuperClusterRef, int> m_gsfEle_superClus;
    
    std::vector <CLHEP::HepLorentzVector> v_gsfEleFromTICL_4mom;
    
    
    for(int iEle = 0; iEle < nEleFromTICL; iEle++)
    {
        reco::GsfElectron gsfEle = v_gsfEleFromTICL->at(iEle);
        
        CLHEP::HepLorentzVector gsfEleFromTICL_4mom;
        gsfEleFromTICL_4mom.setT(gsfEle.energy());
        gsfEleFromTICL_4mom.setX(gsfEle.px());
        gsfEleFromTICL_4mom.setY(gsfEle.py());
        gsfEleFromTICL_4mom.setZ(gsfEle.pz());
        
        v_gsfEleFromTICL_4mom.push_back(gsfEleFromTICL_4mom);
    }
    
    
    // TICL-ele gen-matching
    //TMatrixD mat_gsfEleFromTICL_genEl_deltaR;
    //
    //std::vector <int> v_gsfEleFromTICL_matchedGenEl_idx;
    //
    //std::vector <double> v_gsfEleFromTICL_genEl_minDeltaR = Common::getMinDeltaR(
    //    v_gsfEleFromTICL_4mom,
    //    v_genEl_4mom,
    //    mat_gsfEleFromTICL_genEl_deltaR,
    //    v_gsfEleFromTICL_matchedGenEl_idx
    //);
    
    std::vector <int> v_gsfEleFromTICL_matchedGenEl_idx;
    std::vector <double> v_gsfEleFromTICL_matchedGenEl_deltaR;
    
    Common::getHardestInCone(
        v_gsfEleFromTICL_4mom,
        v_genEl_4mom,
        v_gsfEleFromTICL_matchedGenEl_idx,
        v_gsfEleFromTICL_matchedGenEl_deltaR,
        0.3
    );
    
    for(int iEle = 0; iEle < nEleFromTICL; iEle++)
    {
        reco::GsfElectron gsfEle = v_gsfEleFromTICL->at(iEle);
        CLHEP::HepLorentzVector gsfEleFromTICL_4mom = v_gsfEleFromTICL_4mom.at(iEle);
        
        
        if(gsfEle.pt() < el_minPt || fabs(gsfEle.eta()) < HGCal_minEta || fabs(gsfEle.eta()) > HGCal_maxEta)
        {
            continue;
        }
        
        
        //if(v_gsfEleFromTICL_genEl_minDeltaR.at(iEle) > TICLeleGenMatchDR)
        //{
        //    continue;
        //}
        
        
        //double matchedGenEl_deltaR = v_gsfEleFromTICL_genEl_minDeltaR.at(iEle);
        double matchedGenEl_deltaR = v_gsfEleFromTICL_matchedGenEl_deltaR.at(iEle);
        
        if(matchedGenEl_deltaR > TICLeleGenMatchDR)
        {
            continue;
        }
        
        int matchedGenEl_idx = v_gsfEleFromTICL_matchedGenEl_idx.at(iEle);
        
        treeOutput->v_gsfEleFromTICL_genEl_minDeltaR.push_back(matchedGenEl_deltaR);
        treeOutput->v_gsfEleFromTICL_nearestGenEl_idx.push_back(matchedGenEl_idx);
        
        double matchedGenEl_energy = -99;
        double matchedGenEl_pT = -99;
        double matchedGenEl_eta = -99;
        double matchedGenEl_phi = -99;
        
        if(matchedGenEl_idx >= 0)
        {
            matchedGenEl_energy = v_genEl_4mom.at(matchedGenEl_idx).e();
            matchedGenEl_pT = v_genEl_4mom.at(matchedGenEl_idx).perp();
            matchedGenEl_eta = v_genEl_4mom.at(matchedGenEl_idx).eta();
            matchedGenEl_phi = v_genEl_4mom.at(matchedGenEl_idx).phi();
        }
        
        treeOutput->v_gsfEleFromTICL_matchedGenEl_E.push_back(matchedGenEl_energy);
        treeOutput->v_gsfEleFromTICL_matchedGenEl_pT.push_back(matchedGenEl_pT);
        treeOutput->v_gsfEleFromTICL_matchedGenEl_eta.push_back(matchedGenEl_eta);
        treeOutput->v_gsfEleFromTICL_matchedGenEl_phi.push_back(matchedGenEl_phi);
        
        
        treeOutput->v_gsfEleFromTICL_E.push_back(gsfEle.energy());
        treeOutput->v_gsfEleFromTICL_px.push_back(gsfEle.px());
        treeOutput->v_gsfEleFromTICL_py.push_back(gsfEle.py());
        treeOutput->v_gsfEleFromTICL_pz.push_back(gsfEle.pz());
        
        treeOutput->v_gsfEleFromTICL_pT.push_back(gsfEle.pt());
        treeOutput->v_gsfEleFromTICL_eta.push_back(gsfEle.eta());
        treeOutput->v_gsfEleFromTICL_phi.push_back(gsfEle.phi());
        
        treeOutput->v_gsfEleFromTICL_ET.push_back(gsfEle.et());
        
        treeOutput->gsfEleFromTICL_n++;
        
        treeOutput->v_gsfEleFromTICL_dr03TkSumPt.push_back(gsfEle.dr03TkSumPt());
        
        math::XYZPoint superClus_xyz = gsfEle.superCluster()->position();
        
        treeOutput->v_gsfEleFromTICL_superClus_E.push_back(gsfEle.superCluster()->energy());
        treeOutput->v_gsfEleFromTICL_superClus_ET.push_back(gsfEle.superCluster()->energy() * std::sin(superClus_xyz.theta()));
        treeOutput->v_gsfEleFromTICL_superClus_rawE.push_back(gsfEle.superCluster()->rawEnergy());
        treeOutput->v_gsfEleFromTICL_superClus_rawET.push_back(gsfEle.superCluster()->rawEnergy() * std::sin(superClus_xyz.theta()));
        
        treeOutput->v_gsfEleFromTICL_superClus_theta.push_back(superClus_xyz.theta());
        treeOutput->v_gsfEleFromTICL_superClus_eta.push_back(gsfEle.superCluster()->eta());
        treeOutput->v_gsfEleFromTICL_superClus_phi.push_back(gsfEle.superCluster()->phi());
        treeOutput->v_gsfEleFromTICL_superClus_x.push_back(superClus_xyz.x());
        treeOutput->v_gsfEleFromTICL_superClus_y.push_back(superClus_xyz.y());
        treeOutput->v_gsfEleFromTICL_superClus_z.push_back(superClus_xyz.z());
        treeOutput->v_gsfEleFromTICL_superClus_r.push_back(std::sqrt(superClus_xyz.mag2()));
        
        treeOutput->v_gsfEleFromTICL_superClus_etaWidth.push_back(gsfEle.superCluster()->etaWidth());
        treeOutput->v_gsfEleFromTICL_superClus_phiWidth.push_back(gsfEle.superCluster()->phiWidth());
        
        treeOutput->v_gsfEleFromTICL_superClus_nClus.push_back(gsfEle.superCluster()->clusters().size());
        treeOutput->v_gsfEleFromTICL_superClus_nHit.push_back(gsfEle.superCluster()->hitsAndFractions().size());
        
        std::pair <double, const edm::Ptr <reco::CaloCluster> > p_clusMaxDR = Common::getSCclusterMaxDR2(&(*(gsfEle.superCluster())));
        treeOutput->v_gsfEleFromTICL_superClus_clusMaxDR.push_back(std::sqrt(p_clusMaxDR.first));
        
        std::vector <std::pair <DetId, float> > v_superClus_HandF = gsfEle.superCluster()->hitsAndFractions();
        
        std::vector <std::vector <std::pair <DetId, float> > > vv_superClus_layerHandF = Common::getLayerwiseHandF(
            v_superClus_HandF,
            &recHitTools,
            Constants::HGCalEE_nLayer
        );
        
        printf(
            "[%llu] "
            
            "gsfEleFromTICL %d/%d: "
            "E %0.8f, "
            "pT %0.8f, "
            "eta %+0.8f, "
            "SCrawE %0.10f, "
            "SCeta %0.8f,"
            "\n",
            
            eventNumber,
            
            iEle+1, nEleFromTICL,
            gsfEle.energy(),
            gsfEle.pt(),
            gsfEle.eta(),
            gsfEle.superCluster()->rawEnergy(),
            gsfEle.superCluster()->eta()
        );
        
        // Gsf-track
        treeOutput->v_gsfEleFromTICL_gsfTrack_p.push_back(gsfEle.gsfTrack()->p());
        treeOutput->v_gsfEleFromTICL_gsfTrack_px.push_back(gsfEle.gsfTrack()->px());
        treeOutput->v_gsfEleFromTICL_gsfTrack_py.push_back(gsfEle.gsfTrack()->py());
        treeOutput->v_gsfEleFromTICL_gsfTrack_pz.push_back(gsfEle.gsfTrack()->pz());
        treeOutput->v_gsfEleFromTICL_gsfTrack_pT.push_back(gsfEle.gsfTrack()->pt());
        treeOutput->v_gsfEleFromTICL_gsfTrack_eta.push_back(gsfEle.gsfTrack()->pt());
        treeOutput->v_gsfEleFromTICL_gsfTrack_phi.push_back(gsfEle.gsfTrack()->eta());
        
        // Track at vertex
        treeOutput->v_gsfEleFromTICL_trkAtVtx_p.push_back(gsfEle.trackMomentumAtVtx().r());
        treeOutput->v_gsfEleFromTICL_trkAtVtx_px.push_back(gsfEle.trackMomentumAtVtx().x());
        treeOutput->v_gsfEleFromTICL_trkAtVtx_py.push_back(gsfEle.trackMomentumAtVtx().y());
        treeOutput->v_gsfEleFromTICL_trkAtVtx_pz.push_back(gsfEle.trackMomentumAtVtx().z());
        treeOutput->v_gsfEleFromTICL_trkAtVtx_pT.push_back(gsfEle.trackMomentumAtVtx().rho());
        treeOutput->v_gsfEleFromTICL_trkAtVtx_eta.push_back(gsfEle.trackMomentumAtVtx().eta());
        treeOutput->v_gsfEleFromTICL_trkAtVtx_phi.push_back(gsfEle.trackMomentumAtVtx().phi());
        
        for(auto &key : v_gsfEleFromTICLvar)
        {
            if(m_gsfEleFromTICLvarMap->find(key) == m_gsfEleFromTICLvarMap->emptyRange())
            {
                printf("Error: Cannot find key \"%s\" in gsfEleFromTICL variable map. \n", key.c_str());
                exit(EXIT_FAILURE);
            }
            
            double val = m_gsfEleFromTICLvarMap->find(key)[iEle];
            //printf("Storing \"%s\": %0.4e. \n", key.c_str(), val);
            treeOutput->m_customVarContent.at("gsfEleFromTICL_"+key).push_back(val);
        }
        
        
        const reco::CaloCluster *superClus_seed = gsfEle.superCluster()->seed().get();
        
        CLHEP::Hep3Vector superClus_seed_3vec(
            superClus_seed->x(),
            superClus_seed->y(),
            superClus_seed->z()
        );
        
        treeOutput->v_gsfEleFromTICL_superClusSeed_E.push_back(superClus_seed->energy());
        treeOutput->v_gsfEleFromTICL_superClusSeed_ET.push_back(superClus_seed->energy() * sin(superClus_seed_3vec.theta()));
        treeOutput->v_gsfEleFromTICL_superClusSeed_eta.push_back(superClus_seed->eta());
        treeOutput->v_gsfEleFromTICL_superClusSeed_phi.push_back(superClus_seed->phi());
        
        
        // Sorted SC seed hit energies
        std::vector <std::pair <int, double> > vp_gsfEleFromTICL_superClusSeed_sortedRecHit_index_energy = Common::getSortedHitIndex(
            superClus_seed->hitsAndFractions(),
            m_HGCRecHit
        );
        
        double gsfEleFromTICL_superClusSeed_recHit1_E = 0;
        double gsfEleFromTICL_superClusSeed_recHit2_E = 0;
        
        if(vp_gsfEleFromTICL_superClusSeed_sortedRecHit_index_energy.size() >= 1)
        {
            gsfEleFromTICL_superClusSeed_recHit1_E = vp_gsfEleFromTICL_superClusSeed_sortedRecHit_index_energy.at(0).second;
        }
        
        if(vp_gsfEleFromTICL_superClusSeed_sortedRecHit_index_energy.size() >= 2)
        {
            gsfEleFromTICL_superClusSeed_recHit2_E = vp_gsfEleFromTICL_superClusSeed_sortedRecHit_index_energy.at(1).second;
        }
        
        treeOutput->v_gsfEleFromTICL_superClusSeed_recHit1_E.push_back(gsfEleFromTICL_superClusSeed_recHit1_E);
        treeOutput->v_gsfEleFromTICL_superClusSeed_recHit2_E.push_back(gsfEleFromTICL_superClusSeed_recHit2_E);
    }
    
    
    
    /////////////////////////////////////////////////////////////////
    //////////////////// TICL photons ///////////////////////////////
    /////////////////////////////////////////////////////////////////
    edm::Handle <std::vector <reco::Photon> > v_phoFromTICL;
    iEvent.getByToken(tok_phoFromTICL, v_phoFromTICL);
    
    edm::Handle <edm::MapOfVectors <std::string, double> > m_phoFromTICLvarMap;
    iEvent.getByToken(tok_phoFromTICLvarMap, m_phoFromTICLvarMap);
    
    int nPhoFromTICL = v_phoFromTICL->size();
    
    std::vector <CLHEP::HepLorentzVector> v_phoFromTICL_4mom;
    
    
    for(int iPho = 0; iPho < nPhoFromTICL; iPho++)
    {
        reco::Photon pho = v_phoFromTICL->at(iPho);
        
        CLHEP::HepLorentzVector phoFromTICL_4mom;
        phoFromTICL_4mom.setT(pho.energy());
        phoFromTICL_4mom.setX(pho.px());
        phoFromTICL_4mom.setY(pho.py());
        phoFromTICL_4mom.setZ(pho.pz());
        
        v_phoFromTICL_4mom.push_back(phoFromTICL_4mom);
    }
    
    
    // TICL-ele gen-matching
    //TMatrixD mat_phoFromTICL_genPh_deltaR;
    //
    //std::vector <int> v_phoFromTICL_matchedGenPh_idx;
    //
    //std::vector <double> v_phoFromTICL_genPh_minDeltaR = Common::getMinDeltaR(
    //    v_phoFromTICL_4mom,
    //    v_genPh_4mom,
    //    mat_phoFromTICL_genPh_deltaR,
    //    v_phoFromTICL_matchedGenPh_idx
    //);
    
    std::vector <int> v_phoFromTICL_matchedGenPh_idx;
    std::vector <double> v_phoFromTICL_matchedGenPh_deltaR;
    
    Common::getHardestInCone(
        v_phoFromTICL_4mom,
        v_genPh_4mom,
        v_phoFromTICL_matchedGenPh_idx,
        v_phoFromTICL_matchedGenPh_deltaR,
        0.3
    );
    
    
    for(int iPho = 0; iPho < nPhoFromTICL; iPho++)
    {
        reco::Photon pho = v_phoFromTICL->at(iPho);
        CLHEP::HepLorentzVector phoFromTICL_4mom = v_phoFromTICL_4mom.at(iPho);
        
        if(pho.pt() < el_minPt || fabs(pho.eta()) < HGCal_minEta || fabs(pho.eta()) > HGCal_maxEta)
        {
            continue;
        }
        
        //double matchedGenPh_deltaR = v_phoFromTICL_genPh_minDeltaR.at(iPho);
        double matchedGenPh_deltaR = v_phoFromTICL_matchedGenPh_deltaR.at(iPho);
        
        if(matchedGenPh_deltaR > TICLphoGenMatchDR)
        {
            continue;
        }
        
        printf(
            "[%llu] "
            
            "phoFromTICL %d/%d: "
            "E %0.4f, "
            "pT %0.2f, "
            "eta %+0.2f, "
            "SC E %+0.4f, "
            "SC raw E %+0.4f, "
            
            "\n",
            
            eventNumber,
            
            iPho+1, nPhoFromTICL,
            pho.energy(),
            pho.pt(),
            pho.eta(),
            pho.superCluster()->energy(),
            pho.superCluster()->rawEnergy()
        );
        
        int matchedGenPh_idx = v_phoFromTICL_matchedGenPh_idx.at(iPho);
        
        treeOutput->v_phoFromTICL_genPh_minDeltaR.push_back(matchedGenPh_deltaR);
        treeOutput->v_phoFromTICL_nearestGenPh_idx.push_back(matchedGenPh_idx);
        
        double matchedGenPh_energy = -99;
        double matchedGenPh_pT = -99;
        double matchedGenPh_eta = -99;
        double matchedGenPh_phi = -99;
        
        if(matchedGenPh_idx >= 0)
        {
            matchedGenPh_energy = v_genPh_4mom.at(matchedGenPh_idx).e();
            matchedGenPh_pT = v_genPh_4mom.at(matchedGenPh_idx).perp();
            matchedGenPh_eta = v_genPh_4mom.at(matchedGenPh_idx).eta();
            matchedGenPh_phi = v_genPh_4mom.at(matchedGenPh_idx).phi();
        }
        
        treeOutput->v_phoFromTICL_matchedGenPh_E.push_back(matchedGenPh_energy);
        treeOutput->v_phoFromTICL_matchedGenPh_pT.push_back(matchedGenPh_pT);
        treeOutput->v_phoFromTICL_matchedGenPh_eta.push_back(matchedGenPh_eta);
        treeOutput->v_phoFromTICL_matchedGenPh_phi.push_back(matchedGenPh_phi);
        
        
        treeOutput->v_phoFromTICL_E.push_back(pho.energy());
        treeOutput->v_phoFromTICL_px.push_back(pho.px());
        treeOutput->v_phoFromTICL_py.push_back(pho.py());
        treeOutput->v_phoFromTICL_pz.push_back(pho.pz());
        
        treeOutput->v_phoFromTICL_pT.push_back(pho.pt());
        treeOutput->v_phoFromTICL_eta.push_back(pho.eta());
        treeOutput->v_phoFromTICL_phi.push_back(pho.phi());
        
        treeOutput->v_phoFromTICL_ET.push_back(pho.et());
        
        treeOutput->phoFromTICL_n++;
        
        
        math::XYZPoint superClus_xyz = pho.superCluster()->position();
        
        treeOutput->v_phoFromTICL_superClus_E.push_back(pho.superCluster()->energy());
        treeOutput->v_phoFromTICL_superClus_ET .push_back(pho.superCluster()->energy() * std::sin(superClus_xyz.theta()));
        treeOutput->v_phoFromTICL_superClus_rawE.push_back(pho.superCluster()->rawEnergy());
        treeOutput->v_phoFromTICL_superClus_rawET.push_back(pho.superCluster()->rawEnergy() * std::sin(superClus_xyz.theta()));
        
        treeOutput->v_phoFromTICL_superClus_theta.push_back(superClus_xyz.theta());
        treeOutput->v_phoFromTICL_superClus_eta.push_back(pho.superCluster()->eta());
        treeOutput->v_phoFromTICL_superClus_phi.push_back(pho.superCluster()->phi());
        treeOutput->v_phoFromTICL_superClus_x.push_back(superClus_xyz.x());
        treeOutput->v_phoFromTICL_superClus_y.push_back(superClus_xyz.y());
        treeOutput->v_phoFromTICL_superClus_z.push_back(superClus_xyz.z());
        treeOutput->v_phoFromTICL_superClus_r.push_back(std::sqrt(superClus_xyz.mag2()));
        
        treeOutput->v_phoFromTICL_superClus_etaWidth.push_back(pho.superCluster()->etaWidth());
        treeOutput->v_phoFromTICL_superClus_phiWidth.push_back(pho.superCluster()->phiWidth());
        
        treeOutput->v_phoFromTICL_superClus_nClus.push_back(pho.superCluster()->clusters().size());
        treeOutput->v_phoFromTICL_superClus_nHit.push_back(pho.superCluster()->hitsAndFractions().size());
        
        std::pair <double, const edm::Ptr <reco::CaloCluster> > p_clusMaxDR = Common::getSCclusterMaxDR2(&(*(pho.superCluster())));
        treeOutput->v_phoFromTICL_superClus_clusMaxDR.push_back(std::sqrt(p_clusMaxDR.first));
        
        std::vector <std::pair <DetId, float> > v_superClus_HandF = pho.superCluster()->hitsAndFractions();
        //math::XYZPoint superClus_xyz = pho.superCluster()->position();
        
        treeOutput->v_phoFromTICL_superClus_seed_dEta.push_back(pho.superCluster()->eta() - pho.superCluster()->seed().get()->eta());
        treeOutput->v_phoFromTICL_superClus_seed_dPhi.push_back(TVector2::Phi_mpi_pi(pho.superCluster()->phi() - pho.superCluster()->seed().get()->phi()));
        
        std::vector <double> v_phoFromTICL_recHit_E;
        std::vector <double> v_phoFromTICL_recHit_x;
        std::vector <double> v_phoFromTICL_recHit_y;
        std::vector <double> v_phoFromTICL_recHit_z;
        std::vector <double> v_phoFromTICL_recHit_time;
        std::vector <double> v_phoFromTICL_recHit_timeError;
        std::vector <double> v_phoFromTICL_recHit_eta;
        std::vector <double> v_phoFromTICL_recHit_phi;
        std::vector <double> v_phoFromTICL_recHit_ET;
        std::vector <double> v_phoFromTICL_recHit_detector;
        std::vector <double> v_phoFromTICL_recHit_layer;
        std::vector <double> v_phoFromTICL_recHit_isSimHitMatched;
        std::vector <double> v_phoFromTICL_recHit_SCdEta;
        std::vector <double> v_phoFromTICL_recHit_SCdPhi;
        std::vector <double> v_phoFromTICL_recHit_SCdR;
        
        int nHit = 0;
        
        for(const auto &hnf : v_superClus_HandF)
        {
            DetId hitId = hnf.first;
            double hitEfrac = hnf.second;
            
            int hitLayer = recHitTools.getLayer(hitId);
            
            //if (hitLayer > nLayer_)
            //{
            //    continue;
            //}
        
            //if (hitId.det() != subDet_)
            //{
            //    continue;
            //}
            
            auto pfHitIt = m_PFRecHitHGC.find(hitId);
            
            if (pfHitIt == m_PFRecHitHGC.end())
            {
                continue;
            }
            
            const reco::PFRecHit &pfRecHit = *(pfHitIt->second);
            
            auto hgcHitIt = m_HGCEERecHit.find(hitId);
            
            if (hgcHitIt == m_HGCEERecHit.end())
            {
                continue;
            }
            
            const HGCRecHit &hgcRecHit = *(hgcHitIt->second);
            
            nHit++;
            
            //printf("Photon SC PFRecHit: E %0.4f \n", pfRecHit.energy());
            
            v_phoFromTICL_recHit_E.push_back(pfRecHit.energy());
            v_phoFromTICL_recHit_x.push_back(pfRecHit.position().x());
            v_phoFromTICL_recHit_y.push_back(pfRecHit.position().y());
            v_phoFromTICL_recHit_z.push_back(pfRecHit.position().z());
            v_phoFromTICL_recHit_time.push_back(hgcRecHit.time());
            v_phoFromTICL_recHit_timeError.push_back(hgcRecHit.timeError());
            v_phoFromTICL_recHit_eta.push_back(pfRecHit.positionREP().eta());
            v_phoFromTICL_recHit_phi.push_back(pfRecHit.positionREP().phi());
            v_phoFromTICL_recHit_ET.push_back(pfRecHit.energy() * std::sin(pfRecHit.position().theta()));
            v_phoFromTICL_recHit_detector.push_back(hitId.det());
            v_phoFromTICL_recHit_layer.push_back(recHitTools.getLayer(hitId) - 1); // Start from 0
            v_phoFromTICL_recHit_isSimHitMatched.push_back(m_simHit.find(hitId) != m_simHit.end());
            v_phoFromTICL_recHit_SCdEta.push_back(pfRecHit.positionREP().eta() - pho.superCluster()->eta());
            v_phoFromTICL_recHit_SCdPhi.push_back(reco::deltaPhi(pfRecHit.positionREP().phi(), pho.superCluster()->phi()));
            v_phoFromTICL_recHit_SCdR.push_back(reco::deltaR(pfRecHit.positionREP(), *pho.superCluster()));
        }
        
        treeOutput->v_phoFromTICL_nRecHit.push_back(nHit);
        treeOutput->vv_phoFromTICL_recHit_E.push_back(v_phoFromTICL_recHit_E);
        treeOutput->vv_phoFromTICL_recHit_x.push_back(v_phoFromTICL_recHit_x);
        treeOutput->vv_phoFromTICL_recHit_y.push_back(v_phoFromTICL_recHit_y);
        treeOutput->vv_phoFromTICL_recHit_z.push_back(v_phoFromTICL_recHit_z);
        treeOutput->vv_phoFromTICL_recHit_time.push_back(v_phoFromTICL_recHit_time);
        treeOutput->vv_phoFromTICL_recHit_timeError.push_back(v_phoFromTICL_recHit_timeError);
        treeOutput->vv_phoFromTICL_recHit_eta.push_back(v_phoFromTICL_recHit_eta);
        treeOutput->vv_phoFromTICL_recHit_phi.push_back(v_phoFromTICL_recHit_phi);
        treeOutput->vv_phoFromTICL_recHit_ET.push_back(v_phoFromTICL_recHit_ET);
        treeOutput->vv_phoFromTICL_recHit_detector.push_back(v_phoFromTICL_recHit_detector);
        treeOutput->vv_phoFromTICL_recHit_layer.push_back(v_phoFromTICL_recHit_layer);
        treeOutput->vv_phoFromTICL_recHit_isSimHitMatched.push_back(v_phoFromTICL_recHit_isSimHitMatched);
        treeOutput->vv_phoFromTICL_recHit_SCdEta.push_back(v_phoFromTICL_recHit_SCdEta);
        treeOutput->vv_phoFromTICL_recHit_SCdPhi.push_back(v_phoFromTICL_recHit_SCdPhi);
        treeOutput->vv_phoFromTICL_recHit_SCdR.push_back(v_phoFromTICL_recHit_SCdR);
        
        
        // Layer clusters in cone
        int lcIdx = -1;
        int lcCount = 0;
        std::vector <double> v_phoFromTICL_lcInCone_E;
        std::vector <double> v_phoFromTICL_lcInCone_x;
        std::vector <double> v_phoFromTICL_lcInCone_y;
        std::vector <double> v_phoFromTICL_lcInCone_z;
        std::vector <double> v_phoFromTICL_lcInCone_time;
        std::vector <double> v_phoFromTICL_lcInCone_timeError;
        std::vector <double> v_phoFromTICL_lcInCone_eta;
        std::vector <double> v_phoFromTICL_lcInCone_phi;
        std::vector <double> v_phoFromTICL_lcInCone_ET;
        std::vector <double> v_phoFromTICL_lcInCone_size;
        std::vector <double> v_phoFromTICL_lcInCone_detector;
        std::vector <double> v_phoFromTICL_lcInCone_layer;
        std::vector <double> v_phoFromTICL_lcInCone_SCdEta;
        std::vector <double> v_phoFromTICL_lcInCone_SCdPhi;
        std::vector <double> v_phoFromTICL_lcInCone_SCdR;
        
        double maxDR2 = 0.3*0.3;
        
        for(const auto lc : (*v_HGCALlayerCluster))
        {
            lcIdx++;
            
            edm::Ref <reco::CaloClusterCollection> lcRef(v_HGCALlayerCluster, lcIdx);
            std::pair<float,float> timeAndErr = (*vm_HGCALlayerClusterTime)[lcRef];
            
            if(reco::deltaR2(pho.eta(), pho.phi(), lc.eta(), lc.phi()) > maxDR2)
            {
                continue;
            }
            
            lcCount++;
            
            //printf(
            //    "[phoFromTICL %d/%d, HGCAL layer cluster %d] "
            //    "E %f, size %d, time %f %f "
            //    "\n",
            //    iPho+1, nPhoFromTICL, lcCount,
            //    lc.energy(), (int) lc.size(), timeAndErr.first, timeAndErr.second
            //);
            
            v_phoFromTICL_lcInCone_E.push_back(lc.energy());
            v_phoFromTICL_lcInCone_x.push_back(lc.position().x());
            v_phoFromTICL_lcInCone_y.push_back(lc.position().y());
            v_phoFromTICL_lcInCone_z.push_back(lc.position().z());
            v_phoFromTICL_lcInCone_time.push_back(timeAndErr.first);
            v_phoFromTICL_lcInCone_timeError.push_back(timeAndErr.second);
            v_phoFromTICL_lcInCone_eta.push_back(lc.eta());
            v_phoFromTICL_lcInCone_phi.push_back(lc.phi());
            v_phoFromTICL_lcInCone_ET.push_back(lc.energy() * std::sin(lc.position().theta()));
            v_phoFromTICL_lcInCone_size.push_back(lc.size());
            v_phoFromTICL_lcInCone_detector.push_back(lc.seed().det());
            v_phoFromTICL_lcInCone_layer.push_back(recHitTools.getLayer(lc.seed()) - 1); // Start from 0
            v_phoFromTICL_lcInCone_SCdEta.push_back(lc.eta() - pho.superCluster()->eta());
            v_phoFromTICL_lcInCone_SCdPhi.push_back(reco::deltaPhi(lc.phi(), pho.superCluster()->phi()));
            v_phoFromTICL_lcInCone_SCdR.push_back(reco::deltaR(lc, *pho.superCluster()));
        }
        
        treeOutput->v_phoFromTICL_nlcInCone.push_back(lcCount);
        treeOutput->vv_phoFromTICL_lcInCone_E.push_back(v_phoFromTICL_lcInCone_E);
        treeOutput->vv_phoFromTICL_lcInCone_x.push_back(v_phoFromTICL_lcInCone_x);
        treeOutput->vv_phoFromTICL_lcInCone_y.push_back(v_phoFromTICL_lcInCone_y);
        treeOutput->vv_phoFromTICL_lcInCone_z.push_back(v_phoFromTICL_lcInCone_z);
        treeOutput->vv_phoFromTICL_lcInCone_time.push_back(v_phoFromTICL_lcInCone_time);
        treeOutput->vv_phoFromTICL_lcInCone_timeError.push_back(v_phoFromTICL_lcInCone_timeError);
        treeOutput->vv_phoFromTICL_lcInCone_eta.push_back(v_phoFromTICL_lcInCone_eta);
        treeOutput->vv_phoFromTICL_lcInCone_phi.push_back(v_phoFromTICL_lcInCone_phi);
        treeOutput->vv_phoFromTICL_lcInCone_ET.push_back(v_phoFromTICL_lcInCone_ET);
        treeOutput->vv_phoFromTICL_lcInCone_size.push_back(v_phoFromTICL_lcInCone_size);
        treeOutput->vv_phoFromTICL_lcInCone_detector.push_back(v_phoFromTICL_lcInCone_detector);
        treeOutput->vv_phoFromTICL_lcInCone_layer.push_back(v_phoFromTICL_lcInCone_layer);
        treeOutput->vv_phoFromTICL_lcInCone_SCdEta.push_back(v_phoFromTICL_lcInCone_SCdEta);
        treeOutput->vv_phoFromTICL_lcInCone_SCdPhi.push_back(v_phoFromTICL_lcInCone_SCdPhi);
        treeOutput->vv_phoFromTICL_lcInCone_SCdR.push_back(v_phoFromTICL_lcInCone_SCdR);
        
        
        for(auto &key : v_phoFromTICLvar)
        {
            if(m_phoFromTICLvarMap->find(key) == m_phoFromTICLvarMap->emptyRange())
            {
                printf("Error: Cannot find key \"%s\" in phoFromTICL variable map. \n", key.c_str());
                exit(EXIT_FAILURE);
            }
            
            double val = m_phoFromTICLvarMap->find(key)[iPho];
            //printf("Storing \"%s\": %0.4e. \n", key.c_str(), val);
            treeOutput->m_customVarContent.at("phoFromTICL_"+key).push_back(val);
        }
        
        
        const reco::CaloCluster *superClus_seed = pho.superCluster()->seed().get();
        
        CLHEP::Hep3Vector superClus_seed_3vec(
            superClus_seed->x(),
            superClus_seed->y(),
            superClus_seed->z()
        );
        
        treeOutput->v_phoFromTICL_superClusSeed_E.push_back(superClus_seed->energy());
        treeOutput->v_phoFromTICL_superClusSeed_ET.push_back(superClus_seed->energy() * sin(superClus_seed_3vec.theta()));
        treeOutput->v_phoFromTICL_superClusSeed_eta.push_back(superClus_seed->eta());
        treeOutput->v_phoFromTICL_superClusSeed_phi.push_back(superClus_seed->phi());
    }
    
    
    // Rechit counts for pileup
    treeOutput->nHit_EcalEB = Common::countHits(*v_EcalEBRecHit, 1.0);
    treeOutput->nHit_HGCEE  = Common::countHits(*v_HGCEERecHit , 1.0);
    treeOutput->nHit_HGCHEF = Common::countHits(*v_HGCHEFRecHit, 1.0);
    treeOutput->nHit_HGCHEB = Common::countHits(*v_HGCHEBRecHit, 1.0);
    
    
    // Fill tree
    treeOutput->fill();
    
    //#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
    //ESHandle<SetupData> pSetup;
    //iSetup.get<SetupRecord>().get(pSetup);
    //#endif
    
    printf("\n\n");
    
    fflush(stdout);
    fflush(stderr);
}


// ------------ method called once each job just before starting event loop  ------------
void
TreeMaker::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
TreeMaker::endJob()
{
    
    printf("minLayer = %d, maxLayer = %d \n", minLayer, maxLayer);
    
    
    fflush(stdout);
    fflush(stderr);
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TreeMaker::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TreeMaker);
