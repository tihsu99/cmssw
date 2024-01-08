# ifndef TreeOutputInfo_H
# define TreeOutputInfo_H


# include <iostream>
# include <map>
# include <stdlib.h>
# include <string>
# include <type_traits>
# include <utility>
# include <vector>

# include <TH1F.h>
# include <TH2F.h>
# include <TMatrixD.h>
# include <TROOT.h>
# include <TTree.h> 
# include <TVectorD.h> 

# include "EDAnalyzers/TreeMaker/interface/Constants.h"


namespace TreeOutputInfo
{
    class TreeOutput
    {
        public :
        
        
        TTree *tree;
        
        
        // Run info //
        ULong64_t runNumber;
        ULong64_t eventNumber;
        ULong64_t luminosityNumber;
        ULong64_t bunchCrossingNumber;
        
        
        // Gen event info //
        double genEventWeight;
        
        
        // Gen electron //
        int genEl_n;
        std::vector <double> v_genEl_E;
        std::vector <double> v_genEl_px;
        std::vector <double> v_genEl_py;
        std::vector <double> v_genEl_pz;
        std::vector <double> v_genEl_pT;
        std::vector <double> v_genEl_eta;
        std::vector <double> v_genEl_phi;
        
        std::vector <double> v_genEl_HGCalEEP_EsortedIndex;
        std::vector <double> v_genEl_HGCalEEM_EsortedIndex;
        
        std::vector <double> v_genEl_multiClus_totE;
        
        std::vector <double> v_genEl_multiClus_n;
        std::vector <double> v_genEl_nearestMultiClusEnRatio;
        std::vector <double> v_genEl_multiClusEnRatio;
        
        
        // Gen photon //
        int genPh_n;
        std::vector <double> v_genPh_E;
        std::vector <double> v_genPh_px;
        std::vector <double> v_genPh_py;
        std::vector <double> v_genPh_pz;
        std::vector <double> v_genPh_pT;
        std::vector <double> v_genPh_eta;
        std::vector <double> v_genPh_phi;
        
        std::vector <double> v_genPh_HGCalEEP_EsortedIndex;
        std::vector <double> v_genPh_HGCalEEM_EsortedIndex;
        
        std::vector <double> v_genPh_multiClus_totE;
        
        std::vector <double> v_genPh_HGCalEEP_deltaR;
        std::vector <double> v_genPh_HGCalEEM_deltaR;
        
        
        // Pileup //
        int pileup_n;
        
        
        // Rho //
        double rho;
        
        
        // Hit counts (pileup estimator) //
        double nHit_EcalEB;
        double nHit_HGCEE;
        double nHit_HGCHEF;
        double nHit_HGCHEB;
        
        
        // HGCAL layer clusters //
        int HGCALlayerClus_n;
        std::vector <double> v_HGCALlayerClus_E;
        std::vector <double> v_HGCALlayerClus_x;
        std::vector <double> v_HGCALlayerClus_y;
        std::vector <double> v_HGCALlayerClus_z;
        std::vector <double> v_HGCALlayerClus_eta;
        std::vector <double> v_HGCALlayerClus_phi;
        std::vector <double> v_HGCALlayerClus_ET;
        std::vector <double> v_HGCALlayerClus_time;
        std::vector <double> v_HGCALlayerClus_timeError;
        std::vector <double> v_HGCALlayerClus_detector;
        std::vector <double> v_HGCALlayerClus_layer;
        
        
        // Tracksters //
        int trackster_n;
        std::vector <double> v_trackster_E;
        std::vector <double> v_trackster_x;
        std::vector <double> v_trackster_y;
        std::vector <double> v_trackster_z;
        std::vector <double> v_trackster_eta;
        std::vector <double> v_trackster_phi;
        std::vector <double> v_trackster_ET;
        
        
        // MultiClusters //
        int multiClus_n;
        std::vector <double> v_multiClus_genElIndex;
        std::vector <double> v_multiClus_E;
        std::vector <double> v_multiClus_x;
        std::vector <double> v_multiClus_y;
        std::vector <double> v_multiClus_z;
        std::vector <double> v_multiClus_eta;
        std::vector <double> v_multiClus_phi;
        std::vector <double> v_multiClus_ET;
        
        std::vector <double> v_multiClus_corrE;
        std::vector <double> v_multiClus_corrET;
        
        std::vector <double> v_multiClus_dX;
        std::vector <double> v_multiClus_dY;
        std::vector <double> v_multiClus_dZ;
        
        std::vector <double> v_multiClus_dEta;
        std::vector <double> v_multiClus_dPhi;
        
        std::vector <double> v_multiClus_sigma2rr;
        std::vector <double> v_multiClus_sigma2etaEta;
        std::vector <double> v_multiClus_sigma2phiPhi;
        
        std::vector <double> v_multiClus_sigma2rEta;
        std::vector <double> v_multiClus_sigma2rPhi;
        std::vector <double> v_multiClus_sigma2etaPhi;
        
        std::vector <double> v_multiClus_sigma2diag1;
        std::vector <double> v_multiClus_sigma2diag2;
        std::vector <double> v_multiClus_sigma2diag3;
        
        std::vector <double> v_multiClus_EsortedIndex;
        std::vector <double> v_multiClus_HGCalEEP_EsortedIndex;
        std::vector <double> v_multiClus_HGCalEEM_EsortedIndex;
        
        std::vector <double> v_multiClus_clus_n;
        std::vector <double> v_multiClus_clus_startIndex;
        std::vector <double> v_multiClus_clus_E;
        std::vector <double> v_multiClus_clus_x;
        std::vector <double> v_multiClus_clus_y;
        std::vector <double> v_multiClus_clus_z;
        std::vector <double> v_multiClus_clus_eta;
        std::vector <double> v_multiClus_clus_phi;
        std::vector <double> v_multiClus_clus_ET;
        
        std::vector <double> v_multiClus_clus_layer;
        std::vector <double> v_multiClus_clus_multiplicity;
        std::vector <double> v_multiClus_clus_nHit;
        
        
        //
        std::vector <double> v_multiClus_uniqueClus_n;
        std::vector <double> v_multiClus_uniqueClus_E;
        std::vector <double> v_multiClus_uniqueClus_x;
        std::vector <double> v_multiClus_uniqueClus_y;
        std::vector <double> v_multiClus_uniqueClus_z;
        std::vector <double> v_multiClus_uniqueClus_eta;
        std::vector <double> v_multiClus_uniqueClus_phi;
        std::vector <double> v_multiClus_uniqueClus_ET;
        
        std::vector <double> v_multiClus_uniqueClus_layer;
        std::vector <double> v_multiClus_uniqueClus_multiplicity;
        std::vector <double> v_multiClus_uniqueClus_nHit;
        
        
        //
        std::vector <double> v_multiClus_mc1_dX;
        std::vector <double> v_multiClus_mc1_dY;
        std::vector <double> v_multiClus_mc1_dZ;
        
        std::vector <double> v_multiClus_mc1_dEta;
        std::vector <double> v_multiClus_mc1_dPhi;
        std::vector <double> v_multiClus_mc1_dR;
        
        
        //
        double multiClus_HGCalEEP_meanX;
        double multiClus_HGCalEEP_meanY;
        double multiClus_HGCalEEP_meanZ;
        
        double multiClus_HGCalEEP_meanDx;
        double multiClus_HGCalEEP_meanDy;
        double multiClus_HGCalEEP_meanDz;
        
        double multiClus_HGCalEEP_totE;
        double multiClus_HGCalEEP_totET;
        
        double multiClus_HGCalEEP_diag1;
        double multiClus_HGCalEEP_diag2;
        double multiClus_HGCalEEP_diag3;
        
        //
        double multiClus_HGCalEEM_meanX;
        double multiClus_HGCalEEM_meanY;
        double multiClus_HGCalEEM_meanZ;
        
        double multiClus_HGCalEEM_meanDx;
        double multiClus_HGCalEEM_meanDy;
        double multiClus_HGCalEEM_meanDz;
        
        double multiClus_HGCalEEM_totE;
        double multiClus_HGCalEEM_totET;
        
        double multiClus_HGCalEEM_diag1;
        double multiClus_HGCalEEM_diag2;
        double multiClus_HGCalEEM_diag3;
        
        
        double simHit_n;
        std::vector <double> v_simHit_E;
        std::vector <double> v_simHit_x;
        std::vector <double> v_simHit_y;
        std::vector <double> v_simHit_z;
        std::vector <double> v_simHit_eta;
        std::vector <double> v_simHit_phi;
        std::vector <double> v_simHit_ET;
        std::vector <double> v_simHit_layer;
        std::vector <double> v_simHit_zside;
        std::vector <double> v_simHit_isCaloParticleMatched;
        
        
        double recHit_n;
        std::vector <double> v_recHit_E;
        std::vector <double> v_recHit_x;
        std::vector <double> v_recHit_y;
        std::vector <double> v_recHit_z;
        std::vector <double> v_recHit_eta;
        std::vector <double> v_recHit_phi;
        std::vector <double> v_recHit_ET;
        std::vector <double> v_recHit_layer;
        std::vector <double> v_recHit_zside;
        std::vector <double> v_recHit_matchedSimHitIndex;
        std::vector <double> v_recHit_matchedHGCALlayerClusIndex;
        std::vector <double> v_recHit_isMultiClusMatched;
        std::vector <double> v_recHit_isCaloParticleMatched;
        
        std::vector <double> v_recHit_iType;
        std::vector <double> v_recHit_iCell1;
        std::vector <double> v_recHit_iCell2;
        
        std::vector <double> v_recHit_SiThickness;
        
        std::vector <double> v_simHit_HGCalEEPlayer_totE;
        std::vector <double> v_recHit_HGCalEEPlayer_totE;
        
        std::vector <double> v_simHit_HGCalEEMlayer_totE;
        std::vector <double> v_recHit_HGCalEEMlayer_totE;
        
        
        double gsfEleFromMultiClus_n;
        std::vector <double> v_gsfEleFromMultiClus_E;
        std::vector <double> v_gsfEleFromMultiClus_px;
        std::vector <double> v_gsfEleFromMultiClus_py;
        std::vector <double> v_gsfEleFromMultiClus_pz;
        std::vector <double> v_gsfEleFromMultiClus_pT;
        std::vector <double> v_gsfEleFromMultiClus_eta;
        std::vector <double> v_gsfEleFromMultiClus_phi;
        
        std::vector <double> v_gsfEleFromMultiClus_genEl_minDeltaR;
        std::vector <double> v_gsfEleFromMultiClus_matchedGenEl_E;
        
        
        double gsfEleFromTICL_n;
        std::vector <double> v_gsfEleFromTICL_E;
        std::vector <double> v_gsfEleFromTICL_px;
        std::vector <double> v_gsfEleFromTICL_py;
        std::vector <double> v_gsfEleFromTICL_pz;
        std::vector <double> v_gsfEleFromTICL_pT;
        std::vector <double> v_gsfEleFromTICL_eta;
        std::vector <double> v_gsfEleFromTICL_phi;
        std::vector <double> v_gsfEleFromTICL_ET;
        
        std::vector <double> v_gsfEleFromTICL_genEl_minDeltaR;
        std::vector <double> v_gsfEleFromTICL_nearestGenEl_idx;
        std::vector <double> v_gsfEleFromTICL_matchedGenEl_E;
        std::vector <double> v_gsfEleFromTICL_matchedGenEl_pT;
        std::vector <double> v_gsfEleFromTICL_matchedGenEl_eta;
        std::vector <double> v_gsfEleFromTICL_matchedGenEl_phi;
        
        std::vector <double> v_gsfEleFromTICL_gsfTrack_p;
        std::vector <double> v_gsfEleFromTICL_gsfTrack_px;
        std::vector <double> v_gsfEleFromTICL_gsfTrack_py;
        std::vector <double> v_gsfEleFromTICL_gsfTrack_pz;
        std::vector <double> v_gsfEleFromTICL_gsfTrack_pT;
        std::vector <double> v_gsfEleFromTICL_gsfTrack_eta;
        std::vector <double> v_gsfEleFromTICL_gsfTrack_phi;
        
        std::vector <double> v_gsfEleFromTICL_trkAtVtx_p;
        std::vector <double> v_gsfEleFromTICL_trkAtVtx_px;
        std::vector <double> v_gsfEleFromTICL_trkAtVtx_py;
        std::vector <double> v_gsfEleFromTICL_trkAtVtx_pz;
        std::vector <double> v_gsfEleFromTICL_trkAtVtx_pT;
        std::vector <double> v_gsfEleFromTICL_trkAtVtx_eta;
        std::vector <double> v_gsfEleFromTICL_trkAtVtx_phi;
        
        std::vector <double> v_gsfEleFromTICL_dr03TkSumPt;
        
        std::vector <double> v_gsfEleFromTICL_superClus_E;
        std::vector <double> v_gsfEleFromTICL_superClus_ET;
        std::vector <double> v_gsfEleFromTICL_superClus_rawE;
        std::vector <double> v_gsfEleFromTICL_superClus_rawET;
        std::vector <double> v_gsfEleFromTICL_superClus_theta;
        std::vector <double> v_gsfEleFromTICL_superClus_eta;
        std::vector <double> v_gsfEleFromTICL_superClus_phi;
        std::vector <double> v_gsfEleFromTICL_superClus_x;
        std::vector <double> v_gsfEleFromTICL_superClus_y;
        std::vector <double> v_gsfEleFromTICL_superClus_z;
        std::vector <double> v_gsfEleFromTICL_superClus_r;
        std::vector <double> v_gsfEleFromTICL_superClus_etaWidth;
        std::vector <double> v_gsfEleFromTICL_superClus_phiWidth;
        std::vector <double> v_gsfEleFromTICL_superClus_nClus;
        std::vector <double> v_gsfEleFromTICL_superClus_nHit;
        std::vector <double> v_gsfEleFromTICL_superClus_nearestCellDist;
        std::vector <double> v_gsfEleFromTICL_superClus_cellNeighbour1ringWindow_n;
        std::vector <double> v_gsfEleFromTICL_superClus_cellNeighbour2ringWindow_n;
        std::vector <double> v_gsfEleFromTICL_superClus_clusMaxDR;
        
        std::vector <double> v_gsfEleFromTICL_superClus_seed_dEta;
        std::vector <double> v_gsfEleFromTICL_superClus_seed_dPhi;
        
        std::vector <double> v_gsfEleFromTICL_superClus_recHit1_E;
        std::vector <double> v_gsfEleFromTICL_superClus_recHit2_E;
        
        std::vector <double> v_gsfEleFromTICL_superClus_sigma2etaEtaLW;
        std::vector <double> v_gsfEleFromTICL_superClus_sigma2phiPhiLW;
        
        std::vector <double> v_gsfEleFromTICL_superClus_sigma2rr;
        std::vector <double> v_gsfEleFromTICL_superClus_sigma2etaEta;
        std::vector <double> v_gsfEleFromTICL_superClus_sigma2phiPhi;
        
        std::vector <double> v_gsfEleFromTICL_superClus_sigma2rEta;
        std::vector <double> v_gsfEleFromTICL_superClus_sigma2rPhi;
        std::vector <double> v_gsfEleFromTICL_superClus_sigma2etaPhi;
        
        std::vector <double> v_gsfEleFromTICL_superClus_sigma2diag1;
        std::vector <double> v_gsfEleFromTICL_superClus_sigma2diag2;
        std::vector <double> v_gsfEleFromTICL_superClus_sigma2diag3;
        
        std::vector <std::vector <double> > vv_gsfEleFromTICL_superClus_E_layer;
        
        std::vector <double> v_gsfEleFromTICL_E7;
        std::vector <double> v_gsfEleFromTICL_R7;
        
        std::vector <double> v_gsfEleFromTICL_E19;
        std::vector <double> v_gsfEleFromTICL_R19;
        
        std::vector <double> v_gsfEleFromTICL_R2p0;
        
        std::vector <double> v_gsfEleFromTICL_E2p4;
        std::vector <double> v_gsfEleFromTICL_R2p4;
        std::vector <std::vector <double> > vv_gsfEleFromTICL_E2p4_layer;
        std::vector <std::vector <double> > vv_gsfEleFromTICL_R2p4_layer;
        
        std::vector <double> v_gsfEleFromTICL_E2p6;
        std::vector <double> v_gsfEleFromTICL_R2p6;
        std::vector <std::vector <double> > vv_gsfEleFromTICL_E2p6_layer;
        std::vector <std::vector <double> > vv_gsfEleFromTICL_R2p6_layer;
        
        std::vector <double> v_gsfEleFromTICL_E2p8;
        std::vector <double> v_gsfEleFromTICL_R2p8;
        std::vector <std::vector <double> > vv_gsfEleFromTICL_E2p8_layer;
        std::vector <std::vector <double> > vv_gsfEleFromTICL_R2p8_layer;
        
        std::vector <double> v_gsfEleFromTICL_E3p0;
        std::vector <double> v_gsfEleFromTICL_R3p0;
        std::vector <std::vector <double> > vv_gsfEleFromTICL_E3p0_layer;
        std::vector <std::vector <double> > vv_gsfEleFromTICL_R3p0_layer;
        
        std::vector <double> v_gsfEleFromTICL_E3p5;
        std::vector <double> v_gsfEleFromTICL_R3p5;
        std::vector <std::vector <double> > vv_gsfEleFromTICL_E3p5_layer;
        std::vector <std::vector <double> > vv_gsfEleFromTICL_R3p5_layer;
        
        std::vector <double> v_gsfEleFromTICL_superClusSeed_E;
        std::vector <double> v_gsfEleFromTICL_superClusSeed_ET;
        std::vector <double> v_gsfEleFromTICL_superClusSeed_eta;
        std::vector <double> v_gsfEleFromTICL_superClusSeed_phi;
        
        std::vector <double> v_gsfEleFromTICL_superClusSeed_recHit1_E;
        std::vector <double> v_gsfEleFromTICL_superClusSeed_recHit2_E;
        
        std::vector <std::vector<double> > vv_gsfEleFromTICL_superClusSeed_clus_dEta;
        std::vector <std::vector<double> > vv_gsfEleFromTICL_superClusSeed_clus_dPhi;
        
        std::vector <double> v_gsfEleFromTICL_superClus_TICLclus_n;
        std::vector <std::vector<double> > vv_gsfEleFromTICL_superClus_TICLclus_E;
        std::vector <std::vector<double> > vv_gsfEleFromTICL_superClus_TICLclus_ET;
        std::vector <std::vector<double> > vv_gsfEleFromTICL_superClus_TICLclus_nClus;
        
        std::vector <std::vector<double> > vv_gsfEleFromTICL_superClusSeed_TICLclus_dX;
        std::vector <std::vector<double> > vv_gsfEleFromTICL_superClusSeed_TICLclus_dY;
        std::vector <std::vector<double> > vv_gsfEleFromTICL_superClusSeed_TICLclus_dZ;
        
        std::vector <std::vector<double> > vv_gsfEleFromTICL_superClusSeed_TICLclus_dEta;
        std::vector <std::vector<double> > vv_gsfEleFromTICL_superClusSeed_TICLclus_dPhi;
        std::vector <std::vector<double> > vv_gsfEleFromTICL_superClusSeed_TICLclus_dR;
        
        
        // TDR photons
        double phoFromMultiClus_n;
        std::vector <double> v_phoFromMultiClus_E;
        std::vector <double> v_phoFromMultiClus_px;
        std::vector <double> v_phoFromMultiClus_py;
        std::vector <double> v_phoFromMultiClus_pz;
        std::vector <double> v_phoFromMultiClus_pT;
        std::vector <double> v_phoFromMultiClus_eta;
        std::vector <double> v_phoFromMultiClus_phi;
        
        std::vector <double> v_phoFromMultiClus_genPh_minDeltaR;
        std::vector <double> v_phoFromMultiClus_matchedGenPh_E;
        
        
        // TICL photons
        double phoFromTICL_n;
        std::vector <double> v_phoFromTICL_E;
        std::vector <double> v_phoFromTICL_px;
        std::vector <double> v_phoFromTICL_py;
        std::vector <double> v_phoFromTICL_pz;
        std::vector <double> v_phoFromTICL_pT;
        std::vector <double> v_phoFromTICL_eta;
        std::vector <double> v_phoFromTICL_phi;
        std::vector <double> v_phoFromTICL_ET;
        
        std::vector <double> v_phoFromTICL_genPh_minDeltaR;
        std::vector <double> v_phoFromTICL_nearestGenPh_idx;
        std::vector <double> v_phoFromTICL_matchedGenPh_E;
        std::vector <double> v_phoFromTICL_matchedGenPh_pT;
        std::vector <double> v_phoFromTICL_matchedGenPh_eta;
        std::vector <double> v_phoFromTICL_matchedGenPh_phi;
        
        std::vector <double> v_phoFromTICL_superClus_E;
        std::vector <double> v_phoFromTICL_superClus_ET;
        std::vector <double> v_phoFromTICL_superClus_rawE;
        std::vector <double> v_phoFromTICL_superClus_rawET;
        std::vector <double> v_phoFromTICL_superClus_theta;
        std::vector <double> v_phoFromTICL_superClus_eta;
        std::vector <double> v_phoFromTICL_superClus_phi;
        std::vector <double> v_phoFromTICL_superClus_x;
        std::vector <double> v_phoFromTICL_superClus_y;
        std::vector <double> v_phoFromTICL_superClus_z;
        std::vector <double> v_phoFromTICL_superClus_r;
        std::vector <double> v_phoFromTICL_superClus_etaWidth;
        std::vector <double> v_phoFromTICL_superClus_phiWidth;
        std::vector <double> v_phoFromTICL_superClus_nClus;
        std::vector <double> v_phoFromTICL_superClus_nHit;
        std::vector <double> v_phoFromTICL_superClus_clusMaxDR;
        
        std::vector <double> v_phoFromTICL_superClus_seed_dEta;
        std::vector <double> v_phoFromTICL_superClus_seed_dPhi;
        
        std::vector <double> v_phoFromTICL_superClus_recHit1_E;
        std::vector <double> v_phoFromTICL_superClus_recHit2_E;
        
        std::vector <double> v_phoFromTICL_superClus_sigma2rr;
        std::vector <double> v_phoFromTICL_superClus_sigma2etaEta;
        std::vector <double> v_phoFromTICL_superClus_sigma2phiPhi;
        
        std::vector <double> v_phoFromTICL_superClus_sigma2rEta;
        std::vector <double> v_phoFromTICL_superClus_sigma2rPhi;
        std::vector <double> v_phoFromTICL_superClus_sigma2etaPhi;
        
        std::vector <double> v_phoFromTICL_superClus_sigma2diag1;
        std::vector <double> v_phoFromTICL_superClus_sigma2diag2;
        std::vector <double> v_phoFromTICL_superClus_sigma2diag3;
        
        std::vector <double> v_phoFromTICL_superClusSeed_E;
        std::vector <double> v_phoFromTICL_superClusSeed_ET;
        std::vector <double> v_phoFromTICL_superClusSeed_eta;
        std::vector <double> v_phoFromTICL_superClusSeed_phi;
        
        std::vector <double> v_phoFromTICL_superClusSeed_recHit1_E;
        std::vector <double> v_phoFromTICL_superClusSeed_recHit2_E;
        
        std::vector <double> v_phoFromTICL_nRecHit;
        std::vector <std::vector <double> > vv_phoFromTICL_recHit_E;
        std::vector <std::vector <double> > vv_phoFromTICL_recHit_x;
        std::vector <std::vector <double> > vv_phoFromTICL_recHit_y;
        std::vector <std::vector <double> > vv_phoFromTICL_recHit_z;
        std::vector <std::vector <double> > vv_phoFromTICL_recHit_time;
        std::vector <std::vector <double> > vv_phoFromTICL_recHit_timeError;
        std::vector <std::vector <double> > vv_phoFromTICL_recHit_eta;
        std::vector <std::vector <double> > vv_phoFromTICL_recHit_phi;
        std::vector <std::vector <double> > vv_phoFromTICL_recHit_ET;
        std::vector <std::vector <double> > vv_phoFromTICL_recHit_detector;
        std::vector <std::vector <double> > vv_phoFromTICL_recHit_layer;
        std::vector <std::vector <double> > vv_phoFromTICL_recHit_isSimHitMatched;
        std::vector <std::vector <double> > vv_phoFromTICL_recHit_SCdEta;
        std::vector <std::vector <double> > vv_phoFromTICL_recHit_SCdPhi;
        std::vector <std::vector <double> > vv_phoFromTICL_recHit_SCdR;
        
        std::vector <double> v_phoFromTICL_nlcInCone;
        std::vector <std::vector <double> > vv_phoFromTICL_lcInCone_E;
        std::vector <std::vector <double> > vv_phoFromTICL_lcInCone_x;
        std::vector <std::vector <double> > vv_phoFromTICL_lcInCone_y;
        std::vector <std::vector <double> > vv_phoFromTICL_lcInCone_z;
        std::vector <std::vector <double> > vv_phoFromTICL_lcInCone_time;
        std::vector <std::vector <double> > vv_phoFromTICL_lcInCone_timeError;
        std::vector <std::vector <double> > vv_phoFromTICL_lcInCone_eta;
        std::vector <std::vector <double> > vv_phoFromTICL_lcInCone_phi;
        std::vector <std::vector <double> > vv_phoFromTICL_lcInCone_ET;
        std::vector <std::vector <double> > vv_phoFromTICL_lcInCone_detector;
        std::vector <std::vector <double> > vv_phoFromTICL_lcInCone_size;
        std::vector <std::vector <double> > vv_phoFromTICL_lcInCone_layer;
        std::vector <std::vector <double> > vv_phoFromTICL_lcInCone_SCdEta;
        std::vector <std::vector <double> > vv_phoFromTICL_lcInCone_SCdPhi;
        std::vector <std::vector <double> > vv_phoFromTICL_lcInCone_SCdR;
        
        
        
        struct RvarContent
        {
            std::vector <double> v_Rvar;
            
            
            void clear()
            {
                v_Rvar.clear();
            }
        };
        
        std::map <std::string, RvarContent*> m_RvarContent;
        
        
        struct PCAvarContent
        {
            std::vector <double> v_sigma2uu;
            std::vector <double> v_sigma2vv;
            std::vector <double> v_sigma2ww;
            
            
            void clear()
            {
                v_sigma2uu.clear();
                v_sigma2vv.clear();
                v_sigma2ww.clear();
            }
        };
        
        std::map <std::string, PCAvarContent*> m_PCAvarContent;
        
        
        struct IsoVarContent
        {
            std::vector <double> v_iso_sumETratio;
            std::vector <double> v_iso_sumETratio_charged;
            std::vector <double> v_iso_sumETratio_neutral;
            
            std::vector <double> v_iso_trackSumPt;
            
            std::vector <double> v_HoverE;
            
            void clear()
            {
                v_iso_sumETratio.clear();
                v_iso_sumETratio_charged.clear();
                v_iso_sumETratio_neutral.clear();
                
                v_iso_trackSumPt.clear();
                
                v_HoverE.clear();
            }
        };
        
        std::map <std::string, IsoVarContent*> m_isoVarContent;
        
        
        std::map <std::string, std::vector <double> > m_customVarContent;
        
        std::map <std::string, std::vector <double> > m_hitCount;
        
        
        double caloParticle_n;
        std::vector <double> v_caloParticle_E;
        std::vector <double> v_caloParticle_px;
        std::vector <double> v_caloParticle_py;
        std::vector <double> v_caloParticle_pz;
        std::vector <double> v_caloParticle_pT;
        std::vector <double> v_caloParticle_eta;
        std::vector <double> v_caloParticle_phi;
        std::vector <double> v_caloParticle_pdgid;
        
        char name[500];
        
        
        TreeOutput(std::string details, edm::Service<TFileService> fs)
        {
            //printf("Loading custom ROOT dictionaries. \n");
            //gROOT->ProcessLine(".L EDAnalyzers/TreeMaker/interface/CustomRootDict.cc+");
            //printf("Loaded custom ROOT dictionaries. \n");
            
            tree = fs->make<TTree>(details.c_str(), details.c_str());
            
            
            // Run info //
            tree->Branch("runNumber", &runNumber);
            tree->Branch("eventNumber", &eventNumber);
            tree->Branch("luminosityNumber", &luminosityNumber);
            tree->Branch("bunchCrossingNumber", &bunchCrossingNumber);
            
            
            // Gen event info //
            sprintf(name, "genEventWeight");
            tree->Branch(name, &genEventWeight);
            
            
            // Gen electron //
            sprintf(name, "genEl_n");
            tree->Branch(name, &genEl_n);
            
            sprintf(name, "genEl_E");
            tree->Branch(name, &v_genEl_E);
            
            sprintf(name, "genEl_px");
            tree->Branch(name, &v_genEl_px);
            
            sprintf(name, "genEl_py");
            tree->Branch(name, &v_genEl_py);
            
            sprintf(name, "genEl_pz");
            tree->Branch(name, &v_genEl_pz);
            
            sprintf(name, "genEl_pT");
            tree->Branch(name, &v_genEl_pT);
            
            sprintf(name, "genEl_eta");
            tree->Branch(name, &v_genEl_eta);
            
            sprintf(name, "genEl_phi");
            tree->Branch(name, &v_genEl_phi);
            
            sprintf(name, "genEl_HGCalEEP_EsortedIndex");
            tree->Branch(name, &v_genEl_HGCalEEP_EsortedIndex);
            
            sprintf(name, "genEl_HGCalEEM_EsortedIndex");
            tree->Branch(name, &v_genEl_HGCalEEM_EsortedIndex);
            
            sprintf(name, "genEl_multiClus_totE");
            tree->Branch(name, &v_genEl_multiClus_totE);
            
            sprintf(name, "genEl_multiClus_n");
            tree->Branch(name, &v_genEl_multiClus_n);
            
            sprintf(name, "genEl_nearestMultiClusEnRatio");
            tree->Branch(name, &v_genEl_nearestMultiClusEnRatio);
            
            sprintf(name, "genEl_multiClusEnRatio");
            tree->Branch(name, &v_genEl_multiClusEnRatio);
            
            
            // Gen photon //
            sprintf(name, "genPh_n");
            tree->Branch(name, &genPh_n);
            
            sprintf(name, "genPh_E");
            tree->Branch(name, &v_genPh_E);
            
            sprintf(name, "genPh_px");
            tree->Branch(name, &v_genPh_px);
            
            sprintf(name, "genPh_py");
            tree->Branch(name, &v_genPh_py);
            
            sprintf(name, "genPh_pz");
            tree->Branch(name, &v_genPh_pz);
            
            sprintf(name, "genPh_pT");
            tree->Branch(name, &v_genPh_pT);
            
            sprintf(name, "genPh_eta");
            tree->Branch(name, &v_genPh_eta);
            
            sprintf(name, "genPh_phi");
            tree->Branch(name, &v_genPh_phi);
            
            sprintf(name, "genPh_HGCalEEP_EsortedIndex");
            tree->Branch(name, &v_genPh_HGCalEEP_EsortedIndex);
            
            sprintf(name, "genPh_HGCalEEM_EsortedIndex");
            tree->Branch(name, &v_genPh_HGCalEEM_EsortedIndex);
            
            sprintf(name, "genPh_multiClus_totE");
            tree->Branch(name, &v_genPh_multiClus_totE);
            
            sprintf(name, "genPh_HGCalEEP_deltaR");
            tree->Branch(name, &v_genPh_HGCalEEP_deltaR);
            
            sprintf(name, "genPh_HGCalEEM_deltaR");
            tree->Branch(name, &v_genPh_HGCalEEM_deltaR);
            
            
            // HGCalEE+
            v_simHit_HGCalEEPlayer_totE.resize(Constants::HGCalEE_nLayer, 0.0);
            
            for(int iLayer = 0; iLayer < Constants::HGCalEE_nLayer; iLayer++)
            {
                sprintf(name, "simHit_HGCalEEPlayer%d_totE", iLayer+1);
                tree->Branch(name, &v_simHit_HGCalEEPlayer_totE.at(iLayer));
            }
            
            v_recHit_HGCalEEPlayer_totE.resize(Constants::HGCalEE_nLayer, 0.0);
            
            for(int iLayer = 0; iLayer < Constants::HGCalEE_nLayer; iLayer++)
            {
                sprintf(name, "recHit_HGCalEEPlayer%d_totE", iLayer+1);
                tree->Branch(name, &v_recHit_HGCalEEPlayer_totE.at(iLayer));
            }
            
            // HGCalEE-
            v_simHit_HGCalEEMlayer_totE.resize(Constants::HGCalEE_nLayer, 0.0);
            
            for(int iLayer = 0; iLayer < Constants::HGCalEE_nLayer; iLayer++)
            {
                sprintf(name, "simHit_HGCalEEMlayer%d_totE", iLayer+1);
                tree->Branch(name, &v_simHit_HGCalEEMlayer_totE.at(iLayer));
            }
            
            v_recHit_HGCalEEMlayer_totE.resize(Constants::HGCalEE_nLayer, 0.0);
            
            for(int iLayer = 0; iLayer < Constants::HGCalEE_nLayer; iLayer++)
            {
                sprintf(name, "recHit_HGCalEEMlayer%d_totE", iLayer+1);
                tree->Branch(name, &v_recHit_HGCalEEMlayer_totE.at(iLayer));
            }
            
            
            // Pileup //
            sprintf(name, "pileup_n");
            tree->Branch(name, &pileup_n);
            
            
            // Rho //
            sprintf(name, "rho");
            tree->Branch(name, &rho);
            
            
            sprintf(name, "nHit_EcalEB");
            tree->Branch(name, &nHit_EcalEB);
            
            sprintf(name, "nHit_HGCEE");
            tree->Branch(name, &nHit_HGCEE);
            
            sprintf(name, "nHit_HGCHEF");
            tree->Branch(name, &nHit_HGCHEF);
            
            sprintf(name, "nHit_HGCHEB");
            tree->Branch(name, &nHit_HGCHEB);
            
            
            // HGCAL layer clusters //
            sprintf(name, "HGCALlayerClus_n");
            tree->Branch(name, &HGCALlayerClus_n);
            
            sprintf(name, "HGCALlayerClus_E");
            tree->Branch(name, &v_HGCALlayerClus_E);
            
            sprintf(name, "HGCALlayerClus_x");
            tree->Branch(name, &v_HGCALlayerClus_x);
            
            sprintf(name, "HGCALlayerClus_y");
            tree->Branch(name, &v_HGCALlayerClus_y);
            
            sprintf(name, "HGCALlayerClus_z");
            tree->Branch(name, &v_HGCALlayerClus_z);
            
            sprintf(name, "HGCALlayerClus_eta");
            tree->Branch(name, &v_HGCALlayerClus_eta);
            
            sprintf(name, "HGCALlayerClus_phi");
            tree->Branch(name, &v_HGCALlayerClus_phi);
            
            sprintf(name, "HGCALlayerClus_ET");
            tree->Branch(name, &v_HGCALlayerClus_ET);
            
            sprintf(name, "HGCALlayerClus_time");
            tree->Branch(name, &v_HGCALlayerClus_time);
            
            sprintf(name, "HGCALlayerClus_timeError");
            tree->Branch(name, &v_HGCALlayerClus_timeError);
            
            sprintf(name, "HGCALlayerClus_detector");
            tree->Branch(name, &v_HGCALlayerClus_detector);
            
            sprintf(name, "HGCALlayerClus_layer");
            tree->Branch(name, &v_HGCALlayerClus_layer);
            
            
            // Tracksters //
            sprintf(name, "trackster_n");
            tree->Branch(name, &trackster_n);
            
            sprintf(name, "trackster_E");
            tree->Branch(name, &v_trackster_E);
            
            sprintf(name, "trackster_x");
            tree->Branch(name, &v_trackster_x);
            
            sprintf(name, "trackster_y");
            tree->Branch(name, &v_trackster_y);
            
            sprintf(name, "trackster_z");
            tree->Branch(name, &v_trackster_z);
            
            sprintf(name, "trackster_eta");
            tree->Branch(name, &v_trackster_eta);
            
            sprintf(name, "trackster_phi");
            tree->Branch(name, &v_trackster_phi);
            
            sprintf(name, "trackster_ET");
            tree->Branch(name, &v_trackster_ET);
            
            
            // MultiClusters //
            sprintf(name, "multiClus_n");
            tree->Branch(name, &multiClus_n);
            
            sprintf(name, "multiClus_genElIndex");
            tree->Branch(name, &v_multiClus_genElIndex);
            
            sprintf(name, "multiClus_E");
            tree->Branch(name, &v_multiClus_E);
            
            sprintf(name, "multiClus_x");
            tree->Branch(name, &v_multiClus_x);
            
            sprintf(name, "multiClus_y");
            tree->Branch(name, &v_multiClus_y);
            
            sprintf(name, "multiClus_z");
            tree->Branch(name, &v_multiClus_z);
            
            sprintf(name, "multiClus_eta");
            tree->Branch(name, &v_multiClus_eta);
            
            sprintf(name, "multiClus_phi");
            tree->Branch(name, &v_multiClus_phi);
            
            sprintf(name, "multiClus_ET");
            tree->Branch(name, &v_multiClus_ET);
            
            sprintf(name, "multiClus_corrE");
            tree->Branch(name, &v_multiClus_corrE);
            
            sprintf(name, "multiClus_corrET");
            tree->Branch(name, &v_multiClus_corrET);
            
            sprintf(name, "multiClus_dX");
            tree->Branch(name, &v_multiClus_dX);
            
            sprintf(name, "multiClus_dY");
            tree->Branch(name, &v_multiClus_dY);
            
            sprintf(name, "multiClus_dZ");
            tree->Branch(name, &v_multiClus_dZ);
            
            sprintf(name, "multiClus_dEta");
            tree->Branch(name, &v_multiClus_dEta);
            
            sprintf(name, "multiClus_dPhi");
            tree->Branch(name, &v_multiClus_dPhi);
            
            sprintf(name, "multiClus_sigma2rr");
            tree->Branch(name, &v_multiClus_sigma2rr);
            
            sprintf(name, "multiClus_sigma2etaEta");
            tree->Branch(name, &v_multiClus_sigma2etaEta);
            
            sprintf(name, "multiClus_sigma2phiPhi");
            tree->Branch(name, &v_multiClus_sigma2phiPhi);
            
            sprintf(name, "multiClus_sigma2rEta");
            tree->Branch(name, &v_multiClus_sigma2rEta);
            
            sprintf(name, "multiClus_sigma2rPhi");
            tree->Branch(name, &v_multiClus_sigma2rPhi);
            
            sprintf(name, "multiClus_sigma2etaPhi");
            tree->Branch(name, &v_multiClus_sigma2etaPhi);
            
            sprintf(name, "multiClus_sigma2diag1");
            tree->Branch(name, &v_multiClus_sigma2diag1);
            
            sprintf(name, "multiClus_sigma2diag2");
            tree->Branch(name, &v_multiClus_sigma2diag2);
            
            sprintf(name, "multiClus_sigma2diag3");
            tree->Branch(name, &v_multiClus_sigma2diag3);
            
            sprintf(name, "multiClus_EsortedIndex");
            tree->Branch(name, &v_multiClus_EsortedIndex);
            
            sprintf(name, "multiClus_HGCalEEP_EsortedIndex");
            tree->Branch(name, &v_multiClus_HGCalEEP_EsortedIndex);
            
            sprintf(name, "multiClus_HGCalEEM_EsortedIndex");
            tree->Branch(name, &v_multiClus_HGCalEEM_EsortedIndex);
            
            
            //
            sprintf(name, "multiClus_clus_n");
            tree->Branch(name, &v_multiClus_clus_n);
            
            sprintf(name, "multiClus_clus_startIndex");
            tree->Branch(name, &v_multiClus_clus_startIndex);
            
            sprintf(name, "multiClus_clus_E");
            tree->Branch(name, &v_multiClus_clus_E);
            
            sprintf(name, "multiClus_clus_x");
            tree->Branch(name, &v_multiClus_clus_x);
            
            sprintf(name, "multiClus_clus_y");
            tree->Branch(name, &v_multiClus_clus_y);
            
            sprintf(name, "multiClus_clus_z");
            tree->Branch(name, &v_multiClus_clus_z);
            
            sprintf(name, "multiClus_clus_eta");
            tree->Branch(name, &v_multiClus_clus_eta);
            
            sprintf(name, "multiClus_clus_phi");
            tree->Branch(name, &v_multiClus_clus_phi);
            
            sprintf(name, "multiClus_clus_ET");
            tree->Branch(name, &v_multiClus_clus_ET);
            
            
            sprintf(name, "multiClus_clus_layer");
            tree->Branch(name, &v_multiClus_clus_layer);
            
            sprintf(name, "multiClus_clus_multiplicity");
            tree->Branch(name, &v_multiClus_clus_multiplicity);
            
            sprintf(name, "multiClus_clus_nHit");
            tree->Branch(name, &v_multiClus_clus_nHit);
            
            
            //
            sprintf(name, "multiClus_uniqueClus_n");
            tree->Branch(name, &v_multiClus_uniqueClus_n);
            
            sprintf(name, "multiClus_uniqueClus_E");
            tree->Branch(name, &v_multiClus_uniqueClus_E);
            
            sprintf(name, "multiClus_uniqueClus_x");
            tree->Branch(name, &v_multiClus_uniqueClus_x);
            
            sprintf(name, "multiClus_uniqueClus_y");
            tree->Branch(name, &v_multiClus_uniqueClus_y);
            
            sprintf(name, "multiClus_uniqueClus_z");
            tree->Branch(name, &v_multiClus_uniqueClus_z);
            
            sprintf(name, "multiClus_uniqueClus_eta");
            tree->Branch(name, &v_multiClus_uniqueClus_eta);
            
            sprintf(name, "multiClus_uniqueClus_phi");
            tree->Branch(name, &v_multiClus_uniqueClus_phi);
            
            sprintf(name, "multiClus_uniqueClus_ET");
            tree->Branch(name, &v_multiClus_uniqueClus_ET);
            
            
            sprintf(name, "multiClus_uniqueClus_layer");
            tree->Branch(name, &v_multiClus_uniqueClus_layer);
            
            sprintf(name, "multiClus_uniqueClus_multiplicity");
            tree->Branch(name, &v_multiClus_uniqueClus_multiplicity);
            
            sprintf(name, "multiClus_uniqueClus_nHit");
            tree->Branch(name, &v_multiClus_uniqueClus_nHit);
            
            
            //
            sprintf(name, "multiClus_mc1_dX");
            tree->Branch(name, &v_multiClus_mc1_dX);
            
            sprintf(name, "multiClus_mc1_dY");
            tree->Branch(name, &v_multiClus_mc1_dY);
            
            sprintf(name, "multiClus_mc1_dZ");
            tree->Branch(name, &v_multiClus_mc1_dZ);
            
            sprintf(name, "multiClus_mc1_dEta");
            tree->Branch(name, &v_multiClus_mc1_dEta);
            
            sprintf(name, "multiClus_mc1_dPhi");
            tree->Branch(name, &v_multiClus_mc1_dPhi);
            
            sprintf(name, "multiClus_mc1_dR");
            tree->Branch(name, &v_multiClus_mc1_dR);
            
            
            //
            sprintf(name, "multiClus_HGCalEEP_meanX");
            tree->Branch(name, &multiClus_HGCalEEP_meanX);
            
            sprintf(name, "multiClus_HGCalEEP_meanY");
            tree->Branch(name, &multiClus_HGCalEEP_meanY);
            
            sprintf(name, "multiClus_HGCalEEP_meanZ");
            tree->Branch(name, &multiClus_HGCalEEP_meanZ);
            
            sprintf(name, "multiClus_HGCalEEP_meanDx");
            tree->Branch(name, &multiClus_HGCalEEP_meanDx);
            
            sprintf(name, "multiClus_HGCalEEP_meanDy");
            tree->Branch(name, &multiClus_HGCalEEP_meanDy);
            
            sprintf(name, "multiClus_HGCalEEP_meanDz");
            tree->Branch(name, &multiClus_HGCalEEP_meanDz);
            
            sprintf(name, "multiClus_HGCalEEP_totE");
            tree->Branch(name, &multiClus_HGCalEEP_totE);
            
            sprintf(name, "multiClus_HGCalEEP_totET");
            tree->Branch(name, &multiClus_HGCalEEP_totET);
            
            sprintf(name, "multiClus_HGCalEEP_diag1");
            tree->Branch(name, &multiClus_HGCalEEP_diag1);
            
            sprintf(name, "multiClus_HGCalEEP_diag2");
            tree->Branch(name, &multiClus_HGCalEEP_diag2);
            
            sprintf(name, "multiClus_HGCalEEP_diag3");
            tree->Branch(name, &multiClus_HGCalEEP_diag3);
            
            //
            sprintf(name, "multiClus_HGCalEEM_meanX");
            tree->Branch(name, &multiClus_HGCalEEM_meanX);
            
            sprintf(name, "multiClus_HGCalEEM_meanY");
            tree->Branch(name, &multiClus_HGCalEEM_meanY);
            
            sprintf(name, "multiClus_HGCalEEM_meanZ");
            tree->Branch(name, &multiClus_HGCalEEM_meanZ);
            
            sprintf(name, "multiClus_HGCalEEM_meanDx");
            tree->Branch(name, &multiClus_HGCalEEM_meanDx);
            
            sprintf(name, "multiClus_HGCalEEM_meanDy");
            tree->Branch(name, &multiClus_HGCalEEM_meanDy);
            
            sprintf(name, "multiClus_HGCalEEM_meanDz");
            tree->Branch(name, &multiClus_HGCalEEM_meanDz);
            
            sprintf(name, "multiClus_HGCalEEM_totE");
            tree->Branch(name, &multiClus_HGCalEEM_totE);
            
            sprintf(name, "multiClus_HGCalEEM_totET");
            tree->Branch(name, &multiClus_HGCalEEM_totET);
            
            sprintf(name, "multiClus_HGCalEEM_diag1");
            tree->Branch(name, &multiClus_HGCalEEM_diag1);
            
            sprintf(name, "multiClus_HGCalEEM_diag2");
            tree->Branch(name, &multiClus_HGCalEEM_diag2);
            
            sprintf(name, "multiClus_HGCalEEM_diag3");
            tree->Branch(name, &multiClus_HGCalEEM_diag3);
            
            
            //
            sprintf(name, "simHit_n");
            tree->Branch(name, &simHit_n);
            
            sprintf(name, "simHit_E");
            tree->Branch(name, &v_simHit_E);
            
            sprintf(name, "simHit_x");
            tree->Branch(name, &v_simHit_x);
            
            sprintf(name, "simHit_y");
            tree->Branch(name, &v_simHit_y);
            
            sprintf(name, "simHit_z");
            tree->Branch(name, &v_simHit_z);
            
            sprintf(name, "simHit_eta");
            tree->Branch(name, &v_simHit_eta);
            
            sprintf(name, "simHit_phi");
            tree->Branch(name, &v_simHit_phi);
            
            sprintf(name, "simHit_ET");
            tree->Branch(name, &v_simHit_ET);
            
            sprintf(name, "simHit_layer");
            tree->Branch(name, &v_simHit_layer);
            
            sprintf(name, "simHit_zside");
            tree->Branch(name, &v_simHit_zside);
            
            sprintf(name, "simHit_isCaloParticleMatched");
            tree->Branch(name, &v_simHit_isCaloParticleMatched);
            
            
            //
            sprintf(name, "recHit_n");
            tree->Branch(name, &recHit_n);
            
            sprintf(name, "recHit_E");
            tree->Branch(name, &v_recHit_E);
            
            sprintf(name, "recHit_x");
            tree->Branch(name, &v_recHit_x);
            
            sprintf(name, "recHit_y");
            tree->Branch(name, &v_recHit_y);
            
            sprintf(name, "recHit_z");
            tree->Branch(name, &v_recHit_z);
            
            sprintf(name, "recHit_eta");
            tree->Branch(name, &v_recHit_eta);
            
            sprintf(name, "recHit_phi");
            tree->Branch(name, &v_recHit_phi);
            
            sprintf(name, "recHit_ET");
            tree->Branch(name, &v_recHit_ET);
            
            sprintf(name, "recHit_layer");
            tree->Branch(name, &v_recHit_layer);
            
            sprintf(name, "recHit_zside");
            tree->Branch(name, &v_recHit_zside);
            
            sprintf(name, "recHit_matchedSimHitIndex");
            tree->Branch(name, &v_recHit_matchedSimHitIndex);
            
            sprintf(name, "recHit_matchedHGCALlayerClusIndex");
            tree->Branch(name, &v_recHit_matchedHGCALlayerClusIndex);
            
            sprintf(name, "recHit_isMultiClusMatched");
            tree->Branch(name, &v_recHit_isMultiClusMatched);
            
            sprintf(name, "recHit_isCaloParticleMatched");
            tree->Branch(name, &v_recHit_isCaloParticleMatched);
            
            sprintf(name, "recHit_iType");
            tree->Branch(name, &v_recHit_iType);
            
            sprintf(name, "recHit_iCell1");
            tree->Branch(name, &v_recHit_iCell1);
            
            sprintf(name, "recHit_iCell2");
            tree->Branch(name, &v_recHit_iCell2);
            
            sprintf(name, "recHit_SiThickness");
            tree->Branch(name, &v_recHit_SiThickness);
            
            
            // TDR electrons
            sprintf(name, "gsfEleFromMultiClus_n");
            tree->Branch(name, &gsfEleFromMultiClus_n);
            
            sprintf(name, "gsfEleFromMultiClus_E");
            tree->Branch(name, &v_gsfEleFromMultiClus_E);
            
            sprintf(name, "gsfEleFromMultiClus_px");
            tree->Branch(name, &v_gsfEleFromMultiClus_px);
            
            sprintf(name, "gsfEleFromMultiClus_py");
            tree->Branch(name, &v_gsfEleFromMultiClus_py);
            
            sprintf(name, "gsfEleFromMultiClus_pz");
            tree->Branch(name, &v_gsfEleFromMultiClus_pz);
            
            sprintf(name, "gsfEleFromMultiClus_pT");
            tree->Branch(name, &v_gsfEleFromMultiClus_pT);
            
            sprintf(name, "gsfEleFromMultiClus_eta");
            tree->Branch(name, &v_gsfEleFromMultiClus_eta);
            
            sprintf(name, "gsfEleFromMultiClus_phi");
            tree->Branch(name, &v_gsfEleFromMultiClus_phi);
            
            sprintf(name, "gsfEleFromMultiClus_genEl_minDeltaR");
            tree->Branch(name, &v_gsfEleFromMultiClus_genEl_minDeltaR);
            
            sprintf(name, "gsfEleFromMultiClus_matchedGenEl_E");
            tree->Branch(name, &v_gsfEleFromMultiClus_matchedGenEl_E);
            
            
            // TICL electrons
            sprintf(name, "gsfEleFromTICL_n");
            tree->Branch(name, &gsfEleFromTICL_n);
            
            sprintf(name, "gsfEleFromTICL_E");
            tree->Branch(name, &v_gsfEleFromTICL_E);
            
            sprintf(name, "gsfEleFromTICL_px");
            tree->Branch(name, &v_gsfEleFromTICL_px);
            
            sprintf(name, "gsfEleFromTICL_py");
            tree->Branch(name, &v_gsfEleFromTICL_py);
            
            sprintf(name, "gsfEleFromTICL_pz");
            tree->Branch(name, &v_gsfEleFromTICL_pz);
            
            sprintf(name, "gsfEleFromTICL_pT");
            tree->Branch(name, &v_gsfEleFromTICL_pT);
            
            sprintf(name, "gsfEleFromTICL_eta");
            tree->Branch(name, &v_gsfEleFromTICL_eta);
            
            sprintf(name, "gsfEleFromTICL_phi");
            tree->Branch(name, &v_gsfEleFromTICL_phi);
            
            sprintf(name, "gsfEleFromTICL_ET");
            tree->Branch(name, &v_gsfEleFromTICL_ET);
            
            sprintf(name, "gsfEleFromTICL_genEl_minDeltaR");
            tree->Branch(name, &v_gsfEleFromTICL_genEl_minDeltaR);
            
            sprintf(name, "gsfEleFromTICL_nearestGenEl_idx");
            tree->Branch(name, &v_gsfEleFromTICL_nearestGenEl_idx);
            
            sprintf(name, "gsfEleFromTICL_matchedGenEl_E");
            tree->Branch(name, &v_gsfEleFromTICL_matchedGenEl_E);
            
            sprintf(name, "gsfEleFromTICL_matchedGenEl_pT");
            tree->Branch(name, &v_gsfEleFromTICL_matchedGenEl_pT);
            
            sprintf(name, "gsfEleFromTICL_matchedGenEl_eta");
            tree->Branch(name, &v_gsfEleFromTICL_matchedGenEl_eta);
            
            sprintf(name, "gsfEleFromTICL_matchedGenEl_phi");
            tree->Branch(name, &v_gsfEleFromTICL_matchedGenEl_phi);
            
            
            sprintf(name, "gsfEleFromTICL_gsfTrack_p");
            tree->Branch(name, &v_gsfEleFromTICL_gsfTrack_p);
            
            sprintf(name, "gsfEleFromTICL_gsfTrack_px");
            tree->Branch(name, &v_gsfEleFromTICL_gsfTrack_px);
            
            sprintf(name, "gsfEleFromTICL_gsfTrack_py");
            tree->Branch(name, &v_gsfEleFromTICL_gsfTrack_py);
            
            sprintf(name, "gsfEleFromTICL_gsfTrack_pz");
            tree->Branch(name, &v_gsfEleFromTICL_gsfTrack_pz);
            
            sprintf(name, "gsfEleFromTICL_gsfTrack_pT");
            tree->Branch(name, &v_gsfEleFromTICL_gsfTrack_pT);
            
            sprintf(name, "gsfEleFromTICL_gsfTrack_eta");
            tree->Branch(name, &v_gsfEleFromTICL_gsfTrack_eta);
            
            sprintf(name, "gsfEleFromTICL_gsfTrack_phi");
            tree->Branch(name, &v_gsfEleFromTICL_gsfTrack_phi);
            
            
            sprintf(name, "gsfEleFromTICL_trkAtVtx_p");
            tree->Branch(name, &v_gsfEleFromTICL_trkAtVtx_p);
            
            sprintf(name, "gsfEleFromTICL_trkAtVtx_px");
            tree->Branch(name, &v_gsfEleFromTICL_trkAtVtx_px);
            
            sprintf(name, "gsfEleFromTICL_trkAtVtx_py");
            tree->Branch(name, &v_gsfEleFromTICL_trkAtVtx_py);
            
            sprintf(name, "gsfEleFromTICL_trkAtVtx_pz");
            tree->Branch(name, &v_gsfEleFromTICL_trkAtVtx_pz);
            
            sprintf(name, "gsfEleFromTICL_trkAtVtx_pT");
            tree->Branch(name, &v_gsfEleFromTICL_trkAtVtx_pT);
            
            sprintf(name, "gsfEleFromTICL_trkAtVtx_eta");
            tree->Branch(name, &v_gsfEleFromTICL_trkAtVtx_eta);
            
            sprintf(name, "gsfEleFromTICL_trkAtVtx_phi");
            tree->Branch(name, &v_gsfEleFromTICL_trkAtVtx_phi);
            
            
            sprintf(name, "gsfEleFromTICL_dr03TkSumPt");
            tree->Branch(name, &v_gsfEleFromTICL_dr03TkSumPt);
            
            
            sprintf(name, "gsfEleFromTICL_superClus_E");
            tree->Branch(name, &v_gsfEleFromTICL_superClus_E);
            
            sprintf(name, "gsfEleFromTICL_superClus_ET");
            tree->Branch(name, &v_gsfEleFromTICL_superClus_ET);
            
            sprintf(name, "gsfEleFromTICL_superClus_rawE");
            tree->Branch(name, &v_gsfEleFromTICL_superClus_rawE);
            
            sprintf(name, "gsfEleFromTICL_superClus_rawET");
            tree->Branch(name, &v_gsfEleFromTICL_superClus_rawET);
            
            sprintf(name, "gsfEleFromTICL_superClus_theta");
            tree->Branch(name, &v_gsfEleFromTICL_superClus_theta);
            
            sprintf(name, "gsfEleFromTICL_superClus_eta");
            tree->Branch(name, &v_gsfEleFromTICL_superClus_eta);
            
            sprintf(name, "gsfEleFromTICL_superClus_phi");
            tree->Branch(name, &v_gsfEleFromTICL_superClus_phi);
            
            sprintf(name, "gsfEleFromTICL_superClus_x");
            tree->Branch(name, &v_gsfEleFromTICL_superClus_x);
            
            sprintf(name, "gsfEleFromTICL_superClus_y");
            tree->Branch(name, &v_gsfEleFromTICL_superClus_y);
            
            sprintf(name, "gsfEleFromTICL_superClus_z");
            tree->Branch(name, &v_gsfEleFromTICL_superClus_z);
            
            sprintf(name, "gsfEleFromTICL_superClus_r");
            tree->Branch(name, &v_gsfEleFromTICL_superClus_r);
            
            sprintf(name, "gsfEleFromTICL_superClus_etaWidth");
            tree->Branch(name, &v_gsfEleFromTICL_superClus_etaWidth);
            
            sprintf(name, "gsfEleFromTICL_superClus_phiWidth");
            tree->Branch(name, &v_gsfEleFromTICL_superClus_phiWidth);
            
            sprintf(name, "gsfEleFromTICL_superClus_nClus");
            tree->Branch(name, &v_gsfEleFromTICL_superClus_nClus);
            
            sprintf(name, "gsfEleFromTICL_superClus_nHit");
            tree->Branch(name, &v_gsfEleFromTICL_superClus_nHit);
            
            sprintf(name, "gsfEleFromTICL_superClus_nearestCellDist");
            tree->Branch(name, &v_gsfEleFromTICL_superClus_nearestCellDist);
            
            sprintf(name, "gsfEleFromTICL_superClus_cellNeighbour1ringWindow_n");
            tree->Branch(name, &v_gsfEleFromTICL_superClus_cellNeighbour1ringWindow_n);
            
            sprintf(name, "gsfEleFromTICL_superClus_cellNeighbour2ringWindow_n");
            tree->Branch(name, &v_gsfEleFromTICL_superClus_cellNeighbour2ringWindow_n);
            
            sprintf(name, "gsfEleFromTICL_superClus_clusMaxDR");
            tree->Branch(name, &v_gsfEleFromTICL_superClus_clusMaxDR);
            
            
            sprintf(name, "gsfEleFromTICL_superClus_seed_dEta");
            tree->Branch(name, &v_gsfEleFromTICL_superClus_seed_dEta);
            
            sprintf(name, "gsfEleFromTICL_superClus_seed_dPhi");
            tree->Branch(name, &v_gsfEleFromTICL_superClus_seed_dPhi);
            
            
            sprintf(name, "gsfEleFromTICL_superClus_recHit1_E");
            tree->Branch(name, &v_gsfEleFromTICL_superClus_recHit1_E);
            
            sprintf(name, "gsfEleFromTICL_superClus_recHit2_E");
            tree->Branch(name, &v_gsfEleFromTICL_superClus_recHit2_E);
            
            
            sprintf(name, "gsfEleFromTICL_superClus_sigma2etaEtaLW");
            tree->Branch(name, &v_gsfEleFromTICL_superClus_sigma2etaEtaLW);
            
            sprintf(name, "gsfEleFromTICL_superClus_sigma2phiPhiLW");
            tree->Branch(name, &v_gsfEleFromTICL_superClus_sigma2phiPhiLW);
            
            
            sprintf(name, "gsfEleFromTICL_superClus_sigma2rr");
            tree->Branch(name, &v_gsfEleFromTICL_superClus_sigma2rr);
            
            sprintf(name, "gsfEleFromTICL_superClus_sigma2etaEta");
            tree->Branch(name, &v_gsfEleFromTICL_superClus_sigma2etaEta);
            
            sprintf(name, "gsfEleFromTICL_superClus_sigma2phiPhi");
            tree->Branch(name, &v_gsfEleFromTICL_superClus_sigma2phiPhi);
            
            sprintf(name, "gsfEleFromTICL_superClus_sigma2rEta");
            tree->Branch(name, &v_gsfEleFromTICL_superClus_sigma2rEta);
            
            sprintf(name, "gsfEleFromTICL_superClus_sigma2rPhi");
            tree->Branch(name, &v_gsfEleFromTICL_superClus_sigma2rPhi);
            
            sprintf(name, "gsfEleFromTICL_superClus_sigma2etaPhi");
            tree->Branch(name, &v_gsfEleFromTICL_superClus_sigma2etaPhi);
            
            sprintf(name, "gsfEleFromTICL_superClus_sigma2diag1");
            tree->Branch(name, &v_gsfEleFromTICL_superClus_sigma2diag1);
            
            sprintf(name, "gsfEleFromTICL_superClus_sigma2diag2");
            tree->Branch(name, &v_gsfEleFromTICL_superClus_sigma2diag2);
            
            sprintf(name, "gsfEleFromTICL_superClus_sigma2diag3");
            tree->Branch(name, &v_gsfEleFromTICL_superClus_sigma2diag3);
            
            
            sprintf(name, "gsfEleFromTICL_E7");
            tree->Branch(name, &v_gsfEleFromTICL_E7);
            
            sprintf(name, "gsfEleFromTICL_R7");
            tree->Branch(name, &v_gsfEleFromTICL_R7);
            
            sprintf(name, "gsfEleFromTICL_E19");
            tree->Branch(name, &v_gsfEleFromTICL_E19);
            
            sprintf(name, "gsfEleFromTICL_R19");
            tree->Branch(name, &v_gsfEleFromTICL_R19);
            
            //
            sprintf(name, "gsfEleFromTICL_R2p0");
            tree->Branch(name, &v_gsfEleFromTICL_R2p0);
            
            //
            sprintf(name, "gsfEleFromTICL_E2p4");
            tree->Branch(name, &v_gsfEleFromTICL_E2p4);
            
            sprintf(name, "gsfEleFromTICL_R2p4");
            tree->Branch(name, &v_gsfEleFromTICL_R2p4);
            
            //
            sprintf(name, "gsfEleFromTICL_E2p6");
            tree->Branch(name, &v_gsfEleFromTICL_E2p6);
            
            sprintf(name, "gsfEleFromTICL_R2p6");
            tree->Branch(name, &v_gsfEleFromTICL_R2p6);
            
            //
            sprintf(name, "gsfEleFromTICL_E2p8");
            tree->Branch(name, &v_gsfEleFromTICL_E2p8);
            
            sprintf(name, "gsfEleFromTICL_R2p8");
            tree->Branch(name, &v_gsfEleFromTICL_R2p8);
            
            //
            sprintf(name, "gsfEleFromTICL_E3p0");
            tree->Branch(name, &v_gsfEleFromTICL_E3p0);
            
            sprintf(name, "gsfEleFromTICL_R3p0");
            tree->Branch(name, &v_gsfEleFromTICL_R3p0);
            
            //
            sprintf(name, "gsfEleFromTICL_E3p5");
            tree->Branch(name, &v_gsfEleFromTICL_E3p5);
            
            sprintf(name, "gsfEleFromTICL_R3p5");
            tree->Branch(name, &v_gsfEleFromTICL_R3p5);
            
            //
            sprintf(name, "gsfEleFromTICL_superClusSeed_E");
            tree->Branch(name, &v_gsfEleFromTICL_superClusSeed_E);
            
            sprintf(name, "gsfEleFromTICL_superClusSeed_ET");
            tree->Branch(name, &v_gsfEleFromTICL_superClusSeed_ET);
            
            sprintf(name, "gsfEleFromTICL_superClusSeed_eta");
            tree->Branch(name, &v_gsfEleFromTICL_superClusSeed_eta);
            
            sprintf(name, "gsfEleFromTICL_superClusSeed_phi");
            tree->Branch(name, &v_gsfEleFromTICL_superClusSeed_phi);
            
            sprintf(name, "gsfEleFromTICL_superClusSeed_recHit1_E");
            tree->Branch(name, &v_gsfEleFromTICL_superClusSeed_recHit1_E);
            
            sprintf(name, "gsfEleFromTICL_superClusSeed_recHit2_E");
            tree->Branch(name, &v_gsfEleFromTICL_superClusSeed_recHit2_E);
            
            
            sprintf(name, "gsfEleFromTICL_superClusSeed_clus_dEta");
            tree->Branch(name, &vv_gsfEleFromTICL_superClusSeed_clus_dEta);
            
            sprintf(name, "gsfEleFromTICL_superClusSeed_clus_dPhi");
            tree->Branch(name, &vv_gsfEleFromTICL_superClusSeed_clus_dPhi);
            
            sprintf(name, "gsfEleFromTICL_superClus_TICLclus_n");
            tree->Branch(name, &v_gsfEleFromTICL_superClus_TICLclus_n);
            
            sprintf(name, "gsfEleFromTICL_superClus_TICLclus_E");
            tree->Branch(name, &vv_gsfEleFromTICL_superClus_TICLclus_E);
            
            sprintf(name, "gsfEleFromTICL_superClus_TICLclus_ET");
            tree->Branch(name, &vv_gsfEleFromTICL_superClus_TICLclus_ET);
            
            sprintf(name, "gsfEleFromTICL_superClus_TICLclus_nClus");
            tree->Branch(name, &vv_gsfEleFromTICL_superClus_TICLclus_nClus);
            
            sprintf(name, "gsfEleFromTICL_superClusSeed_TICLclus_dX");
            tree->Branch(name, &vv_gsfEleFromTICL_superClusSeed_TICLclus_dX);
            
            sprintf(name, "gsfEleFromTICL_superClusSeed_TICLclus_dY");
            tree->Branch(name, &vv_gsfEleFromTICL_superClusSeed_TICLclus_dY);
            
            sprintf(name, "gsfEleFromTICL_superClusSeed_TICLclus_dZ");
            tree->Branch(name, &vv_gsfEleFromTICL_superClusSeed_TICLclus_dZ);
            
            sprintf(name, "gsfEleFromTICL_superClusSeed_TICLclus_dEta");
            tree->Branch(name, &vv_gsfEleFromTICL_superClusSeed_TICLclus_dEta);
            
            sprintf(name, "gsfEleFromTICL_superClusSeed_TICLclus_dPhi");
            tree->Branch(name, &vv_gsfEleFromTICL_superClusSeed_TICLclus_dPhi);
            
            sprintf(name, "gsfEleFromTICL_superClusSeed_TICLclus_dR");
            tree->Branch(name, &vv_gsfEleFromTICL_superClusSeed_TICLclus_dR);
            
            
            vv_gsfEleFromTICL_superClus_E_layer.resize(Constants::HGCalEE_nLayer, {});
            
            vv_gsfEleFromTICL_E2p4_layer.resize(Constants::HGCalEE_nLayer, {});
            vv_gsfEleFromTICL_R2p4_layer.resize(Constants::HGCalEE_nLayer, {});
            
            vv_gsfEleFromTICL_E2p6_layer.resize(Constants::HGCalEE_nLayer, {});
            vv_gsfEleFromTICL_R2p6_layer.resize(Constants::HGCalEE_nLayer, {});
            
            vv_gsfEleFromTICL_E2p8_layer.resize(Constants::HGCalEE_nLayer, {});
            vv_gsfEleFromTICL_R2p8_layer.resize(Constants::HGCalEE_nLayer, {});
            
            vv_gsfEleFromTICL_E3p0_layer.resize(Constants::HGCalEE_nLayer, {});
            vv_gsfEleFromTICL_R3p0_layer.resize(Constants::HGCalEE_nLayer, {});
            
            vv_gsfEleFromTICL_E3p5_layer.resize(Constants::HGCalEE_nLayer, {});
            vv_gsfEleFromTICL_R3p5_layer.resize(Constants::HGCalEE_nLayer, {});
            
            for(int iLayer = 0; iLayer < Constants::HGCalEE_nLayer; iLayer++)
            {
                sprintf(name, "gsfEleFromTICL_superClus_E_layer%d", iLayer+1);
                tree->Branch(name, &vv_gsfEleFromTICL_superClus_E_layer.at(iLayer));
                
                
                //
                //sprintf(name, "gsfEleFromTICL_E2p4_layer%d", iLayer+1);
                //tree->Branch(name, &vv_gsfEleFromTICL_E2p4_layer.at(iLayer));
                //
                //sprintf(name, "gsfEleFromTICL_R2p4_layer%d", iLayer+1);
                //tree->Branch(name, &vv_gsfEleFromTICL_R2p4_layer.at(iLayer));
                //
                ////
                //sprintf(name, "gsfEleFromTICL_E2p6_layer%d", iLayer+1);
                //tree->Branch(name, &vv_gsfEleFromTICL_E2p6_layer.at(iLayer));
                //
                //sprintf(name, "gsfEleFromTICL_R2p6_layer%d", iLayer+1);
                //tree->Branch(name, &vv_gsfEleFromTICL_R2p6_layer.at(iLayer));
                
                //
                sprintf(name, "gsfEleFromTICL_E2p8_layer%d", iLayer+1);
                tree->Branch(name, &vv_gsfEleFromTICL_E2p8_layer.at(iLayer));
                
                sprintf(name, "gsfEleFromTICL_R2p8_layer%d", iLayer+1);
                tree->Branch(name, &vv_gsfEleFromTICL_R2p8_layer.at(iLayer));
                
                //
                //sprintf(name, "gsfEleFromTICL_E3p0_layer%d", iLayer+1);
                //tree->Branch(name, &vv_gsfEleFromTICL_E3p0_layer.at(iLayer));
                //
                //sprintf(name, "gsfEleFromTICL_R3p0_layer%d", iLayer+1);
                //tree->Branch(name, &vv_gsfEleFromTICL_R3p0_layer.at(iLayer));
                //
                ////
                //sprintf(name, "gsfEleFromTICL_E3p5_layer%d", iLayer+1);
                //tree->Branch(name, &vv_gsfEleFromTICL_E3p5_layer.at(iLayer));
                //
                //sprintf(name, "gsfEleFromTICL_R3p5_layer%d", iLayer+1);
                //tree->Branch(name, &vv_gsfEleFromTICL_R3p5_layer.at(iLayer));
            }
            
            
            // TDR photons
            sprintf(name, "phoFromMultiClus_n");
            tree->Branch(name, &phoFromMultiClus_n);
            
            sprintf(name, "phoFromMultiClus_E");
            tree->Branch(name, &v_phoFromMultiClus_E);
            
            sprintf(name, "phoFromMultiClus_px");
            tree->Branch(name, &v_phoFromMultiClus_px);
            
            sprintf(name, "phoFromMultiClus_py");
            tree->Branch(name, &v_phoFromMultiClus_py);
            
            sprintf(name, "phoFromMultiClus_pz");
            tree->Branch(name, &v_phoFromMultiClus_pz);
            
            sprintf(name, "phoFromMultiClus_pT");
            tree->Branch(name, &v_phoFromMultiClus_pT);
            
            sprintf(name, "phoFromMultiClus_eta");
            tree->Branch(name, &v_phoFromMultiClus_eta);
            
            sprintf(name, "phoFromMultiClus_phi");
            tree->Branch(name, &v_phoFromMultiClus_phi);
            
            sprintf(name, "phoFromMultiClus_genPh_minDeltaR");
            tree->Branch(name, &v_phoFromMultiClus_genPh_minDeltaR);
            
            sprintf(name, "phoFromMultiClus_matchedGenPh_E");
            tree->Branch(name, &v_phoFromMultiClus_matchedGenPh_E);
            
            
            // TICL photons
            sprintf(name, "phoFromTICL_n");
            tree->Branch(name, &phoFromTICL_n);
            
            sprintf(name, "phoFromTICL_E");
            tree->Branch(name, &v_phoFromTICL_E);
            
            sprintf(name, "phoFromTICL_px");
            tree->Branch(name, &v_phoFromTICL_px);
            
            sprintf(name, "phoFromTICL_py");
            tree->Branch(name, &v_phoFromTICL_py);
            
            sprintf(name, "phoFromTICL_pz");
            tree->Branch(name, &v_phoFromTICL_pz);
            
            sprintf(name, "phoFromTICL_pT");
            tree->Branch(name, &v_phoFromTICL_pT);
            
            sprintf(name, "phoFromTICL_eta");
            tree->Branch(name, &v_phoFromTICL_eta);
            
            sprintf(name, "phoFromTICL_phi");
            tree->Branch(name, &v_phoFromTICL_phi);
            
            sprintf(name, "phoFromTICL_ET");
            tree->Branch(name, &v_phoFromTICL_ET);
            
            sprintf(name, "phoFromTICL_genPh_minDeltaR");
            tree->Branch(name, &v_phoFromTICL_genPh_minDeltaR);
            
            sprintf(name, "phoFromTICL_nearestGenPh_idx");
            tree->Branch(name, &v_phoFromTICL_nearestGenPh_idx);
            
            sprintf(name, "phoFromTICL_matchedGenPh_E");
            tree->Branch(name, &v_phoFromTICL_matchedGenPh_E);
            
            sprintf(name, "phoFromTICL_matchedGenPh_pT");
            tree->Branch(name, &v_phoFromTICL_matchedGenPh_pT);
            
            sprintf(name, "phoFromTICL_matchedGenPh_eta");
            tree->Branch(name, &v_phoFromTICL_matchedGenPh_eta);
            
            sprintf(name, "phoFromTICL_matchedGenPh_phi");
            tree->Branch(name, &v_phoFromTICL_matchedGenPh_phi);
            
            
            sprintf(name, "phoFromTICL_superClus_E");
            tree->Branch(name, &v_phoFromTICL_superClus_E);
            
            sprintf(name, "phoFromTICL_superClus_ET");
            tree->Branch(name, &v_phoFromTICL_superClus_ET);
            
            sprintf(name, "phoFromTICL_superClus_rawE");
            tree->Branch(name, &v_phoFromTICL_superClus_rawE);
            
            sprintf(name, "phoFromTICL_superClus_rawET");
            tree->Branch(name, &v_phoFromTICL_superClus_rawET);
            
            sprintf(name, "phoFromTICL_superClus_theta");
            tree->Branch(name, &v_phoFromTICL_superClus_theta);
            
            sprintf(name, "phoFromTICL_superClus_eta");
            tree->Branch(name, &v_phoFromTICL_superClus_eta);
            
            sprintf(name, "phoFromTICL_superClus_phi");
            tree->Branch(name, &v_phoFromTICL_superClus_phi);
            
            sprintf(name, "phoFromTICL_superClus_x");
            tree->Branch(name, &v_phoFromTICL_superClus_x);
            
            sprintf(name, "phoFromTICL_superClus_y");
            tree->Branch(name, &v_phoFromTICL_superClus_y);
            
            sprintf(name, "phoFromTICL_superClus_z");
            tree->Branch(name, &v_phoFromTICL_superClus_z);
            
            sprintf(name, "phoFromTICL_superClus_r");
            tree->Branch(name, &v_phoFromTICL_superClus_r);
            
            sprintf(name, "phoFromTICL_superClus_etaWidth");
            tree->Branch(name, &v_phoFromTICL_superClus_etaWidth);
            
            sprintf(name, "phoFromTICL_superClus_phiWidth");
            tree->Branch(name, &v_phoFromTICL_superClus_phiWidth);
            
            sprintf(name, "phoFromTICL_superClus_nClus");
            tree->Branch(name, &v_phoFromTICL_superClus_nClus);
            
            sprintf(name, "phoFromTICL_superClus_nHit");
            tree->Branch(name, &v_phoFromTICL_superClus_nHit);
            
            sprintf(name, "phoFromTICL_superClus_clusMaxDR");
            tree->Branch(name, &v_phoFromTICL_superClus_clusMaxDR);
            
            
            sprintf(name, "phoFromTICL_superClus_seed_dEta");
            tree->Branch(name, &v_phoFromTICL_superClus_seed_dEta);
            
            sprintf(name, "phoFromTICL_superClus_seed_dPhi");
            tree->Branch(name, &v_phoFromTICL_superClus_seed_dPhi);
            
            
            sprintf(name, "phoFromTICL_superClus_recHit1_E");
            tree->Branch(name, &v_phoFromTICL_superClus_recHit1_E);
            
            sprintf(name, "phoFromTICL_superClus_recHit2_E");
            tree->Branch(name, &v_phoFromTICL_superClus_recHit2_E);
            
            
            sprintf(name, "phoFromTICL_superClus_sigma2rr");
            tree->Branch(name, &v_phoFromTICL_superClus_sigma2rr);
            
            sprintf(name, "phoFromTICL_superClus_sigma2etaEta");
            tree->Branch(name, &v_phoFromTICL_superClus_sigma2etaEta);
            
            sprintf(name, "phoFromTICL_superClus_sigma2phiPhi");
            tree->Branch(name, &v_phoFromTICL_superClus_sigma2phiPhi);
            
            sprintf(name, "phoFromTICL_superClus_sigma2rEta");
            tree->Branch(name, &v_phoFromTICL_superClus_sigma2rEta);
            
            sprintf(name, "phoFromTICL_superClus_sigma2rPhi");
            tree->Branch(name, &v_phoFromTICL_superClus_sigma2rPhi);
            
            sprintf(name, "phoFromTICL_superClus_sigma2etaPhi");
            tree->Branch(name, &v_phoFromTICL_superClus_sigma2etaPhi);
            
            sprintf(name, "phoFromTICL_superClus_sigma2diag1");
            tree->Branch(name, &v_phoFromTICL_superClus_sigma2diag1);
            
            sprintf(name, "phoFromTICL_superClus_sigma2diag2");
            tree->Branch(name, &v_phoFromTICL_superClus_sigma2diag2);
            
            sprintf(name, "phoFromTICL_superClus_sigma2diag3");
            tree->Branch(name, &v_phoFromTICL_superClus_sigma2diag3);
            
            
            sprintf(name, "phoFromTICL_superClusSeed_E");
            tree->Branch(name, &v_phoFromTICL_superClusSeed_E);
            
            sprintf(name, "phoFromTICL_superClusSeed_ET");
            tree->Branch(name, &v_phoFromTICL_superClusSeed_ET);
            
            sprintf(name, "phoFromTICL_superClusSeed_eta");
            tree->Branch(name, &v_phoFromTICL_superClusSeed_eta);
            
            sprintf(name, "phoFromTICL_superClusSeed_phi");
            tree->Branch(name, &v_phoFromTICL_superClusSeed_phi);
            
            sprintf(name, "phoFromTICL_superClusSeed_recHit1_E");
            tree->Branch(name, &v_phoFromTICL_superClusSeed_recHit1_E);
            
            sprintf(name, "phoFromTICL_superClusSeed_recHit2_E");
            tree->Branch(name, &v_phoFromTICL_superClusSeed_recHit2_E);
            
            
            //
            sprintf(name, "phoFromTICL_nRecHit");
            tree->Branch(name, &v_phoFromTICL_nRecHit);
            
            sprintf(name, "phoFromTICL_recHit_E");
            tree->Branch(name, &vv_phoFromTICL_recHit_E);
            
            sprintf(name, "phoFromTICL_recHit_x");
            tree->Branch(name, &vv_phoFromTICL_recHit_x);
            
            sprintf(name, "phoFromTICL_recHit_y");
            tree->Branch(name, &vv_phoFromTICL_recHit_y);
            
            sprintf(name, "phoFromTICL_recHit_z");
            tree->Branch(name, &vv_phoFromTICL_recHit_z);
            
            sprintf(name, "phoFromTICL_recHit_time");
            tree->Branch(name, &vv_phoFromTICL_recHit_time);
            
            sprintf(name, "phoFromTICL_recHit_timeError");
            tree->Branch(name, &vv_phoFromTICL_recHit_timeError);
            
            sprintf(name, "phoFromTICL_recHit_eta");
            tree->Branch(name, &vv_phoFromTICL_recHit_eta);
            
            sprintf(name, "phoFromTICL_recHit_phi");
            tree->Branch(name, &vv_phoFromTICL_recHit_phi);
            
            sprintf(name, "phoFromTICL_recHit_ET");
            tree->Branch(name, &vv_phoFromTICL_recHit_ET);
            
            sprintf(name, "phoFromTICL_recHit_detector");
            tree->Branch(name, &vv_phoFromTICL_recHit_detector);
            
            sprintf(name, "phoFromTICL_recHit_layer");
            tree->Branch(name, &vv_phoFromTICL_recHit_layer);
            
            sprintf(name, "phoFromTICL_recHit_isSimHitMatched");
            tree->Branch(name, &vv_phoFromTICL_recHit_isSimHitMatched);
            
            sprintf(name, "phoFromTICL_recHit_SCdEta");
            tree->Branch(name, &vv_phoFromTICL_recHit_SCdEta);
            
            sprintf(name, "phoFromTICL_recHit_SCdPhi");
            tree->Branch(name, &vv_phoFromTICL_recHit_SCdPhi);
            
            sprintf(name, "phoFromTICL_recHit_SCdR");
            tree->Branch(name, &vv_phoFromTICL_recHit_SCdR);
            
            
            //
            sprintf(name, "phoFromTICL_nlcInCone");
            tree->Branch(name, &v_phoFromTICL_nlcInCone);
            
            sprintf(name, "phoFromTICL_lcInCone_E");
            tree->Branch(name, &vv_phoFromTICL_lcInCone_E);
            
            sprintf(name, "phoFromTICL_lcInCone_x");
            tree->Branch(name, &vv_phoFromTICL_lcInCone_x);
            
            sprintf(name, "phoFromTICL_lcInCone_y");
            tree->Branch(name, &vv_phoFromTICL_lcInCone_y);
            
            sprintf(name, "phoFromTICL_lcInCone_z");
            tree->Branch(name, &vv_phoFromTICL_lcInCone_z);
            
            sprintf(name, "phoFromTICL_lcInCone_time");
            tree->Branch(name, &vv_phoFromTICL_lcInCone_time);
            
            sprintf(name, "phoFromTICL_lcInCone_timeError");
            tree->Branch(name, &vv_phoFromTICL_lcInCone_timeError);
            
            sprintf(name, "phoFromTICL_lcInCone_eta");
            tree->Branch(name, &vv_phoFromTICL_lcInCone_eta);
            
            sprintf(name, "phoFromTICL_lcInCone_phi");
            tree->Branch(name, &vv_phoFromTICL_lcInCone_phi);
            
            sprintf(name, "phoFromTICL_lcInCone_ET");
            tree->Branch(name, &vv_phoFromTICL_lcInCone_ET);
            
            sprintf(name, "phoFromTICL_lcInCone_size");
            tree->Branch(name, &vv_phoFromTICL_lcInCone_size);
            
            sprintf(name, "phoFromTICL_lcInCone_detector");
            tree->Branch(name, &vv_phoFromTICL_lcInCone_detector);
            
            sprintf(name, "phoFromTICL_lcInCone_layer");
            tree->Branch(name, &vv_phoFromTICL_lcInCone_layer);
            
            sprintf(name, "phoFromTICL_lcInCone_SCdEta");
            tree->Branch(name, &vv_phoFromTICL_lcInCone_SCdEta);
            
            sprintf(name, "phoFromTICL_lcInCone_SCdPhi");
            tree->Branch(name, &vv_phoFromTICL_lcInCone_SCdPhi);
            
            sprintf(name, "phoFromTICL_lcInCone_SCdR");
            tree->Branch(name, &vv_phoFromTICL_lcInCone_SCdR);
            
            //
            sprintf(name, "caloParticle_n");
            tree->Branch(name, &caloParticle_n);
            
            sprintf(name, "caloParticle_E");
            tree->Branch(name, &v_caloParticle_E);
            
            sprintf(name, "caloParticle_px");
            tree->Branch(name, &v_caloParticle_px);
            
            sprintf(name, "caloParticle_py");
            tree->Branch(name, &v_caloParticle_py);
            
            sprintf(name, "caloParticle_pz");
            tree->Branch(name, &v_caloParticle_pz);
            
            sprintf(name, "caloParticle_pT");
            tree->Branch(name, &v_caloParticle_pT);
            
            sprintf(name, "caloParticle_eta");
            tree->Branch(name, &v_caloParticle_eta);
            
            sprintf(name, "caloParticle_phi");
            tree->Branch(name, &v_caloParticle_phi);
            
            sprintf(name, "caloParticle_pdgid");
            tree->Branch(name, &v_caloParticle_pdgid);
        }
        
        
        void init_RvarContent(std::string objName, std::string suffix)
        {
            if(!objName.length())
            {
                printf("Error in TreeOutput::init_RvarContent(...): Argument \"objName\" cannot be empty. \n");
                printf("Exiting... \n");
                
                exit(EXIT_FAILURE);
            }
            
            if(!suffix.length())
            {
                printf("Error in TreeOutput::init_RvarContent(...): Argument \"suffix\" cannot be empty. \n");
                printf("Exiting... \n");
                
                exit(EXIT_FAILURE);
            }
            
            RvarContent *rVarContent = new RvarContent();
            
            
            std::string key = objName + "_" + suffix;
            
            
            //
            sprintf(name, "%s_Rvar_%s", objName.c_str(), suffix.c_str());
            tree->Branch(name, &rVarContent->v_Rvar);
            
            
            // Add to the map
            m_RvarContent[key] = rVarContent;
        }
        
        
        void init_PCAvarContent(std::string objName, std::string suffix)
        {
            if(!objName.length())
            {
                printf("Error in TreeOutput::init_PCAvarContent(...): Argument \"objName\" cannot be empty. \n");
                printf("Exiting... \n");
                
                exit(EXIT_FAILURE);
            }
            
            if(!suffix.length())
            {
                printf("Error in TreeOutput::init_PCAvarContent(...): Argument \"suffix\" cannot be empty. \n");
                printf("Exiting... \n");
                
                exit(EXIT_FAILURE);
            }
            
            PCAvarContent *pcaVarContent = new PCAvarContent();
            
            
            std::string key = objName + "_" + suffix;
            
            
            //
            sprintf(name, "%s_sigma2uu_%s", objName.c_str(), suffix.c_str());
            tree->Branch(name, &pcaVarContent->v_sigma2uu);
            
            sprintf(name, "%s_sigma2vv_%s", objName.c_str(), suffix.c_str());
            tree->Branch(name, &pcaVarContent->v_sigma2vv);
            
            sprintf(name, "%s_sigma2ww_%s", objName.c_str(), suffix.c_str());
            tree->Branch(name, &pcaVarContent->v_sigma2ww);
            
            
            // Add to the map
            m_PCAvarContent[key] = pcaVarContent;
        }
        
        
        void init_isoVarContent(std::string objName, std::string suffix)
        {
            if(!objName.length())
            {
                printf("Error in TreeOutput::init_isoVarContent(...): Argument \"objName\" cannot be empty. \n");
                printf("Exiting... \n");
                
                exit(EXIT_FAILURE);
            }
            
            if(!suffix.length())
            {
                printf("Error in TreeOutput::init_isoVarContent(...): Argument \"suffix\" cannot be empty. \n");
                printf("Exiting... \n");
                
                exit(EXIT_FAILURE);
            }
            
            IsoVarContent *isoVarContent = new IsoVarContent();
            
            
            std::string key = objName + "_" + suffix;
            
            
            //
            sprintf(name, "%s_iso_sumETratio_%s", objName.c_str(), suffix.c_str());
            tree->Branch(name, &isoVarContent->v_iso_sumETratio);
            
            
            sprintf(name, "%s_iso_trackSumPt_%s", objName.c_str(), suffix.c_str());
            tree->Branch(name, &isoVarContent->v_iso_trackSumPt);
            
            
            //
            sprintf(name, "%s_HoverE_%s", objName.c_str(), suffix.c_str());
            tree->Branch(name, &isoVarContent->v_HoverE);
            
            
            // Add to the map
            m_isoVarContent[key] = isoVarContent;
        }
        
        
        void addToCustomVarMap(std::string objName, std::vector <std::string> v_key)
        {
            if(!objName.length())
            {
                printf("Error in TreeOutput::addToCustomVarMap(...): Argument \"objName\" cannot be empty. \n");
                printf("Exiting... \n");
                
                exit(EXIT_FAILURE);
            }
            
            for(auto &key : v_key)
            {
                std::string key_mod = objName + "_" + key;
                printf("TreeOutput::addToCustomVarMap(...): Adding key \"%s\". \n", key_mod.c_str());
                
                if(m_customVarContent.find(key_mod) != m_customVarContent.end())
                {
                    printf("Error in TreeOutput::addToCustomVarMap(...): Key \"%s\" already exists in map. \n", key_mod.c_str());
                    printf("Exiting... \n");
                }
                
                m_customVarContent[key_mod] = {};
                
                sprintf(name, "%s", key_mod.c_str());
                tree->Branch(name, &m_customVarContent[key_mod]);
            }
        }
        
        
        void fill()
        {
            tree->Fill();
        }
        
        
        void clear()
        {
            // Gen electron //
            genEl_n = 0;
            v_genEl_E.clear();
            v_genEl_px.clear();
            v_genEl_py.clear();
            v_genEl_pz.clear();
            v_genEl_pT.clear();
            v_genEl_eta.clear();
            v_genEl_phi.clear();
            
            v_genEl_HGCalEEP_EsortedIndex.clear();
            v_genEl_HGCalEEM_EsortedIndex.clear();
            
            v_genEl_multiClus_totE.clear();
            
            v_genEl_multiClus_n.clear();
            v_genEl_nearestMultiClusEnRatio.clear();
            v_genEl_multiClusEnRatio.clear();
            
            
            // Gen photon //
            genPh_n = 0;
            v_genPh_E.clear();
            v_genPh_px.clear();
            v_genPh_py.clear();
            v_genPh_pz.clear();
            v_genPh_pT.clear();
            v_genPh_eta.clear();
            v_genPh_phi.clear();
            
            v_genPh_HGCalEEP_EsortedIndex.clear();
            v_genPh_HGCalEEM_EsortedIndex.clear();
            
            v_genPh_multiClus_totE.clear();
            
            v_genPh_HGCalEEP_deltaR.clear();
            v_genPh_HGCalEEM_deltaR.clear();
            
            
            // Pileup //
            pileup_n = 0;
            
            
            // Rho //
            rho = 0;
            
            
            nHit_EcalEB = 0;
            nHit_HGCEE = 0;
            nHit_HGCHEF = 0;
            nHit_HGCHEB = 0;
            
            
            // HGCAL layer clusters //
            HGCALlayerClus_n = 0;
            v_HGCALlayerClus_E.clear();
            v_HGCALlayerClus_x.clear();
            v_HGCALlayerClus_y.clear();
            v_HGCALlayerClus_z.clear();
            v_HGCALlayerClus_eta.clear();
            v_HGCALlayerClus_phi.clear();
            v_HGCALlayerClus_ET.clear();
            v_HGCALlayerClus_time.clear();
            v_HGCALlayerClus_timeError.clear();
            v_HGCALlayerClus_detector.clear();
            v_HGCALlayerClus_layer.clear();
            
            
            // Tracksters
            trackster_n = 0;
            v_trackster_E.clear();
            v_trackster_x.clear();
            v_trackster_y.clear();
            v_trackster_z.clear();
            v_trackster_eta.clear();
            v_trackster_phi.clear();
            v_trackster_ET.clear();
            
            
            // MultiClusters
            multiClus_n = 0;
            v_multiClus_genElIndex.clear();
            v_multiClus_E.clear();
            v_multiClus_x.clear();
            v_multiClus_y.clear();
            v_multiClus_z.clear();
            v_multiClus_eta.clear();
            v_multiClus_phi.clear();
            v_multiClus_ET.clear();
            
            v_multiClus_corrE.clear();
            v_multiClus_corrET.clear();
            
            v_multiClus_dX.clear();
            v_multiClus_dY.clear();
            v_multiClus_dZ.clear();
            
            v_multiClus_dEta.clear();
            v_multiClus_dPhi.clear();
            
            v_multiClus_sigma2rr.clear();
            v_multiClus_sigma2etaEta.clear();
            v_multiClus_sigma2phiPhi.clear();
            
            v_multiClus_sigma2rEta.clear();
            v_multiClus_sigma2rPhi.clear();
            v_multiClus_sigma2etaPhi.clear();
            
            v_multiClus_sigma2diag1.clear();
            v_multiClus_sigma2diag2.clear();
            v_multiClus_sigma2diag3.clear();
            
            v_multiClus_EsortedIndex.clear();
            v_multiClus_HGCalEEP_EsortedIndex.clear();
            v_multiClus_HGCalEEM_EsortedIndex.clear();
            
            
            v_multiClus_clus_n.clear();
            v_multiClus_clus_startIndex.clear();
            v_multiClus_clus_E.clear();
            v_multiClus_clus_x.clear();
            v_multiClus_clus_y.clear();
            v_multiClus_clus_z.clear();
            v_multiClus_clus_eta.clear();
            v_multiClus_clus_phi.clear();
            v_multiClus_clus_ET.clear();
            
            v_multiClus_clus_layer.clear();
            v_multiClus_clus_multiplicity.clear();
            v_multiClus_clus_nHit.clear();
            
            
            v_multiClus_uniqueClus_n.clear();
            v_multiClus_uniqueClus_E.clear();
            v_multiClus_uniqueClus_x.clear();
            v_multiClus_uniqueClus_y.clear();
            v_multiClus_uniqueClus_z.clear();
            v_multiClus_uniqueClus_eta.clear();
            v_multiClus_uniqueClus_phi.clear();
            v_multiClus_uniqueClus_ET.clear();
            
            v_multiClus_uniqueClus_layer.clear();
            v_multiClus_uniqueClus_multiplicity.clear();
            v_multiClus_uniqueClus_nHit.clear();
            
            v_multiClus_mc1_dX.clear();
            v_multiClus_mc1_dY.clear();
            v_multiClus_mc1_dZ.clear();
            
            v_multiClus_mc1_dEta.clear();
            v_multiClus_mc1_dPhi.clear();
            v_multiClus_mc1_dR.clear();
            
            
            multiClus_HGCalEEP_meanX = 0;
            multiClus_HGCalEEP_meanY = 0;
            multiClus_HGCalEEP_meanZ = 0;
            multiClus_HGCalEEP_totE = 0;
            multiClus_HGCalEEP_totET = 0;
            
            multiClus_HGCalEEP_diag1 = 0;
            multiClus_HGCalEEP_diag2 = 0;
            multiClus_HGCalEEP_diag3 = 0;
            
            
            multiClus_HGCalEEM_meanX = 0;
            multiClus_HGCalEEM_meanY = 0;
            multiClus_HGCalEEM_meanZ = 0;
            multiClus_HGCalEEM_totE = 0;
            multiClus_HGCalEEM_totET = 0;
            
            multiClus_HGCalEEM_diag1 = 0;
            multiClus_HGCalEEM_diag2 = 0;
            multiClus_HGCalEEM_diag3 = 0;
            
            
            for(int iLayer = 0; iLayer < Constants::HGCalEE_nLayer; iLayer++)
            {
                v_simHit_HGCalEEPlayer_totE.at(iLayer) = 0;
                v_recHit_HGCalEEPlayer_totE.at(iLayer) = 0;
                
                v_simHit_HGCalEEMlayer_totE.at(iLayer) = 0;
                v_recHit_HGCalEEMlayer_totE.at(iLayer) = 0;
            }
            
            
            //
            simHit_n = 0;
            v_simHit_E.clear();
            v_simHit_x.clear();
            v_simHit_y.clear();
            v_simHit_z.clear();
            v_simHit_eta.clear();
            v_simHit_phi.clear();
            v_simHit_ET.clear();
            v_simHit_layer.clear();
            v_simHit_zside.clear();
            v_simHit_isCaloParticleMatched.clear();
            
            
            //
            recHit_n = 0;
            v_recHit_E.clear();
            v_recHit_x.clear();
            v_recHit_y.clear();
            v_recHit_z.clear();
            v_recHit_eta.clear();
            v_recHit_phi.clear();
            v_recHit_ET.clear();
            v_recHit_layer.clear();
            v_recHit_zside.clear();
            v_recHit_matchedSimHitIndex.clear();
            v_recHit_matchedHGCALlayerClusIndex.clear();
            v_recHit_isMultiClusMatched.clear();
            v_recHit_isCaloParticleMatched.clear();
            v_recHit_iType.clear();
            v_recHit_iCell1.clear();
            v_recHit_iCell2.clear();
            
            v_recHit_SiThickness.clear();
            
            
            //
            gsfEleFromMultiClus_n = 0;
            v_gsfEleFromMultiClus_E.clear();
            v_gsfEleFromMultiClus_px.clear();
            v_gsfEleFromMultiClus_py.clear();
            v_gsfEleFromMultiClus_pz.clear();
            v_gsfEleFromMultiClus_pT.clear();
            v_gsfEleFromMultiClus_eta.clear();
            v_gsfEleFromMultiClus_phi.clear();
            
            v_gsfEleFromMultiClus_genEl_minDeltaR.clear();
            v_gsfEleFromMultiClus_matchedGenEl_E.clear();
            
            
            //
            gsfEleFromTICL_n = 0;
            v_gsfEleFromTICL_E.clear();
            v_gsfEleFromTICL_px.clear();
            v_gsfEleFromTICL_py.clear();
            v_gsfEleFromTICL_pz.clear();
            v_gsfEleFromTICL_pT.clear();
            v_gsfEleFromTICL_eta.clear();
            v_gsfEleFromTICL_phi.clear();
            v_gsfEleFromTICL_ET.clear();
            
            v_gsfEleFromTICL_genEl_minDeltaR.clear();
            v_gsfEleFromTICL_nearestGenEl_idx.clear();
            v_gsfEleFromTICL_matchedGenEl_E.clear();
            v_gsfEleFromTICL_matchedGenEl_pT.clear();
            v_gsfEleFromTICL_matchedGenEl_eta.clear();
            v_gsfEleFromTICL_matchedGenEl_phi.clear();
            
            v_gsfEleFromTICL_gsfTrack_p.clear();
            v_gsfEleFromTICL_gsfTrack_px.clear();
            v_gsfEleFromTICL_gsfTrack_py.clear();
            v_gsfEleFromTICL_gsfTrack_pz.clear();
            v_gsfEleFromTICL_gsfTrack_pT.clear();
            v_gsfEleFromTICL_gsfTrack_eta.clear();
            v_gsfEleFromTICL_gsfTrack_phi.clear();
            
            v_gsfEleFromTICL_trkAtVtx_p.clear();
            v_gsfEleFromTICL_trkAtVtx_px.clear();
            v_gsfEleFromTICL_trkAtVtx_py.clear();
            v_gsfEleFromTICL_trkAtVtx_pz.clear();
            v_gsfEleFromTICL_trkAtVtx_pT.clear();
            v_gsfEleFromTICL_trkAtVtx_eta.clear();
            v_gsfEleFromTICL_trkAtVtx_phi.clear();
            
            v_gsfEleFromTICL_dr03TkSumPt.clear();
            
            v_gsfEleFromTICL_superClus_E.clear();
            v_gsfEleFromTICL_superClus_ET.clear();
            v_gsfEleFromTICL_superClus_rawE.clear();
            v_gsfEleFromTICL_superClus_rawET.clear();
            v_gsfEleFromTICL_superClus_theta.clear();
            v_gsfEleFromTICL_superClus_eta.clear();
            v_gsfEleFromTICL_superClus_phi.clear();
            v_gsfEleFromTICL_superClus_x.clear();
            v_gsfEleFromTICL_superClus_y.clear();
            v_gsfEleFromTICL_superClus_z.clear();
            v_gsfEleFromTICL_superClus_r.clear();
            v_gsfEleFromTICL_superClus_etaWidth.clear();
            v_gsfEleFromTICL_superClus_phiWidth.clear();
            v_gsfEleFromTICL_superClus_nClus.clear();
            v_gsfEleFromTICL_superClus_nHit.clear();
            v_gsfEleFromTICL_superClus_nearestCellDist.clear();
            v_gsfEleFromTICL_superClus_cellNeighbour1ringWindow_n.clear();
            v_gsfEleFromTICL_superClus_cellNeighbour2ringWindow_n.clear();
            v_gsfEleFromTICL_superClus_clusMaxDR.clear();
            
            v_gsfEleFromTICL_superClus_seed_dEta.clear();
            v_gsfEleFromTICL_superClus_seed_dPhi.clear();
            
            v_gsfEleFromTICL_superClus_recHit1_E.clear();
            v_gsfEleFromTICL_superClus_recHit2_E.clear();
            
            v_gsfEleFromTICL_superClus_sigma2etaEtaLW.clear();
            v_gsfEleFromTICL_superClus_sigma2phiPhiLW.clear();
            
            v_gsfEleFromTICL_superClus_sigma2rr.clear();
            v_gsfEleFromTICL_superClus_sigma2etaEta.clear();
            v_gsfEleFromTICL_superClus_sigma2phiPhi.clear();
            
            v_gsfEleFromTICL_superClus_sigma2rEta.clear();
            v_gsfEleFromTICL_superClus_sigma2rPhi.clear();
            v_gsfEleFromTICL_superClus_sigma2etaPhi.clear();
            
            v_gsfEleFromTICL_superClus_sigma2diag1.clear();
            v_gsfEleFromTICL_superClus_sigma2diag2.clear();
            v_gsfEleFromTICL_superClus_sigma2diag3.clear();
            
            v_gsfEleFromTICL_E7.clear();
            v_gsfEleFromTICL_R7.clear();
            
            v_gsfEleFromTICL_E19.clear();
            v_gsfEleFromTICL_R19.clear();
            
            v_gsfEleFromTICL_E2p4.clear();
            v_gsfEleFromTICL_R2p4.clear();
            
            v_gsfEleFromTICL_E2p6.clear();
            v_gsfEleFromTICL_R2p6.clear();
            
            v_gsfEleFromTICL_E2p8.clear();
            v_gsfEleFromTICL_R2p8.clear();
            
            v_gsfEleFromTICL_E3p0.clear();
            v_gsfEleFromTICL_R3p0.clear();
            
            v_gsfEleFromTICL_E3p5.clear();
            v_gsfEleFromTICL_R3p5.clear();
            
            v_gsfEleFromTICL_superClusSeed_E.clear();
            v_gsfEleFromTICL_superClusSeed_ET.clear();
            v_gsfEleFromTICL_superClusSeed_eta.clear();
            v_gsfEleFromTICL_superClusSeed_phi.clear();
            
            v_gsfEleFromTICL_superClusSeed_recHit1_E.clear();
            v_gsfEleFromTICL_superClusSeed_recHit2_E.clear();
            
            vv_gsfEleFromTICL_superClusSeed_clus_dEta.clear();
            vv_gsfEleFromTICL_superClusSeed_clus_dPhi.clear();
            
            v_gsfEleFromTICL_superClus_TICLclus_n.clear();
            
            vv_gsfEleFromTICL_superClus_TICLclus_E.clear();
            vv_gsfEleFromTICL_superClus_TICLclus_ET.clear();
            vv_gsfEleFromTICL_superClus_TICLclus_nClus.clear();
            
            vv_gsfEleFromTICL_superClusSeed_TICLclus_dX.clear();
            vv_gsfEleFromTICL_superClusSeed_TICLclus_dY.clear();
            vv_gsfEleFromTICL_superClusSeed_TICLclus_dZ.clear();
            
            vv_gsfEleFromTICL_superClusSeed_TICLclus_dEta.clear();
            vv_gsfEleFromTICL_superClusSeed_TICLclus_dPhi.clear();
            vv_gsfEleFromTICL_superClusSeed_TICLclus_dR.clear();
            
            
            for(int iLayer = 0; iLayer < Constants::HGCalEE_nLayer; iLayer++)
            {
                vv_gsfEleFromTICL_superClus_E_layer.at(iLayer).clear();
                
                vv_gsfEleFromTICL_E2p4_layer.at(iLayer).clear();
                vv_gsfEleFromTICL_R2p4_layer.at(iLayer).clear();
                
                vv_gsfEleFromTICL_E2p6_layer.at(iLayer).clear();
                vv_gsfEleFromTICL_R2p6_layer.at(iLayer).clear();
                
                vv_gsfEleFromTICL_E2p8_layer.at(iLayer).clear();
                vv_gsfEleFromTICL_R2p8_layer.at(iLayer).clear();
                
                vv_gsfEleFromTICL_E3p0_layer.at(iLayer).clear();
                vv_gsfEleFromTICL_R3p0_layer.at(iLayer).clear();
                
                vv_gsfEleFromTICL_E3p5_layer.at(iLayer).clear();
                vv_gsfEleFromTICL_R3p5_layer.at(iLayer).clear();
            }
            
            
            for(auto iter = m_RvarContent.begin(); iter != m_RvarContent.end(); iter++)
            {
                iter->second->clear();
            }
            
            for(auto iter = m_PCAvarContent.begin(); iter != m_PCAvarContent.end(); iter++)
            {
                iter->second->clear();
            }
            
            for(auto iter = m_isoVarContent.begin(); iter != m_isoVarContent.end(); iter++)
            {
                iter->second->clear();
            }
            
            
            // TDR photons
            phoFromMultiClus_n = 0;
            v_phoFromMultiClus_E.clear();
            v_phoFromMultiClus_px.clear();
            v_phoFromMultiClus_py.clear();
            v_phoFromMultiClus_pz.clear();
            v_phoFromMultiClus_pT.clear();
            v_phoFromMultiClus_eta.clear();
            v_phoFromMultiClus_phi.clear();
            
            v_phoFromMultiClus_genPh_minDeltaR.clear();
            v_phoFromMultiClus_matchedGenPh_E.clear();
            
            
            // TICL photons
            phoFromTICL_n = 0;
            v_phoFromTICL_E.clear();
            v_phoFromTICL_px.clear();
            v_phoFromTICL_py.clear();
            v_phoFromTICL_pz.clear();
            v_phoFromTICL_pT.clear();
            v_phoFromTICL_eta.clear();
            v_phoFromTICL_phi.clear();
            v_phoFromTICL_ET.clear();
            
            v_phoFromTICL_genPh_minDeltaR.clear();
            v_phoFromTICL_nearestGenPh_idx.clear();
            v_phoFromTICL_matchedGenPh_E.clear();
            v_phoFromTICL_matchedGenPh_pT.clear();
            v_phoFromTICL_matchedGenPh_eta.clear();
            v_phoFromTICL_matchedGenPh_phi.clear();
            
            v_phoFromTICL_superClus_E.clear();
            v_phoFromTICL_superClus_ET.clear();
            v_phoFromTICL_superClus_rawE.clear();
            v_phoFromTICL_superClus_rawET.clear();
            v_phoFromTICL_superClus_theta.clear();
            v_phoFromTICL_superClus_eta.clear();
            v_phoFromTICL_superClus_phi.clear();
            v_phoFromTICL_superClus_x.clear();
            v_phoFromTICL_superClus_y.clear();
            v_phoFromTICL_superClus_z.clear();
            v_phoFromTICL_superClus_r.clear();
            v_phoFromTICL_superClus_etaWidth.clear();
            v_phoFromTICL_superClus_phiWidth.clear();
            v_phoFromTICL_superClus_nClus.clear();
            v_phoFromTICL_superClus_nHit.clear();
            v_phoFromTICL_superClus_clusMaxDR.clear();
            
            v_phoFromTICL_superClus_seed_dEta.clear();
            v_phoFromTICL_superClus_seed_dPhi.clear();
            
            v_phoFromTICL_superClus_recHit1_E.clear();
            v_phoFromTICL_superClus_recHit2_E.clear();
            
            v_phoFromTICL_superClus_sigma2rr.clear();
            v_phoFromTICL_superClus_sigma2etaEta.clear();
            v_phoFromTICL_superClus_sigma2phiPhi.clear();
            
            v_phoFromTICL_superClus_sigma2rEta.clear();
            v_phoFromTICL_superClus_sigma2rPhi.clear();
            v_phoFromTICL_superClus_sigma2etaPhi.clear();
            
            v_phoFromTICL_superClus_sigma2diag1.clear();
            v_phoFromTICL_superClus_sigma2diag2.clear();
            v_phoFromTICL_superClus_sigma2diag3.clear();
            
            v_phoFromTICL_superClusSeed_E.clear();
            v_phoFromTICL_superClusSeed_ET.clear();
            v_phoFromTICL_superClusSeed_eta.clear();
            v_phoFromTICL_superClusSeed_phi.clear();
            
            v_phoFromTICL_superClusSeed_recHit1_E.clear();
            v_phoFromTICL_superClusSeed_recHit2_E.clear();
            
            v_phoFromTICL_nRecHit.clear();
            vv_phoFromTICL_recHit_E.clear();
            vv_phoFromTICL_recHit_x.clear();
            vv_phoFromTICL_recHit_y.clear();
            vv_phoFromTICL_recHit_z.clear();
            vv_phoFromTICL_recHit_time.clear();
            vv_phoFromTICL_recHit_timeError.clear();
            vv_phoFromTICL_recHit_eta.clear();
            vv_phoFromTICL_recHit_phi.clear();
            vv_phoFromTICL_recHit_ET.clear();
            vv_phoFromTICL_recHit_detector.clear();
            vv_phoFromTICL_recHit_layer.clear();
            vv_phoFromTICL_recHit_isSimHitMatched.clear();
            vv_phoFromTICL_recHit_SCdEta.clear();
            vv_phoFromTICL_recHit_SCdPhi.clear();
            vv_phoFromTICL_recHit_SCdR.clear();
            
            v_phoFromTICL_nlcInCone.clear();
            vv_phoFromTICL_lcInCone_E.clear();
            vv_phoFromTICL_lcInCone_x.clear();
            vv_phoFromTICL_lcInCone_y.clear();
            vv_phoFromTICL_lcInCone_z.clear();
            vv_phoFromTICL_lcInCone_time.clear();
            vv_phoFromTICL_lcInCone_timeError.clear();
            vv_phoFromTICL_lcInCone_eta.clear();
            vv_phoFromTICL_lcInCone_phi.clear();
            vv_phoFromTICL_lcInCone_ET.clear();
            vv_phoFromTICL_lcInCone_size.clear();
            vv_phoFromTICL_lcInCone_detector.clear();
            vv_phoFromTICL_lcInCone_layer.clear();
            vv_phoFromTICL_lcInCone_SCdEta.clear();
            vv_phoFromTICL_lcInCone_SCdPhi.clear();
            vv_phoFromTICL_lcInCone_SCdR.clear();
            
            
            //
            caloParticle_n = 0;
            v_caloParticle_E.clear();
            v_caloParticle_px.clear();
            v_caloParticle_py.clear();
            v_caloParticle_pz.clear();
            v_caloParticle_pT.clear();
            v_caloParticle_eta.clear();
            v_caloParticle_phi.clear();
            v_caloParticle_pdgid.clear();
            
            
            for(auto &v : m_customVarContent)
            {
                v.second.clear();
            }
        }
    };
}


# endif
