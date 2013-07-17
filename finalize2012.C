#include <TTree.h>
#include <TFile.h>
#include <TLegend.h>
#include <TF1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TCanvas.h>
#include <TMath.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <TLorentzVector.h>
#include <fillPlot2012.h>

inline double delta_phi(double phi1, double phi2) {
  
  double dphi = TMath::Abs(phi1 - phi2);
  return (dphi <= TMath::Pi())? dphi : TMath::TwoPi() - dphi;
}

double probability_disc(double obs, double exp){
  
  double prob=0.;
  TF1 *mypoisson = new TF1("mypoisson","TMath::Poisson(x,[0])",0,10000);
  mypoisson->SetParameter(0, exp);
  if (exp<100)
    prob= mypoisson->Integral(obs-0.5,999);
  else
    prob= mypoisson->Integral(obs-0.5,9999);
  return TMath::NormQuantile(1-prob);
}

vector <double> finalize2012(double int_exp_2010, double int_exp_2012, double pt1=50, double pt2=30, double pthiggsmin = -100, double pthiggsmax = -100, double ptj1=20, double ptj2=15, double misset=70, double deltae=2.5, double zep=2.5, double mjj=300, double deltap = 2.6, double jetmet=0., double p1met=0., double p2met=0., double hmet=0., double phigg=300., int eb = 1, int r9 = 1, int cic = 4, bool usepfcic = true, bool thirdcat = 0, bool leptontag = 0, bool leptonveto = 0, bool vbfveto = 0, bool inclusive = 0, string variableMC = "massgg", string variableData = "massgg", int nbin = 200, double min = 90, double max = 190, string axis = "m(#gamma#gamma)[GeV]"){
  
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  gStyle->SetOptStat(1111111);  
  gStyle->SetOptFit(111110); 
  gStyle->SetOptFile(1); 
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerSize(.3);
  gStyle->SetMarkerColor(1);
  
  TCanvas* c0 = new TCanvas("c0"," ",200,10,500,500);
  c0->Clear();
  
  // input files
  string mcnames[33];
  mcnames[0]  = "box, 10-25";
  mcnames[1]  = "box, 25-250";
  mcnames[2]  = "box, >250";
  mcnames[3]  = "diphoton, 10-25";
  mcnames[4]  = "diphoton, 25-250";
  mcnames[5]  = "diphoton, >250";
  mcnames[6]  = "gjet, 20-40";
  mcnames[7]  = "gjet, >40";
  mcnames[8]  = "qcdpt>40";
  mcnames[9]  = "qcd30<pt<40";
  mcnames[10] = "dy";
  mcnames[11] = "WenuG";
  mcnames[12] = "WmnuG";
  mcnames[13] = "WtnuG";
  mcnames[14] = "ZeeG";
  mcnames[15] = "ZmmG";
  mcnames[16] = "ZttG";
  mcnames[17] = "ggWminus";
  mcnames[18] = "ggWplus";
  mcnames[19] = "ggZ";
  mcnames[20] = "ggtt";
  mcnames[21] = "ttjets";
  mcnames[22] = "Wjets";
  mcnames[23] = "WW";
  mcnames[24] = "WZ";
  mcnames[25] = "ZZ";
  mcnames[26] = "higgsgluglu";
  mcnames[27] = "higgsVBF";
  mcnames[28] = "higgsWH, Wleptonic";
  mcnames[29] = "higgsWH, Whadronic";
  mcnames[30] = "higgsZH, Zleptonic";
  mcnames[31] = "higgsZH, Zhadronic";
  mcnames[32] = "higgsZH, Zneutrinos";

  TFile* mc_2012[33];            // [33] = MC samples
  TFile* mc_gluglu_2012[11];     // [11] = masses
  TFile* mc_vbf_2012[11];
  TFile* mc_wzh_2012[11];
  int h_masses[11] = {100,105,110,115,120,125,130,135,140,145,150};

  TString redntpDir= "/Users/crovelli/Work/dati/Hgg";
  
  TString preselectionLevel;
  TString preselectionLevelData;
  // preselectionLevel="preselectionMVA";
  preselectionLevel="cicpfloose";
  preselectionLevelData="cicpfloose";
  // preselectionLevel="cicloosenoeleveto";
  TString preselectionLevelCS="preselectionCS";

  // full 2012
  TFile* data = TFile::Open(redntpDir+"/redntp.52xv5_data."+preselectionLevelData+".paper-Hgg-scale.fixPerEleMva/merged/fullIchep.root");   
  // TFile* datacs = TFile::Open(redntpDir+"/redntp.52xv5_data.preselectionCS.paper-Hgg-scale.v1/merged/redntp_Photon-Run2012A-DoublePhoton-Run2012B-PromptReco-v1.root");
  
  // mc signal and backgrounds   
  if(int_exp_2012>0){
    
    // box samples
    mc_2012[0] = TFile::Open(redntpDir+"/redntp.52xv5."+preselectionLevel+".paper-Hgg-scale.fixPerEleMva/merged/redntp_DiPhotonBox_Pt-10To25_8TeV-pythia6.root");
    mc_2012[1] = TFile::Open(redntpDir+"/redntp.52xv5."+preselectionLevel+".paper-Hgg-scale.fixPerEleMva/merged/redntp_DiPhotonBox_Pt-25To250_8TeV-pythia6.root");
    mc_2012[2] = TFile::Open(redntpDir+"/redntp.52xv5."+preselectionLevel+".paper-Hgg-scale.fixPerEleMva/merged/redntp_DiPhotonBox_Pt-250ToInf_8TeV-pythia6.root");

    // diphoton jets samples 
    // mc_2012[3] = TFile::Open(redntpDir+"/redntp.52xv5."+preselectionLevel+".paper-Hgg-scale.fixPerEleMva/merged/redntp_DiPhotonBorn_Pt-10To25_8TeV-pythia6.root");
    // mc_2012[4] = TFile::Open(redntpDir+"/redntp.52xv5."+preselectionLevel+".paper-Hgg-scale.fixPerEleMva/merged/redntp_DiPhotonBorn_Pt-25To250_8TeV-pythia6.root");
    // mc_2012[5] = TFile::Open(redntpDir+"/redntp.52xv5."+preselectionLevel+".paper-Hgg-scale.fixPerEleMva/merged/redntp_DiPhotonBorn_Pt-250ToInf_8TeV-pythia6.root");
    // chiara: per non riscrivere tutto metto tre volte lo stesso sample e poi la sezione d'urto a zero per 2 dei 3
    mc_2012[3] = TFile::Open(redntpDir+"/redntp.52xv5."+preselectionLevel+".paper-Hgg-scale.fixPerEleMva/merged/redntp_DiPhotonJets_8TeV-madgraph-tarball-v2.root");
    mc_2012[4] = TFile::Open(redntpDir+"/redntp.52xv5."+preselectionLevel+".paper-Hgg-scale.fixPerEleMva/merged/redntp_DiPhotonJets_8TeV-madgraph-tarball-v2.root");
    mc_2012[5] = TFile::Open(redntpDir+"/redntp.52xv5."+preselectionLevel+".paper-Hgg-scale.fixPerEleMva/merged/redntp_DiPhotonJets_8TeV-madgraph-tarball-v2.root");

    // gjet samples
    mc_2012[6] = TFile::Open(redntpDir+"/redntp.52xv5."+preselectionLevel+".paper-Hgg-scale.fixPerEleMva/merged/redntp_GJet_Pt-20to40_doubleEMEnriched_TuneZ2star_8TeV-pythia6.root");
    mc_2012[7] = TFile::Open(redntpDir+"/redntp.52xv5."+preselectionLevel+".paper-Hgg-scale.fixPerEleMva/merged/redntp_GJet_Pt40_doubleEMEnriched_TuneZ2star_8TeV-pythia6.root");

    // qcd pt>40 samples
    mc_2012[8] = TFile::Open(redntpDir+"/redntp.52xv5."+preselectionLevel+".paper-Hgg-scale.fixPerEleMva/merged/redntp_QCD_Pt-40_doubleEMEnriched_TuneZ2star_8TeV-pythia6.root");
    // qcd 30<pt<40 samples
    mc_2012[9] = TFile::Open(redntpDir+"/redntp.52xv5."+preselectionLevel+".paper-Hgg-scale.fixPerEleMva/merged/redntp_QCD_Pt-30to40_doubleEMEnriched_TuneZ2star_8TeV-pythia6.root");

    // drell yan samples                                                                         
    mc_2012[10] = TFile::Open(redntpDir+"/redntp.52xv5."+preselectionLevel+".paper-Hgg-scale.fixPerEleMva/merged/redntp_DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball.root");

    // W+gamma, W->lnu 
    // chiara: quest'anno WG->lnuG sono tutti assieme. Per non rifare tutto, metto 3 volte lo stesso e le altre x-sec a zero    
    // prendo la x-sec da prep
    mc_2012[11] = TFile::Open(redntpDir+"/redntp.52xv5."+preselectionLevel+".paper-Hgg-scale.fixPerEleMva/merged/redntp_WGToLNuG_TuneZ2star_8TeV-madgraph-tauola.root");
    // W+gamma - dummy 
    mc_2012[12] = TFile::Open(redntpDir+"/redntp.52xv5."+preselectionLevel+".paper-Hgg-scale.fixPerEleMva/merged/redntp_WGToLNuG_TuneZ2star_8TeV-madgraph-tauola.root");
    // W+gamma - dummy
    mc_2012[13] = TFile::Open(redntpDir+"/redntp.52xv5."+preselectionLevel+".paper-Hgg-scale.fixPerEleMva/merged/redntp_WGToLNuG_TuneZ2star_8TeV-madgraph-tauola.root");

    // Z+gamma, Z->ll
    // chiara: quest'anno ZG->inclusive Per non rifare tutto, metto 3 volte lo stesso e le altre x-sec a zero
    // prendo la x-sec da prep
    mc_2012[14]  = TFile::Open(redntpDir+"/redntp.52xv5."+preselectionLevel+".paper-Hgg-scale.fixPerEleMva/merged/redntp_ZG_Inclusive_8TeV-madgraph.root");
    // Z+gamma, dummy
    mc_2012[15] = TFile::Open(redntpDir+"/redntp.52xv5."+preselectionLevel+".paper-Hgg-scale.fixPerEleMva/merged/redntp_ZG_Inclusive_8TeV-madgraph.root");
    // Z+gamma, dummy
    mc_2012[16] = TFile::Open(redntpDir+"/redntp.52xv5."+preselectionLevel+".paper-Hgg-scale.fixPerEleMva/merged/redntp_ZG_Inclusive_8TeV-madgraph.root");

    // W-gg                
    mc_2012[17] = TFile::Open(redntpDir+"/redntp.52xv5."+preselectionLevel+".paper-Hgg-scale.fixPerEleMva/merged/redntp_WmGG_cmkuo.root");
    // W+ gg               
    mc_2012[18] = TFile::Open(redntpDir+"/redntp.52xv5."+preselectionLevel+".paper-Hgg-scale.fixPerEleMva/merged/redntp_WpGG-cmkuo.root");
    // Zgg                 
    mc_2012[19] = TFile::Open(redntpDir+"/redntp.52xv5."+preselectionLevel+".paper-Hgg-scale.fixPerEleMva/merged/redntp_ZGG-cmkuo.root");
    // ttgg                
    mc_2012[20] = TFile::Open(redntpDir+"/redntp.52xv5."+preselectionLevel+".paper-Hgg-scale.fixPerEleMva/merged/redntp_TTbarGG_0Jet_S1-cmkuo-TTGG_525_RECO_s46_v1.root");

    // TTjets                                                                                              
    mc_2012[21] = TFile::Open(redntpDir+"/redntp.52xv5."+preselectionLevel+".paper-Hgg-scale.fixPerEleMva/merged/redntp_TTJets_TuneZ2star_8TeV-madgraph-tauola.root");
    // Wjets                                                                                                      
    mc_2012[22] = TFile::Open(redntpDir+"/redntp.52xv5."+preselectionLevel+".paper-Hgg-scale.fixPerEleMva/merged/redntp_WJetsToLNu_TuneZ2Star_8TeV-madgraph.root");

    // WW                                                                                                                
    mc_2012[23] = TFile::Open(redntpDir+"/redntp.52xv5."+preselectionLevel+".paper-Hgg-scale.fixPerEleMva/merged/redntp_WWJetsTo2L2Nu_TuneZ2star_8TeV-madgraph-tauola.root");  // chiara
    // WZ                                                                                                                
    mc_2012[24] = TFile::Open(redntpDir+"/redntp.52xv5."+preselectionLevel+".paper-Hgg-scale.fixPerEleMva/merged/redntp_WZTo3LNu_TuneZ2star_8TeV_pythia6_tauola.root");
    // ZZ
    mc_2012[25] = TFile::Open(redntpDir+"/redntp.52xv5."+preselectionLevel+".paper-Hgg-scale.fixPerEleMva/merged/redntp_ZZTo2L2Nu_TuneZ2star_8TeV_pythia6_tauola.root");

    // gluglu higgs samples                                                                                              
    mc_2012[26] = TFile::Open(redntpDir+"/redntp.52xv5."+preselectionLevel+".paper-Hgg-scale.fixPerEleMva/merged/redntp_GluGluToHToGG_M-125_8TeV-powheg-pythia6.root");
    // vbf higgs samples                                                                                                 
    mc_2012[27] = TFile::Open(redntpDir+"/redntp.52xv5."+preselectionLevel+".paper-Hgg-scale.fixPerEleMva/merged/redntp_VBF_HToGG_M-125_8TeV-powheg-pythia6.root");
    // W/Z H higgs samples: only W->lnu                                                                                  
    mc_2012[28] = TFile::Open(redntpDir+"/redntp.52xv5."+preselectionLevel+".paper-Hgg-scale.fixPerEleMva/merged/redntp_WH_ZH_HToGG_M-125_8TeV-pythia6.root");
    // W/Z H higgs samples: only W->qq
    mc_2012[29] = TFile::Open(redntpDir+"/redntp.52xv5."+preselectionLevel+".paper-Hgg-scale.fixPerEleMva/merged/redntp_WH_ZH_HToGG_M-125_8TeV-pythia6.root");
    // W/Z H higgs samples: only Z->ll 
    mc_2012[30] = TFile::Open(redntpDir+"/redntp.52xv5."+preselectionLevel+".paper-Hgg-scale.fixPerEleMva/merged/redntp_WH_ZH_HToGG_M-125_8TeV-pythia6.root");
    // W/Z H higgs samples: only Z->qq
    mc_2012[31] = TFile::Open(redntpDir+"/redntp.52xv5."+preselectionLevel+".paper-Hgg-scale.fixPerEleMva/merged/redntp_WH_ZH_HToGG_M-125_8TeV-pythia6.root");
    // W/Z H higgs samples: only Z->nunu
    mc_2012[32] = TFile::Open(redntpDir+"/redntp.52xv5."+preselectionLevel+".paper-Hgg-scale.fixPerEleMva/merged/redntp_WH_ZH_HToGG_M-125_8TeV-pythia6.root");
  }
  
  // signal mc for different masses
  for (int i=0;i<11;i++){  
    char hmass[100]; sprintf(hmass,"%d",h_masses[i]);

    // gluglu samples                                                                                                          
    mc_gluglu_2012[i] = TFile::Open(redntpDir+"/redntp.52xv5."+preselectionLevel+".paper-Hgg-scale.fixPerEleMva/merged/redntp_GluGluToHToGG_M-"+ hmass + "_8TeV-powheg-pythia6.root");

    // vbf higgs samples                 
    mc_vbf_2012[i] = TFile::Open(redntpDir+"/redntp.52xv5."+preselectionLevel+".paper-Hgg-scale.fixPerEleMva/merged/redntp_VBF_HToGG_M-"+ hmass +"_8TeV-powheg-pythia6.root");

    // W/Z H higgs samples  - all together
    if (i!=6 && i!=7) mc_wzh_2012[i] = TFile::Open(redntpDir+"/redntp.52xv5."+preselectionLevel+".paper-Hgg-scale.fixPerEleMva/merged/redntp_WH_ZH_HToGG_M-"+ hmass +"_8TeV-pythia6.root");
    if (i==6 || i==7) mc_wzh_2012[i] = TFile::Open(redntpDir+"/redntp.52xv5."+preselectionLevel+".paper-Hgg-scale.fixPerEleMva/merged/redntp_WH_ZH_HToGG_M-"+ hmass +"_8TeV-pythia6.root");
  }

  // k factors - same scale factors for 7 TeV and 8 TeV
  double kfactordiphot = 1.3;                // for prompt and box 
  double kfactordiphotmadgraph = 1.15;       // for madgraph di-jets
  double kfactorgamjet = 1.3;
  double kfactorqcd = 1;
  double kfactordy = 1.15;
  double kfactorwg = 1;       
  double kfactorzg = 1;       
  double kfactorwmgg = 1;     
  double kfactorwpgg = 1;     
  double kfactorzgg = 1;      
  double kfactorttgg = 1;     
  double kfactorttjets = 1;   
  double kfactorwjets = 1;    
  double kfactorww = 1;       
  double kfactorwz = 1;       
  double kfactorzz = 1;       

  // cross sections at 8 TeV and scaling
  double boosthiggs(1);
  double cross_mc[33];

  cross_mc[0]  = 424.8 * kfactordiphot;                 // box
  cross_mc[1]  = 15.54 * kfactordiphot;                 // box
  cross_mc[2]  = 0.029038 * kfactordiphot;              // box
  // cross_mc[3]  = 264.5 * kfactordiphot;              // born
  // cross_mc[4]  = 25.48 * kfactordiphot;              // born
  // cross_mc[5]  = 0.029038 * kfactordiphot;           // born
  cross_mc[3]  = 81. * kfactordiphotmadgraph;           // madgraph
  cross_mc[4]  = 0. * kfactordiphotmadgraph;            // madgraph; chiara: per non riscrivere tutto metto 3 volte lo stesso sample e sigma a zero
  cross_mc[5]  = 0. * kfactordiphotmadgraph;            // madgraph; chiara: per non riscrivere tutto metto 3 volte lo stesso sample e sigma a zero
  cross_mc[6]  = 0.001835 * 81930.0 * kfactorgamjet;    // gjet
  cross_mc[7]  = 0.05387 * 8884.0 * kfactorgamjet;      // gjet
  cross_mc[8]  = 0.000235 * 5.195e+07 * kfactorqcd;     // qcd pt>40
  cross_mc[9]  = 0.002175 * 2.365e+07 * kfactorqcd;     // qcd 30<pt<40 
  cross_mc[10] = 2950. * kfactordy;                     // drell yan          
  cross_mc[11] = 322.356 * kfactorwg;       // Wgamma, W->lnu - da AN FP    (461.6 secondo prep) 
  cross_mc[12] = 0. * kfactorwg;            // Wgamma, dummy
  cross_mc[13] = 0. * kfactorwg;            // Wgamma, dummy
  cross_mc[14] = 181.338 * kfactorzg;       // Zgamma, inclusive - da AN FP (132.6 secondo prep)
  cross_mc[15] = 0. * kfactorzg;            // Zgamma, dummy
  cross_mc[16] = 0. * kfactorzg;            // Zgamma, dummy
  cross_mc[17] = 0.0504 * kfactorwmgg;      // W-gg  - da AN FP
  cross_mc[18] = 0.0667 * kfactorwpgg;      // W+gg  - da AN FP
  cross_mc[19] = 0.068 * kfactorzgg;        // Z+gg  - da AN FP
  cross_mc[20] = 0.001316 * kfactorttgg;    // tt+gg - da AN FP
  cross_mc[21] = 225.197 * kfactorttjets;   // tt+jets; da twiki https://twiki.cern.ch/twiki/bin/viewauth/CMS/StandardModelCrossSectionsat8TeV?skin=drupal   
  cross_mc[22] = 37509. * kfactorwjets;     // W+jets; da spreadsheet di HWW
  cross_mc[23] = 5.99515 * kfactorww;       // WW; da spreadsheet di HWW
  cross_mc[24] = 0.7346 * kfactorwz;        // WZ; da spreadsheet di HWW
  cross_mc[25] = 0.3649 * kfactorzz;        // ZZ; da spreadsheet di HWW

  /*
  cross_mc[0]  = 0.;           // box
  cross_mc[1]  = 0.;           // box
  cross_mc[2]  = 0.;           // box
  cross_mc[3]  = 0.;           // madgraph
  cross_mc[4]  = 0.;           // madgraph
  cross_mc[5]  = 0.;           // madgraph
  cross_mc[6]  = 0.;           // gjet
  cross_mc[7]  = 0.;           // gjet
  cross_mc[8]  = 0.;           // qcd pt>40
  cross_mc[9]  = 0.;           // qcd 30<pt<40 
  cross_mc[10] = 0.;           // drell yan          
  cross_mc[11] = 0.;           // Wgamma, W->lnu - da AN FP    (461.6 secondo prep) 
  cross_mc[12] = 0.;           // Wgamma, dummy
  cross_mc[13] = 0.;           // Wgamma, dummy
  cross_mc[14] = 0.;           // Zgamma, inclusive - da AN FP (132.6 secondo prep)
  cross_mc[15] = 0. * kfactorzg;            // Zgamma, dummy
  cross_mc[16] = 0. * kfactorzg;            // Zgamma, dummy
  cross_mc[17] = 0.0504 * kfactorwmgg;      // W-gg  - da AN FP
  cross_mc[18] = 0.0667 * kfactorwpgg;      // W+gg  - da AN FP
  cross_mc[19] = 0.068 * kfactorzgg;        // Z+gg  - da AN FP
  cross_mc[20] = 0.001316 * kfactorttgg;    // tt+gg - da AN FP
  cross_mc[21] = 0.;                        // tt+jets; 
  cross_mc[22] = 0.;                        // W+jets; da spreadsheet di HWW
  cross_mc[23] = 0. * kfactorww;       // WW; da spreadsheet di HWW
  cross_mc[24] = 0. * kfactorwz;        // WZ; da spreadsheet di HWW
  cross_mc[25] = 0. * kfactorzz;        // ZZ; da spreadsheet di HWW
  */


  // sigmaxBR SM @ 125 GeV
  cross_mc[26] = 19.52  * 2.13e-03 * boosthiggs;           // glu glu higgs SM 125,    
  cross_mc[27] = 1.578  * 2.13e-03 * boosthiggs;           // vbf higgs SM 125,        
  cross_mc[28] = 0.6966 * 2.13e-03 * boosthiggs * 0.3257;  // WH, W->lnu, higgs SM 125  
  cross_mc[29] = 0.6966 * 2.13e-03 * boosthiggs * 0.6743;  // WH, W->qq,  higgs SM 125  
  cross_mc[30] = 0.3943 * 2.13e-03 * boosthiggs * 0.10096; // ZH, Z->ll,  higgs SM 125  
  cross_mc[31] = 0.3943 * 2.13e-03 * boosthiggs * 0.6991;  // ZH, Z->qq,  higgs SM 125  
  cross_mc[32] = 0.3943 * 2.13e-03 * boosthiggs * 0.20;    // ZH, Z->nn,  higgs SM 125  
  // sigmaxBR FB   @ 125 GeV
  /*
  cross_mc[26]  = 0; // glu glu FP higgs
  cross_mc[27]  = 1.578  * 0.0154 * boosthiggs;             // vbf FP higgs, from YellowReport @ 8 TeV
  cross_mc[28]  = 0.6966 * 0.0154 * boosthiggs * 0.3257;    // WH, W->lnu, higgs SM 125   
  cross_mc[29]  = 0.6966 * 0.0154 * boosthiggs * 0.6743;    // WH, W->qq,  higgs SM 125   
  cross_mc[30]  = 0.3943 * 0.0154 * boosthiggs * 0.10096;   // ZH, Z->ll,  higgs SM 125   
  cross_mc[31]  = 0.3943 * 0.0154 * boosthiggs * 0.6991;    // ZH, Z->qq,  higgs SM 125   
  cross_mc[32]  = 0.3943 * 0.0154 * boosthiggs * 0.20;      // ZH, Z->nn,  higgs SM 125   
  */

  // ===================================
  // sigmaxBR SM per il segnale alle diverse masse 
  double cross_gluglu_mc[11], cross_vbf_mc[11], cross_wzh_mc[11];
  cross_gluglu_mc[0]  = 30.12 * 1.58e-03 * boosthiggs;    // 100    // from YellowReport @ 8 TeV
  cross_gluglu_mc[1]  = 27.39 * 1.77e-03 * boosthiggs;    // 105    // " " 
  cross_gluglu_mc[2]  = 25.04 * 1.95e-03 * boosthiggs;    // 110    // " " 
  cross_gluglu_mc[3]  = 22.96 * 2.11e-03 * boosthiggs;    // 115    // " "
  cross_gluglu_mc[4]  = 21.13 * 2.23e-03 * boosthiggs;    // 120    // " "
  cross_gluglu_mc[5]  = 19.52 * 2.28e-03 * boosthiggs;    // 125    // " "
  cross_gluglu_mc[6]  = 18.07 * 2.25e-03 * boosthiggs;    // 130    // " "  
  cross_gluglu_mc[7]  = 16.79 * 2.13e-03 * boosthiggs;    // 135    // " "  
  cross_gluglu_mc[8]  = 15.63 * 1.93e-03 * boosthiggs;    // 140    // " "  
  cross_gluglu_mc[9]  = 14.59 * 1.68e-03 * boosthiggs;    // 145    // " " 
  cross_gluglu_mc[10] = 13.65 * 1.37e-03 * boosthiggs;    // 150    // " " 
  // 
  cross_vbf_mc[0]  = 1.988 * 1.58e-03 * boosthiggs;    // 100    // from YellowReport @ 8 TeV
  cross_vbf_mc[1]  = 1.897 * 1.77e-03 * boosthiggs;    // 105    // 
  cross_vbf_mc[2]  = 1.809 * 1.95e-03 * boosthiggs;    // 110    // 
  cross_vbf_mc[3]  = 1.729 * 2.11e-03 * boosthiggs;    // 115    // 
  cross_vbf_mc[4]  = 1.649 * 2.23e-03 * boosthiggs;    // 120    // 
  cross_vbf_mc[5]  = 1.578 * 2.28e-03 * boosthiggs;    // 125    // 
  cross_vbf_mc[6]  = 1.511 * 2.25e-03 * boosthiggs;    // 130    // 
  cross_vbf_mc[7]  = 1.448 * 2.13e-03 * boosthiggs;    // 135    // 
  cross_vbf_mc[8]  = 1.389 * 1.93e-03 * boosthiggs;    // 140    // 
  cross_vbf_mc[9]  = 1.333 * 1.68e-03 * boosthiggs;    // 145    // 
  cross_vbf_mc[10] = 1.280 * 1.37e-03 * boosthiggs;    // 150    // 
  // 
  cross_wzh_mc[0]  = (1.432 + 0.7807)  * 1.58e-03 * boosthiggs;    // 100    // from YellowReport @ 8TeV
  cross_wzh_mc[1]  = (1.229 + 0.6750)  * 1.77e-03 * boosthiggs;    // 105    // ""
  cross_wzh_mc[2]  = (1.060 + 0.5869)  * 1.95e-03 * boosthiggs;    // 110    // ""
  cross_wzh_mc[3]  = (0.9165 + 0.5117) * 2.11e-03 * boosthiggs;    // 115    // ""
  cross_wzh_mc[4]  = (0.7966 + 0.4483) * 2.23e-03 * boosthiggs;    // 120    // ""
  cross_wzh_mc[5]  = (0.6966 + 0.3943) * 2.28e-03 * boosthiggs;    // 125    // ""
  cross_wzh_mc[6]  = (0.6095 + 0.3473) * 2.25e-03 * boosthiggs;    // 130    // ""
  cross_wzh_mc[7]  = (0.5351 + 0.3074) * 2.13e-03 * boosthiggs;    // 135    // ""
  cross_wzh_mc[8]  = (0.4713 + 0.2728) * 1.93e-03 * boosthiggs;    // 140    // ""
  cross_wzh_mc[9]  = (0.4164 + 0.2424) * 1.68e-03 * boosthiggs;    // 145    // ""
  cross_wzh_mc[10] = (0.3681 + 0.2159) * 1.37e-03 * boosthiggs;    // 150    // ""
  //
  //
  /*
  // sigmaxBR FP per il segnale alle diverse masse     
  cross_gluglu_mc[0]  = 0.;    // 100    // from YellowReport @ 8 TeV
  cross_gluglu_mc[1]  = 0.;
  cross_gluglu_mc[2]  = 0.;
  cross_gluglu_mc[3]  = 0.;
  cross_gluglu_mc[4]  = 0.;
  cross_gluglu_mc[5]  = 0.;
  cross_gluglu_mc[6]  = 0.;
  cross_gluglu_mc[7]  = 0.;
  cross_gluglu_mc[8]  = 0.;
  cross_gluglu_mc[9]  = 0.;
  cross_gluglu_mc[10] = 0.;
  // 
  cross_vbf_mc[0]  = 1.988 * 0.1824 * boosthiggs;    // 100    // from YellowReport @ 8 TeV
  cross_vbf_mc[1]  = 1.897 * 0.1029 * boosthiggs;    // 105    // 
  cross_vbf_mc[2]  = 1.809 * 5.95e-02 * boosthiggs;    // 110    // 
  cross_vbf_mc[3]  = 1.729 * 3.61e-02 * boosthiggs;    // 115    // 
  cross_vbf_mc[4]  = 1.649 * 2.31e-02 * boosthiggs;    // 120    // 
  cross_vbf_mc[5]  = 1.578 * 1.54e-02 * boosthiggs;    // 125    // 
  cross_vbf_mc[6]  = 1.511 * 1.06e-02 * boosthiggs;    // 130    // 
  cross_vbf_mc[7]  = 1.448 * 7.50e-03 * boosthiggs;    // 135    // 
  cross_vbf_mc[8]  = 1.389 * 5.38e-03 * boosthiggs;    // 140    // 
  cross_vbf_mc[9]  = 1.333 * 3.86e-03 * boosthiggs;    // 145    // 
  cross_vbf_mc[10] = 1.280 * 2.70e-03 * boosthiggs;    // 150    // 
  //
  cross_wzh_mc[0]  = (1.432 + 0.7807)  * 0.1824 * boosthiggs;      // 100    // from YellowReport @ 8TeV
  cross_wzh_mc[1]  = (1.229 + 0.6750)  * 0.1029 * boosthiggs;      // 105    // ""
  cross_wzh_mc[2]  = (1.060 + 0.5869)  * 5.95e-02 * boosthiggs;    // 110    // ""
  cross_wzh_mc[3]  = (0.9165 + 0.5117) * 3.61e-02 * boosthiggs;    // 115    // ""
  cross_wzh_mc[4]  = (0.7966 + 0.4483) * 2.31e-02 * boosthiggs;    // 120    // ""
  cross_wzh_mc[5]  = (0.6966 + 0.3943) * 1.54e-02 * boosthiggs;    // 125    // ""
  cross_wzh_mc[6]  = (0.6095 + 0.3473) * 1.06e-02 * boosthiggs;    // 130    // ""
  cross_wzh_mc[7]  = (0.5351 + 0.3074) * 7.50e-03 * boosthiggs;    // 135    // ""
  cross_wzh_mc[8]  = (0.4713 + 0.2728) * 5.38e-03 * boosthiggs;    // 140    // ""
  cross_wzh_mc[9]  = (0.4164 + 0.2424) * 3.86e-03 * boosthiggs;    // 145    // ""
  cross_wzh_mc[10] = (0.3681 + 0.2159) * 2.70e-03 * boosthiggs;    // 150    // ""
  */

  // getting the number of original events in each sample (processed with CMSSW)
  int n_mc_2012[33],n_gluglu_2012[11],n_vbf_2012[11],n_wzh_2012[11];  
  for(int i=0; i<28; i++){
    n_mc_2012[i] = 0;
    if(int_exp_2012>0) n_mc_2012[i] = ((TH1D*)mc_2012[i]->Get("ptphotgen1"))->GetEntries();  
  }
  if(int_exp_2012>0) { 
    n_mc_2012[28]  = ((TH1D*)mc_2012[28]->Get("ptphotgen1wl"))->GetEntries();  
    n_mc_2012[29]  = ((TH1D*)mc_2012[29]->Get("ptphotgen1wh"))->GetEntries();  
    n_mc_2012[30]  = ((TH1D*)mc_2012[30]->Get("ptphotgen1zl"))->GetEntries();  
    n_mc_2012[31]  = ((TH1D*)mc_2012[31]->Get("ptphotgen1zh"))->GetEntries();  
    n_mc_2012[32]  = ((TH1D*)mc_2012[32]->Get("ptphotgen1zn"))->GetEntries();  
  }

  for(int i=0; i<11; i++){ 
    n_gluglu_2012[i] = n_vbf_2012[i] = n_wzh_2012[i] = 0;
    n_gluglu_2012[i] = ((TH1D*)mc_gluglu_2012[i]->Get("ptphotgen1"))->GetEntries();
    n_vbf_2012[i]    = ((TH1D*)mc_vbf_2012[i]->Get("ptphotgen1"))->GetEntries();
    n_wzh_2012[i]    = ((TH1D*)mc_wzh_2012[i]->Get("ptphotgen1"))->GetEntries();
  }

  // setting the scaling factor to actual lumi 
  double scale_mc_2012[33];
  for(int i=0; i<33; i++){
    scale_mc_2012[i] = 0; 
    if(int_exp_2012>0) { 
      scale_mc_2012[i] = cross_mc[i] * int_exp_2012 / n_mc_2012[i];
    }
  }

  // setting the scaling factor to actual lumi for the different signal masses 
  double scale_gluglu_mc_2012[11], scale_vbf_mc_2012[11], scale_wzh_mc_2012[11];  
  for(int i=0; i<11; i++){
    scale_gluglu_mc_2012[i] = 0; 
    scale_vbf_mc_2012[i]    = 0; 
    scale_wzh_mc_2012[i]    = 0; 
    if(int_exp_2012>0) { 
      scale_gluglu_mc_2012[i] = cross_gluglu_mc[i] * int_exp_2012 / n_gluglu_2012[i];
      scale_vbf_mc_2012[i]    = cross_vbf_mc[i] * int_exp_2012 / n_vbf_2012[i];
      scale_wzh_mc_2012[i]    = cross_wzh_mc[i] * int_exp_2012 / n_wzh_2012[i];
    }
  }

  // char for output name
  char name[1000];
  char allcut[3000];
  sprintf(allcut,"%3.1f%s%3.1f%s%3.1f%s%3.1f%s%3.1f%s%3.1f%s%3.1f%s%3.1f%s%3.1f%s%3.1f%s%3.1f%s%3.1f%s%3.1f%s%3.1f%s%3.1f%s%3.1f%s%3.1f%s%d%s%d%s%d%s%d%s%d",pt1,"_",pt2,"_",pthiggsmin,"_",pthiggsmax,"_",ptj1,"_",ptj2,"_",misset,"_",deltae,"_",zep,"_",mjj,"_",deltap,"_",jetmet,"_",p1met,"_",p2met,"_",hmet,"_",phigg,"_",eb,"_",r9,"_",thirdcat,"_",leptontag,"_",leptonveto,"_",cic);
  
  // output root file
  sprintf(name,"%s%s%s%s%s","results_gg/histo_",variableData.c_str(),"_",allcut,".root");
  TFile * hOutputFile = new TFile(name, "RECREATE" ) ;

  // histograms needed by the machinery
  TH1D* vardata   = new TH1D("vardata",  "vardata",  nbin,min,max);
  // TH1D* vardatacs = new TH1D("vardatacs","vardatacs",nbin,min,max);
  TH1D* var_mc_2012[33];            // 33 = MC samples we have 
  TH1D* var_gluglu_2012[11];        // 11 = number of masses
  TH1D* var_vbf_2012[11];
  TH1D* var_wzh_2012[11];
  TH1D* var[14];                   // 14 = number of species

  for (int i=0; i<14; i++) {
    sprintf(name,"%s%d","var",i);
    var[i] = new TH1D(name,name,nbin,min,max);
  }

  for (int i=0; i<33; i++) {
    sprintf(name,"%s%d","var_mc_2012_",i);
    var_mc_2012[i] = new TH1D(name,name,nbin,min,max);
  }

  for (int i=0; i<11; i++) { 
    sprintf(name,"%s%d","var_gluglu_2012_",h_masses[i]);
    var_gluglu_2012[i] = new TH1D(name,name,nbin,min,max);
    sprintf(name,"%s%d","var_vbf_2012_",h_masses[i]);
    var_vbf_2012[i] = new TH1D(name,name,nbin,min,max);
    sprintf(name,"%s%d","var_wzh_2012_",h_masses[i]);
    var_wzh_2012[i] = new TH1D(name,name,nbin,min,max);
  }

  // creating the fillers and setting cuts
  fillPlot2012 data_fill((TTree*)data->Get("AnaTree"), 1);
  // fillPlot2012 datacs_fill((TTree*)datacs->Get("AnaTree"), 1);
  data_fill.Setcuts(pt1,pt2,pthiggsmin,pthiggsmax,ptj1,ptj2,misset,deltae,zep,mjj,deltap,jetmet,p1met,p2met,hmet,phigg,eb,r9,thirdcat,leptontag,leptonveto,vbfveto,inclusive);
  // datacs_fill.Setcuts(pt1,pt2,pthiggsmin,pthiggsmax,ptj1,ptj2,misset,deltae,zep,mjj,deltap,jetmet,p1met,p2met,hmet,phigg,eb,r9,thirdcat,leptontag,leptonveto,vbfveto,inclusive);
  if(cic>0) {
    data_fill.setCic(cic);
    // datacs_fill.setCic(cic);
    data_fill.usePFCic(usepfcic);  
    // datacs_fill.usePFCic(usepfcic);
  }
  
  fillPlot2012* mc_2012_fill[33];
  fillPlot2012* mc_gluglu_2012_fill[11];    
  fillPlot2012* mc_vbf_2012_fill[11];  
  fillPlot2012* mc_wzh_2012_fill[11];  

  for (int i=0; i<33; i++){
    if(int_exp_2012>0) mc_2012_fill[i] = new fillPlot2012((TTree*)mc_2012[i]->Get("AnaTree"), 1);   
    if(int_exp_2012>0) mc_2012_fill[i]->Setcuts(pt1,pt2,pthiggsmin,pthiggsmax,ptj1,ptj2,misset,deltae,zep,mjj,deltap,jetmet,p1met,p2met,hmet,phigg,eb,r9,thirdcat,leptontag,leptonveto,vbfveto,inclusive);
    
    mc_2012_fill[i]->DoPuReweight();
    // mc_2012_fill[i]->DoPtReweight();     
    // if(i==3 || i==4 || i==5) mc_2012_fill[i]->DoRemoveDiphot();
    
    if(cic>0){
      if(int_exp_2012>0) { 
	mc_2012_fill[i]->setCic(cic);
	mc_2012_fill[i]->usePFCic(usepfcic);
      }
    }
  }

  for (int i=0; i<11; i++){  
    mc_gluglu_2012_fill[i] = new fillPlot2012((TTree*)mc_gluglu_2012[i]->Get("AnaTree"), 1);
    mc_vbf_2012_fill[i]    = new fillPlot2012((TTree*)mc_vbf_2012[i]->Get("AnaTree"), 1);
    mc_wzh_2012_fill[i]    = new fillPlot2012((TTree*)mc_wzh_2012[i]->Get("AnaTree"), 1);
    //
    mc_gluglu_2012_fill[i]->Setcuts(pt1,pt2,pthiggsmin,pthiggsmax,ptj1,ptj2,misset,deltae,zep,mjj,deltap,jetmet,p1met,p2met,hmet,phigg,eb,r9,thirdcat,leptontag,leptonveto,vbfveto,inclusive);
    mc_vbf_2012_fill[i]->Setcuts(pt1,pt2,pthiggsmin,pthiggsmax,ptj1,ptj2,misset,deltae,zep,mjj,deltap,jetmet,p1met,p2met,hmet,phigg,eb,r9,thirdcat,leptontag,leptonveto,vbfveto,inclusive);
    mc_wzh_2012_fill[i]->Setcuts(pt1,pt2,pthiggsmin,pthiggsmax,ptj1,ptj2,misset,deltae,zep,mjj,deltap,jetmet,p1met,p2met,hmet,phigg,eb,r9,thirdcat,leptontag,leptonveto,vbfveto,inclusive);
    //
    mc_gluglu_2012_fill[i]->DoPuReweight();
    mc_vbf_2012_fill[i]->DoPuReweight();
    mc_wzh_2012_fill[i]->DoPuReweight();
    //
    // mc_gluglu_2012_fill[i]->DoPtReweight();    
    // mc_vbf_2012_fill[i]->DoPtReweight();       
    // mc_wzh_2012_fill[i]->DoPtReweight();       
    //
    mc_gluglu_2012_fill[i]->setCic(cic);
    mc_vbf_2012_fill[i]->setCic(cic);
    mc_wzh_2012_fill[i]->setCic(cic);
    //
    mc_gluglu_2012_fill[i]->usePFCic(usepfcic);
    mc_vbf_2012_fill[i]->usePFCic(usepfcic);
    mc_wzh_2012_fill[i]->usePFCic(usepfcic);
  }

  // smear mc
  // for (int i=0; i<33; i++){
  //    if(int_exp_2012>0) mc_2012_fill[i]->DoSmearing(1.,0.0001);
  //  }  


  // filling histograms
  std::cout << " ++++++++++++++ DATA ++++++++++++++++" << std::endl;
  cout << "running over " << ((TTree*)data->Get("AnaTree"))->GetEntries("") << " data events" <<  endl;
  if (variableData == "massgg") {  
    sprintf(name,"%s%s%s","results_gg/events_",allcut,".txt");
    data_fill.Writetxt(name);
    sprintf(name,"%s%s%s","results_gg/events_",allcut,".root");
    data_fill.WriteRoot(name);
  }
  vardata->Add(data_fill.Plot(variableData,"data", nbin, min, max, 0, 100)); 
  std::cout << "Selected events on data " << vardata->GetEntries() << std::endl;

  // cout << "running over " << ((TTree*)datacs->Get("AnaTree"))->GetEntries("") << " data events (for cs)" <<  endl; 
  //if (variableData == "massgg") {  
  //  sprintf(name,"%s%s%s","results_gg/events_",allcut,"_cs.txt");
  //  datacs_fill.Writetxt(name);
  //  sprintf(name,"%s%s%s","results_gg/events_",allcut,"_cs.root");
  //  datacs_fill.WriteRoot(name);
  // }
  //vardatacs->Add(datacs_fill.Plot(variableData,"datacs", nbin, min, max, 1, 110)); 
  //std::cout << "Selected events on data cs " << vardatacs->GetEntries() << std::endl;

  std::cout << " ++++++++++++++ MC ++++++++++++++++" << std::endl;
  for (int i=0; i<33; i++){ 
    sprintf(name,"%s%s",mcnames[i].c_str()," 2012");
    if(int_exp_2012>0) {
      cout << "running over " << ((TTree*)mc_2012[i]->Get("AnaTree"))->GetEntries("") << " " << name << " events" <<  endl;   
      if (i>=28 && i<33) { 
	TH1F *myRecordH = new TH1F("myRecordH","",100,0.,100000000.);
	if (i==28) ((TTree*)mc_2012[i]->Get("AnaTree"))->Project("myRecordH","gen_custom_processId","gen_custom_processId==3101 || gen_custom_processId==3201 || gen_custom_processId==3301");
	if (i==29) ((TTree*)mc_2012[i]->Get("AnaTree"))->Project("myRecordH","gen_custom_processId","gen_custom_processId==3401");
	if (i==30) ((TTree*)mc_2012[i]->Get("AnaTree"))->Project("myRecordH","gen_custom_processId","gen_custom_processId==4101 || gen_custom_processId==4201 || gen_custom_processId==4301");
	if (i==31) ((TTree*)mc_2012[i]->Get("AnaTree"))->Project("myRecordH","gen_custom_processId","gen_custom_processId==4401");
	if (i==32) ((TTree*)mc_2012[i]->Get("AnaTree"))->Project("myRecordH","gen_custom_processId","gen_custom_processId==4501");
	cout << "of them " << myRecordH->GetEntries() << " are from the interesting production mechanism" << endl;
	delete myRecordH;
      } 
      sprintf(name,"%s%s%s%s%s","results_gg/events_",mcnames[i].c_str(),"_2012_",allcut,".root");
      if (variableMC == "massgg") mc_2012_fill[i]->WriteRoot(name);             
      var_mc_2012[i]->Add( mc_2012_fill[i]->Plot(variableMC, name, nbin, min, max, 0, i) );  
      std::cout << "Selected events on mc2012 " << name << " " << var_mc_2012[i]->GetEntries() << std::endl;
    }
  }

  std::cout << " ++++++++++++++ signal MC ++++++++++++++++" << std::endl;
  /*
  for (int i=0; i<11; i++){          
    cout << "running over " << ((TTree*)mc_gluglu_2012[i]->Get("AnaTree"))->GetEntries("") << " gluglu M=" << h_masses[i] << " events" <<  endl; 
    sprintf(name,"%s%d%s%s%s","results_gg/events_gluglu",h_masses[i],"_2012_",allcut,".root");
    var_gluglu_2012[i]->Add(mc_gluglu_2012_fill[i]->Plot(variableMC, name, nbin, min, max, 0, 50));
    std::cout << "Selected events on mc2012 gluglu " << h_masses[i] << " " << var_gluglu_2012[i]->GetEntries() << std::endl;
 
    cout << "running over " << ((TTree*)mc_vbf_2012[i]->Get("AnaTree"))->GetEntries("") << " vbf M=" << h_masses[i] << " events" <<  endl; 
    sprintf(name,"%s%d%s%s%s","results_gg/events_vbf",h_masses[i],"_2012_",allcut,".root");
    var_vbf_2012[i]->Add(mc_vbf_2012_fill[i]->Plot(variableMC, name, nbin, min, max, 0, 50));
    std::cout << "Selected events on mc2012 vbf " << h_masses[i] << " " << var_vbf_2012[i]->GetEntries() << std::endl;

    cout << "running over " << ((TTree*)mc_wzh_2012[i]->Get("AnaTree"))->GetEntries("") << " wzh M=" << h_masses[i] << " events" <<  endl; 
    sprintf(name,"%s%d%s%s%s","results_gg/events_wzh",h_masses[i],"_2012_",allcut,".root");  
    var_wzh_2012[i]->Add(mc_wzh_2012_fill[i]->Plot(variableMC, name, nbin, min, max, 0, 50));
    std::cout << "Selected events on mc2012 wzh " << h_masses[i] << " " << var_wzh_2012[i]->GetEntries() << std::endl;
  }

  // scale mc to equivalent lumi
  for (int i=0; i<33; i++){
    if(int_exp_2012>0) var_mc_2012[i]->Scale(scale_mc_2012[i]);  
  }

  // scale signals to equivalent lumi
  for (int i=0; i<11; i++){ 
    if(int_exp_2012>0) {
      var_gluglu_2012[i] -> Scale(scale_gluglu_mc_2012[i]);  
      var_vbf_2012[i]    -> Scale(scale_vbf_mc_2012[i]);  
      var_wzh_2012[i]    -> Scale(scale_wzh_mc_2012[i]);  
    }
  }
  */

  // counting number of events passing selection (scaled)   
  double num_mc_2012[33],num_gluglu_2012[11],num_vbf_2012[11],num_wzh_2012[11];
  double num_uns_mc_2012[33],num_uns_gluglu_2012[11],num_uns_vbf_2012[11],num_uns_wzh_2012[11];

  for (int i=0; i<33; i++){ 
    num_mc_2012[i] = 0;
    if(int_exp_2012>0) num_mc_2012[i] = var_mc_2012[i]->Integral();        
    // 
    num_uns_mc_2012[i] = 0;
    if(int_exp_2012>0) num_uns_mc_2012[i] = var_mc_2012[i]->GetEntries();  
  }

  for(int i=0; i<11; i++){
    num_gluglu_2012[i] = num_vbf_2012[i] = num_wzh_2012[i] = 0;
    num_gluglu_2012[i] = var_gluglu_2012[i]->Integral(); 
    num_vbf_2012[i]    = var_vbf_2012[i]->Integral(); 
    num_wzh_2012[i]    = var_wzh_2012[i]->Integral(); 
    // 
    num_uns_gluglu_2012[i] = num_uns_vbf_2012[i] = num_uns_wzh_2012[i] = 0;
    num_uns_gluglu_2012[i] = var_gluglu_2012[i]->GetEntries(); 
    num_uns_vbf_2012[i]    = var_vbf_2012[i]->GetEntries(); 
    num_uns_wzh_2012[i]    = var_wzh_2012[i]->GetEntries(); 
  }

  // scale control sample
  vardata->Sumw2();
  //vardatacs->Sumw2();   
  double num_data    = vardata->Integral();
  //double num_data_cs = vardatacs->Integral();  
  //vardatacs->Scale(num_data/num_data_cs); 

  // stack histograms - chiara: ancora da controllare
  for (int i=1; i<nbin+1; i++){      
    for (int j=0; j<14; j++){            
      int offset(0);
      if(j>0) offset = 6;    // to add higgs contributions up 
      if(j>1) offset = 8;    // to add dibosons contribuions up
      if(j>6) offset = 9;    // to add Wgg contributions up  
      if(j>7) offset = 11;   // to add Zg contributions up   
      if(j>8) offset = 13;   // to add Wg contributions up   
      if(j>10) offset =14;   // to add QCD contributions up  
      if(j>11) offset =15;   // to add gamma+jets contributions up
      if(j>12) offset =17;   // to add born contributions up
      for (int k=0 ; k<33-j-offset; k++){   
	if(int_exp_2012>0) var[j]->SetBinContent(i,var_mc_2012[k]->GetBinContent(i) + var[j]->GetBinContent(i));
      }	
    }    
  }

  // final plots
  char ytitle[100];
  sprintf(ytitle,"%s%d%s","N_{ev}/",int(int_exp_2012),"pb^{-1}");
  for(int i=0; i<14; i++){
    var[i]->SetTitle("");
    var[i]->SetStats(0);
    var[i]->SetTitleOffset(1.25,"Y");
    var[i]->SetYTitle(ytitle);
    var[i]->SetXTitle(axis.c_str());
    var[i]->SetLineColor(kBlack);
    var[i]->SetLineWidth(2);
  }
  
  // legenda
  TLegendEntry *legge;
  TLegend *leg;
  leg = new TLegend(0.6,0.6,0.85,0.85);
  leg->SetFillStyle(0); leg->SetBorderSize(0); leg->SetTextSize(0.05);
  leg->SetFillColor(0);
  sprintf(name,"H M(125)");
  if(boosthiggs!=1) sprintf(name,"%s%2.0f", "H M(125) x ", boosthiggs);

  legge = leg->AddEntry(var[0], name, "f");
  legge = leg->AddEntry(var[1], "VV", "f");
  legge = leg->AddEntry(var[2], "Wjets", "f");
  legge = leg->AddEntry(var[3], "ttjets", "f");
  legge = leg->AddEntry(var[4], "ttgg", "f");
  legge = leg->AddEntry(var[5], "Zgg", "f");
  legge = leg->AddEntry(var[6], "Wgg", "f");
  legge = leg->AddEntry(var[7], "Zgamma", "f");
  legge = leg->AddEntry(var[8], "Wgamma", "f");  
  legge = leg->AddEntry(var[9], "DY", "f");
  legge = leg->AddEntry(var[10], "QCD", "f");
  legge = leg->AddEntry(var[11], "#gamma + jets", "f");
  legge = leg->AddEntry(var[12], "di-#gamma + jets", "f");
  legge = leg->AddEntry(var[13], "di-#gamma box", "f");
  
  // mc only plot: data vs background
  var[0]->SetFillColor(kYellow);
  var[0]->Draw();
  var[1]->SetFillColor(6);
  var[1]->Draw("same");
  var[2]->SetFillColor(2);
  var[2]->Draw("same");
  var[3]->SetFillColor(38);
  var[3]->Draw("same");
  var[4]->SetFillColor(46);
  var[4]->Draw("same");
  var[5]->SetFillColor(16);
  var[5]->Draw("same");
  var[6]->SetFillColor(4);
  var[6]->Draw("same");
  var[7]->SetFillColor(3);   
  var[7]->Draw("same");
  var[8]->SetFillColor(64);  
  var[8]->Draw("same");
  var[9]->SetFillColor(40);
  var[9]->Draw("same");
  var[10]->SetFillColor(95);
  var[10]->Draw("same");
  var[11]->SetFillColor(48);
  var[11]->Draw("same");
  var[12]->SetFillColor(52);
  var[12]->Draw("same");
  var[13]->SetFillColor(30);
  var[13]->Draw("same");

  leg->Draw();
  sprintf(name,"%s%s%s%s%s","results_gg/mc_",variableData.c_str(),"_",allcut,".png");
  c0->SaveAs(name);
  sprintf(name,"%s%s%s%s%s","results_gg/mc_",variableData.c_str(),"_",allcut,".root");
  c0->SaveAs(name);

  // higgs only plot 
  TH1D* myClone = (TH1D*) var_mc_2012[26]->Clone("myClone");  
  myClone->Add(var_mc_2012[27]);   
  myClone->Add(var_mc_2012[28]);   
  myClone->Add(var_mc_2012[29]);   
  myClone->Add(var_mc_2012[30]);   
  myClone->Add(var_mc_2012[31]);   
  myClone->Add(var_mc_2012[32]);   
  sprintf(ytitle,"%s%d%s","N_{ev}/",int(int_exp_2012),"pb^{-1}");
  myClone->SetXTitle(axis.c_str());
  myClone->SetYTitle(ytitle);
  myClone->SetTitle("");
  myClone->SetLineColor(kBlack);
  myClone->SetLineWidth(2);
  myClone->SetFillColor(kYellow);
  myClone->Draw();
  sprintf(name,"%s%s%s%s%s","results_gg/higgs_",variableData.c_str(),"_",allcut,".png");
  c0->SaveAs(name);

  // only VH
  TH1D* myClone2 = (TH1D*) var_mc_2012[28]->Clone("myClone");  
  myClone2->Add(var_mc_2012[29]);   
  myClone2->Add(var_mc_2012[30]);   
  myClone2->Add(var_mc_2012[31]);   
  myClone2->Add(var_mc_2012[32]);   

  // some calculation for s/sqrt(b)
  double integralhiggs ;
  double entrieshiggs ;
  double integralbkg ;
  if (variableData == "massgg") {
    sprintf(name,"%s%s%s%d%s","results_gg/optimalcut_",variableData.c_str(),"_",eb,".txt");
    integralhiggs = myClone->Integral(43,56);
    entrieshiggs  = myClone->Integral(0,201);
    integralbkg   = var[1]->Integral(30,71)/3.;  // chiara: perche' /3.?

    cout << endl;
    cout << "Analysis considering all signal samples" << endl;
    cout << "Analysis done in " << myClone->GetBinLowEdge(43) << " - " << myClone->GetBinLowEdge(56)+myClone->GetBinWidth(56) << " for signal" << endl;
    cout << "Analysis done in " << var[1]->GetBinLowEdge(30)  << " - " << var[1]->GetBinLowEdge(71)+var[1]->GetBinWidth(71)   << " for background" << endl;
    cout << "Fraction in signal box "  << integralhiggs/entrieshiggs << endl;
    cout << "Number of signal events " << integralhiggs/boosthiggs << endl;
    cout << "Number of bkg events "    << integralbkg << endl;
    cout << "S/sqrt(B) "   << integralhiggs/boosthiggs/sqrt(integralbkg) << endl;
    cout << "S/B "         << integralhiggs/boosthiggs/integralbkg << endl;
    cout << "probability " << probability_disc(integralhiggs/boosthiggs+integralbkg,integralbkg) << endl;
    cout << endl;
    cout << endl;
    cout << "chiaraaaaaaa: expected VH signal events in " 
	 << myClone2->GetBinLowEdge(4) << "-" << (myClone2->GetBinLowEdge(5)+myClone2->GetBinWidth(5)) 
	 << " = " << myClone2->Integral(4,5)/boosthiggs << endl;  
    cout << endl;

    // here only signal HW, W->lnu
    TH1D* myClone3 = (TH1D*) var_mc_2012[28]->Clone("myClone3");  
    integralhiggs = myClone3->Integral(43,56);
    entrieshiggs  = myClone3->Integral(0,201);
    cout << endl;
    cout << "Analysis considering only W->lnu" << endl;
    cout << "Number of signal events " << integralhiggs/boosthiggs << endl;
    cout << "S/sqrt(B) "   << integralhiggs/boosthiggs/sqrt(integralbkg) << endl;
    cout << "S/B "         << integralhiggs/boosthiggs/integralbkg << endl;
    cout << "probability " << probability_disc(integralhiggs/boosthiggs+integralbkg,integralbkg) << endl;

    // here only signal HZ, Z->nn
    TH1D* myClone4 = (TH1D*) var_mc_2012[32]->Clone("myClone4");  
    integralhiggs = myClone4->Integral(43,56);
    entrieshiggs  = myClone4->Integral(0,201);
    cout << endl;
    cout << "Analysis considering only Z->nn" << endl;
    cout << "Number of signal events " << integralhiggs/boosthiggs << endl;
    cout << "S/sqrt(B) "   << integralhiggs/boosthiggs/sqrt(integralbkg) << endl;
    cout << "S/B "         << integralhiggs/boosthiggs/integralbkg << endl;
    cout << "probability " << probability_disc(integralhiggs/boosthiggs+integralbkg,integralbkg) << endl;
  }  

  // data only plot
  vardata->SetXTitle(axis.c_str());
  vardata->SetTitle("");
  vardata->SetStats(0);
  vardata->SetMarkerStyle(8);
  vardata->SetMarkerSize(.9);
  vardata->SetTitleOffset(1.25,"Y");
  vardata->Draw("pe");
  sprintf(name,"%s%s%s%s%s","results_gg/data_",variableData.c_str(),"_",allcut,".png");
  c0->SaveAs(name);

  // data overlaid to mc
  legge = leg->AddEntry(vardata, "data", "p");

  double themax =   vardata->GetMaximum();
  if(var[0]->GetMaximum()>themax) themax = var[0]->GetMaximum();
  if (
      variableData == "etaphot1" || variableData == "etaphot2" ||
      variableData == "phiphot1" || variableData == "phiphot2" ||
      variableData == "etajet1" || variableData == "etajet2" ||
      variableData == "phijet1" || variableData == "phijet2"
      )
    vardata->SetMaximum(themax*2.0);
  else
    vardata->SetMaximum(themax*1.1);
  
  vardata->SetMinimum(0.000001);
  for(int ii=0; ii<14; ii++) var[ii]->SetMinimum(0.000001);
  var[0]->Draw("same");
  var[1]->Draw("same");
  var[2]->Draw("same");
  var[3]->Draw("same");
  var[4]->Draw("same");
  var[5]->Draw("same");
  var[6]->Draw("same");
  var[7]->Draw("same");
  var[8]->Draw("same");
  var[9]->Draw("same");
  var[10]->Draw("same");
  var[11]->Draw("same");
  var[12]->Draw("same");
  var[13]->Draw("same");
  leg->Draw();
  vardata->Draw("pesame");
  gPad->RedrawAxis();
  sprintf(name,"%s%s%s%s%s","results_gg/data-mc_",variableData.c_str(),"_",allcut,".root");
  c0->SaveAs(name);
  sprintf(name,"%s%s%s%s%s","results_gg/data-mc_",variableData.c_str(),"_",allcut,".png");
  c0->SaveAs(name);

  // data with control sample
  vardata->Draw("pe");
  //vardatacs->SetLineColor(46);
  //vardatacs->SetFillColor(42);
  //vardatacs->SetFillColor(7);
  //vardatacs->SetLineWidth(3);
  //vardatacs->Draw("hsame");

  vardata->Draw("pesame");
  // sprintf(name,"%s%s%s%s%s","results_gg/datacs_",variableData.c_str(),"_",allcut,".png");
  gPad->RedrawAxis();
  TLegendEntry *legge2;
  TLegend *leg2;
  leg2 = new TLegend(0.6,0.65,0.9,0.85);
  leg2->SetFillStyle(0); leg2->SetBorderSize(0); leg2->SetTextSize(0.035);
  leg2->SetFillColor(0);
  legge2 = leg2->AddEntry(vardata, "default sel.", "p");
  //legge2 = leg2->AddEntry(vardatacs, "control sample", "f");
  leg2->Draw();
  c0->SaveAs(name);

  // additional plot with larger binnin (only for invariant mass)
  if (variableData == "massgg"){
    //     var[0]->Rebin(4);
    //     themax = var[0]->GetMaximum();
    //     vardata->Rebin(4);
    //     cout << vardata->GetMaximum() << "   " << var[0]->GetMaximum() << endl;
    //     if(themax < vardata->GetMaximum()) themax = vardata->GetMaximum();
    //     var[0]->SetMaximum(themax*1.1);
    //     var[0]->SetMinimum(0.);
    //     var[0]->Draw();
    //     var[1]->Rebin(4);
    //     var[1]->Draw("same");
    //     var[2]->Rebin(4);
    //     var[2]->Draw("same");
    //     var[3]->Rebin(4);
    //     var[3]->Draw("same");
    //     var[4]->Rebin(4);
    //     var[4]->Draw("same");
    //     var[5]->Rebin(4);
    //     var[5]->Draw("same");
    //     leg->Draw();
    //     gPad->RedrawAxis();
    //     sprintf(name,"%s%s%s%s%s","results_gg/mc_rebin_",variableData.c_str(),"_",allcut,".png");
    //     c0->SaveAs(name);
    
    //     vardata->Draw("pesame");
    //     sprintf(name,"%s%s%s%s%s","results_gg/data-mc_rebin_",variableData.c_str(),"_",allcut,".png");
    //     c0->SaveAs(name);
    
    //     vardatacs->Rebin(4);
    //     vardata->SetMaximum(themax*1.1);
    //     vardata->Draw("pe");
    //     vardatacs->Draw("same");
    //     vardata->Draw("pesame");
    //     leg2->Draw();
    //     gPad->RedrawAxis();
    //     sprintf(name,"%s%s%s%s%s","results_gg/datacs_rebin_",variableData.c_str(),"_",allcut,".png");
    //     c0->SaveAs(name);
    
    //     double newmax(0);
    //     for (int i=5;i<50;i++){
    //       double tempmax = vardata->GetBinContent(i);
    //       if(tempmax > newmax) newmax = tempmax;
    //     }
    //     vardata->SetAxisRange(100.,150.);
    //     vardata->SetMaximum(newmax*1.3);
    //     vardata->Draw("pe");
    //     sprintf(name,"%s%s%s%s%s","results_gg/data_rebin_resize",variableData.c_str(),"_",allcut,".png");
    //     c0->SaveAs(name);

    sprintf(name,"%s%s%s","results_gg/yields_",allcut,".txt");
    ofstream outfile(name);  
    outfile << "####################################" << endl;
    outfile << "CUTS " << endl;
    outfile << "####################################" << endl;
    outfile << "ptphot1 : " << pt1 << endl;
    outfile << "ptphot2 : " << pt2 << endl;
    outfile << "ptjet1 : " << ptj1 << endl;
    outfile << "ptjet2 : " << ptj2 << endl;
    outfile << "met : " << misset  << endl;
    outfile << "deltaetacut : " << deltae << endl;
    outfile << "deltaphicut : " << deltap << endl;
    outfile << "deltaPhi(met,jet)cut : "   << jetmet << endl;
    outfile << "deltaPhi(met,phot1)cut : " << p1met << endl;
    outfile << "deltaPhi(met,phot2)cut : " << p2met << endl;
    outfile << "deltaPhi(met,higgs)cut : " << hmet << endl;
    outfile << "deltaPhi(phot1,phot2)cut : " << phigg << endl;
    outfile << "zeppencut : " << zep << endl;
    outfile << "invmassjetcut : " << mjj << endl;
    outfile << "CiC level : " << cic << endl;
    outfile << "Use PF cic : " << usepfcic << endl;
    outfile << "ebcat : " << eb << endl;
    outfile << "r9cat : " << r9 << endl;
    outfile << "thirdcat : " << thirdcat << endl;
    outfile << "leptontag : " << leptontag << endl;
    outfile << "leptonveto : " << leptonveto << endl;
    outfile << "vbfveto : " << vbfveto << endl;
    outfile << "inclusive analysis : " << inclusive << endl;
    outfile << endl;
    outfile << "####################################" << endl;
    outfile << "N of generated events" << endl;
    outfile << "####################################" << endl;
    outfile << "# events hig_gg2012       =    "  << n_mc_2012[26] << endl;
    outfile << "# events hig_vbf2012      =    "  << n_mc_2012[27] << endl;
    outfile << "# events hig_w2012, wlept =    "  << n_mc_2012[28] << endl;
    outfile << "# events hig_w2012, whad  =    "  << n_mc_2012[29] << endl;
    outfile << "# events hig_z2012, zlept =    "  << n_mc_2012[30] << endl;
    outfile << "# events hig_z2012, zhad  =    "  << n_mc_2012[31] << endl;
    outfile << "# events hig_z2012, zneut =    "  << n_mc_2012[32] << endl;
    outfile << "# events dy_2012 =        "       << n_mc_2012[10] << endl;
    outfile << "# events box_2012, 10-25 =     "  << n_mc_2012[0] << endl;
    outfile << "# events box_2012, 25-250 =    "  << n_mc_2012[1] << endl;
    outfile << "# events box_2012, >250 =      "  << n_mc_2012[2] << endl;
    outfile << "# events diphotjet_2012, 10-25 = "   << n_mc_2012[3] << endl;
    outfile << "# events diphotjet_2012, 25-250 = "  << n_mc_2012[4] << endl;
    outfile << "# events diphotjet_2012, >250 = "    << n_mc_2012[5] << endl;
    outfile << "# events gjet_2012, 20-40 =      "   << n_mc_2012[6] << endl;
    outfile << "# events gjet_2012, >40 =      "     << n_mc_2012[7] << endl;
    outfile << "# events qcd_2012 =       "  << n_mc_2012[8]  << endl;
    outfile << "# events qcd2_2012 =      "  << n_mc_2012[9]  << endl;
    outfile << "# events wgamma_ele_2012 = " << n_mc_2012[11] << endl;
    outfile << "# events wgamma_mu_2012  = " << n_mc_2012[12] << endl;
    outfile << "# events wgamma_tau_2012 = " << n_mc_2012[13] << endl;
    outfile << "# events zgamma_ele_2012 = " << n_mc_2012[14] << endl;
    outfile << "# events zgamma_mu_2012  = " << n_mc_2012[15] << endl;
    outfile << "# events zgamma_tau_2012 = " << n_mc_2012[16] << endl;
    outfile << "# events ggWminus_2012 = "   << n_mc_2012[17] << endl;
    outfile << "# events ggWplus_2012 = "    << n_mc_2012[18] << endl;
    outfile << "# events ggZ_2012 = "        << n_mc_2012[19] << endl;
    outfile << "# events ggtt_2012 = "       << n_mc_2012[20] << endl;
    outfile << "# events ttjets_2012 = "     << n_mc_2012[21] << endl;
    outfile << "# events Wjets_2012 = "      << n_mc_2012[22] << endl;
    outfile << "# events WW_2012 = "         << n_mc_2012[23] << endl;
    outfile << "# events WZ_2012 = "         << n_mc_2012[24] << endl;
    outfile << "# events ZZ_2012 = "         << n_mc_2012[25] << endl;
    outfile << endl;
    outfile << "N of selected events (not scaled for lumi)" << endl;  
    outfile << "##########################################" << endl;
    outfile << "# events hig_gg2012       =    " << num_uns_mc_2012[26] << endl;
    outfile << "# events hig_vbf2012      =    " << num_uns_mc_2012[27] << endl;
    outfile << "# events hig_w2012, wlept =    " << num_uns_mc_2012[28] << endl;
    outfile << "# events hig_w2012, whad  =    " << num_uns_mc_2012[29] << endl;
    outfile << "# events hig_z2012, zlept =    " << num_uns_mc_2012[30] << endl;
    outfile << "# events hig_z2012, zhad  =    " << num_uns_mc_2012[31] << endl;
    outfile << "# events hig_z2012, zneut =    " << num_uns_mc_2012[32] << endl;
    outfile << "# events dy_2012 =        "      << num_uns_mc_2012[10] << endl;
    outfile << "# events box_2012, 10-25 =     " << num_uns_mc_2012[0]  << endl;
    outfile << "# events box_2012, 25-250 =    " << num_uns_mc_2012[1]  << endl;
    outfile << "# events box_2012, >250 =      " << num_uns_mc_2012[2]  << endl;
    outfile << "# events diphotjet_2012, 10-25 = "   << num_uns_mc_2012[3] << endl;
    outfile << "# events diphotjet_2012, 25-250 = "  << num_uns_mc_2012[4] << endl;
    outfile << "# events diphotjet_2012, >250 = "    << num_uns_mc_2012[5] << endl;
    outfile << "# events gjet_2012, 20-40 =      "   << num_uns_mc_2012[6] << endl;
    outfile << "# events gjet_2012, >40 =      "     << num_uns_mc_2012[7] << endl;
    outfile << "# events qcd_2012 =       "  << num_uns_mc_2012[8] << endl;
    outfile << "# events qcd2_2012 =      "  << num_uns_mc_2012[9] << endl;
    outfile << "# events wgamma_ele_2012 = " << num_uns_mc_2012[11] << endl;
    outfile << "# events wgamma_mu_2012  = " << num_uns_mc_2012[12] << endl;
    outfile << "# events wgamma_tau_2012 = " << num_uns_mc_2012[13] << endl;
    outfile << "# events zgamma_ele_2012 = " << num_uns_mc_2012[14] << endl;
    outfile << "# events zgamma_mu_2012  = " << num_uns_mc_2012[15] << endl;
    outfile << "# events zgamma_tau_2012 = " << num_uns_mc_2012[16] << endl;
    outfile << "# events ggWminus_2012 = "   << num_uns_mc_2012[17] << endl;
    outfile << "# events ggWplus_2012 = "    << num_uns_mc_2012[18] << endl;
    outfile << "# events ggZ_2012 = "        << num_uns_mc_2012[19] << endl;
    outfile << "# events ggtt_2012 = "       << num_uns_mc_2012[20] << endl;
    outfile << "# events ttjets_2012 = "     << num_uns_mc_2012[21] << endl;
    outfile << "# events Wjets_2012 = "      << num_uns_mc_2012[22] << endl;
    outfile << "# events WW_2012 = "         << num_uns_mc_2012[23] << endl;
    outfile << "# events WZ_2012 = "         << num_uns_mc_2012[24] << endl;
    outfile << "# events ZZ_2012 = "         << num_uns_mc_2012[25] << endl;
    outfile << endl;
    outfile << "####################################" << endl;
    outfile << "N of selected events and eff." << endl;
    outfile << "####################################" << endl; 
    outfile << "ndata      = " << num_data << endl;
    outfile << endl;

    double num_bkg(0), err_num_bkg(0);
    double num_mc_total[33],num_uns_mc_total[33], n_mc_total[33];
    double err_num_mc_total[33],err_num_uns_mc_total[33];
    for (int i=0; i<33; i++){
      if(i<26){
	num_bkg += num_mc_2012[i];
      }
      num_mc_total[i] = num_mc_2012[i];
      num_uns_mc_total[i] = num_uns_mc_2012[i];
      err_num_uns_mc_total[i] = sqrt(num_uns_mc_total[i]);
      if(num_uns_mc_total[i]) err_num_mc_total[i] = err_num_uns_mc_total[i] * num_mc_total[i]/num_uns_mc_total[i];
      else err_num_mc_total[i] = 0;
      n_mc_total[i] = n_mc_2012[i] * scale_mc_2012[i];
    }
    for (int i=0; i<26; i++){
      err_num_bkg = sqrt(err_num_bkg*err_num_bkg + err_num_mc_total[i]*err_num_mc_total[i]);
    }

    outfile << "nallbkg      = " << num_bkg << " +/- " << err_num_bkg << endl;
    outfile << "nhig glu     = " << num_mc_total[26] << " +/- " << err_num_mc_total[26] << endl;
    outfile << "nhig vbf     = " << num_mc_total[27] << " +/- " << err_num_mc_total[27] << endl;
    outfile << "nhig wlept   = " << num_mc_total[28] << " +/- " << err_num_mc_total[28] << endl;
    outfile << "nhig whad    = " << num_mc_total[29] << " +/- " << err_num_mc_total[29] << endl;
    outfile << "nhig zlept   = " << num_mc_total[30] << " +/- " << err_num_mc_total[30] << endl;
    outfile << "nhig zhad    = " << num_mc_total[31] << " +/- " << err_num_mc_total[31] << endl;
    outfile << "nhig zneut   = " << num_mc_total[32] << " +/- " << err_num_mc_total[32] << endl;
    outfile << "ndy          = " << num_mc_total[10] << " +/- " << err_num_mc_total[10] << endl;
    outfile << "nbox,10-25   = " << num_mc_total[0]   << " +/- " << err_num_mc_total[0]  << endl;
    outfile << "nbox,25-250  = " << num_mc_total[1]   << " +/- " << err_num_mc_total[1]  << endl;
    outfile << "nbox,>250    = " << num_mc_total[2]   << " +/- " << err_num_mc_total[2]  << endl;
    outfile << "ndiphot,10-25  = " << num_mc_total[3] << " +/- " << err_num_mc_total[3]  << endl;
    outfile << "ndiphot,25-250 = " << num_mc_total[4] << " +/- " << err_num_mc_total[4]  << endl;
    outfile << "ndiphot,>250   = " << num_mc_total[5] << " +/- " << err_num_mc_total[5]  << endl;
    outfile << "ngjet,20-40    = " << num_mc_total[6]  << " +/- " << err_num_mc_total[6]  << endl;
    outfile << "ngjet,>40      = " << num_mc_total[7]  << " +/- " << err_num_mc_total[7]  << endl;
    outfile << "nqcd40       = " << num_mc_total[8]   << " +/- " << err_num_mc_total[8]   << endl;
    outfile << "nqcd30-40    = " << num_mc_total[9]   << " +/- " << err_num_mc_total[9]   << endl;
    outfile << "nwg-ele      = " << num_mc_total[11]  << " +/- " << err_num_mc_total[11]  << endl;
    outfile << "nwg-mu       = " << num_mc_total[12]  << " +/- " << err_num_mc_total[12]  << endl;
    outfile << "nwg-tau      = " << num_mc_total[13]  << " +/- " << err_num_mc_total[13]  << endl;
    outfile << "nzg-ele      = " << num_mc_total[14]  << " +/- " << err_num_mc_total[14]  << endl;
    outfile << "nzg-mu       = " << num_mc_total[15] << " +/- " << err_num_mc_total[15] << endl;
    outfile << "nzg-tau      = " << num_mc_total[16] << " +/- " << err_num_mc_total[16] << endl;
    outfile << "nggWminus    = " << num_mc_total[17] << " +/- " << err_num_mc_total[17] << endl;
    outfile << "nggWplus     = " << num_mc_total[18] << " +/- " << err_num_mc_total[18] << endl;
    outfile << "nggZ         = " << num_mc_total[19] << " +/- " << err_num_mc_total[19] << endl;
    outfile << "nggtt        = " << num_mc_total[20] << " +/- " << err_num_mc_total[20] << endl;
    outfile << "nttjets      = " << num_mc_total[21] << " +/- " << err_num_mc_total[21] << endl;
    outfile << "nWjets       = " << num_mc_total[22] << " +/- " << err_num_mc_total[22] << endl;
    outfile << "nWW          = " << num_mc_total[23] << " +/- " << err_num_mc_total[23] << endl;
    outfile << "nWZ          = " << num_mc_total[24] << " +/- " << err_num_mc_total[24] << endl;
    outfile << "nZZ          = " << num_mc_total[25] << " +/- " << err_num_mc_total[25] << endl;
    outfile << endl;              
    outfile << "eff nhig     = " 
	    << (num_mc_total[26] + num_mc_total[27] + num_mc_total[28] + num_mc_total[29] + num_mc_total[30] + num_mc_total[31] + num_mc_total[32]) / (n_mc_total[26] + n_mc_total[27] + n_mc_total[28] + n_mc_total[29] + n_mc_total[30] + n_mc_total[31] + n_mc_total[32]) << endl;

    outfile << "eff nhig glu    = " << num_mc_total[26] / n_mc_total[26] << endl;
    outfile << "eff nhig vbf    = " << num_mc_total[27] / n_mc_total[27] << endl;
    outfile << "eff nhig wlept  = " << num_mc_total[28] / n_mc_total[28] << endl;
    outfile << "eff nhig whad   = " << num_mc_total[29] / n_mc_total[29] << endl;
    outfile << "eff nhig zlept  = " << num_mc_total[30] / n_mc_total[30] << endl;
    outfile << "eff nhig zhad   = " << num_mc_total[31] / n_mc_total[31] << endl;
    outfile << "eff nhig zneut  = " << num_mc_total[32] / n_mc_total[32] << endl;
    outfile << "eff ndy       = " << num_mc_total[10] / n_mc_total[10] << endl;
    outfile << "eff nbox 1     = " << num_mc_total[0] / n_mc_total[0] << endl;
    outfile << "eff nbox 2     = " << num_mc_total[1] / n_mc_total[1] << endl;
    outfile << "eff nbox 3     = " << num_mc_total[2] / n_mc_total[2] << endl;
    outfile << "eff ndiphot 1  = " << num_mc_total[3] / n_mc_total[3] << endl;
    outfile << "eff ndiphot 2  = " << num_mc_total[4] / n_mc_total[4] << endl;
    outfile << "eff ndiphot 3  = " << num_mc_total[5] / n_mc_total[5] << endl;
    outfile << "eff ngjet 1    = " << num_mc_total[6] / n_mc_total[6] << endl;
    outfile << "eff ngjet 2    = " << num_mc_total[7] / n_mc_total[7] << endl;
    outfile << "eff nqcd      = " << num_mc_total[8] / n_mc_total[8] << endl;
    outfile << "eff nqcd30-40 = " << num_mc_total[9] / n_mc_total[9] << endl;
    outfile << "eff wgamma-ele = " << num_mc_total[11] / n_mc_total[11] << endl;
    outfile << "eff wgamma-mu = "  << num_mc_total[12] / n_mc_total[12] << endl;
    outfile << "eff wgamma-tau = " << num_mc_total[13] / n_mc_total[13] << endl;
    outfile << "eff zgamma-ele = " << num_mc_total[14] / n_mc_total[14] << endl;
    outfile << "eff zgamma-mu = "  << num_mc_total[15] / n_mc_total[15] << endl;
    outfile << "eff zgamma-tau = " << num_mc_total[16] / n_mc_total[16] << endl;
    outfile << "eff ggWminus = " << num_mc_total[17] / n_mc_total[17] << endl;    
    outfile << "eff ggWplus  = " << num_mc_total[18] / n_mc_total[18] << endl;    
    outfile << "eff ggZ  = "     << num_mc_total[19] / n_mc_total[19] << endl;    
    outfile << "eff ggtt  = "    << num_mc_total[20] / n_mc_total[20] << endl;    
    outfile << "eff ttjets  = "  << num_mc_total[21] / n_mc_total[21] << endl;    
    outfile << "eff Wjets  = "   << num_mc_total[22] / n_mc_total[22] << endl;    
    outfile << "eff WW = "       << num_mc_total[23] / n_mc_total[23] << endl;    
    outfile << "eff WZ = "       << num_mc_total[24] / n_mc_total[24] << endl;    
    outfile << "eff ZZ = "       << num_mc_total[25] / n_mc_total[25] << endl;    

    outfile << endl;
    for(int i=0; i<11; i++){  
      outfile << "eff nhig gluglu "<< h_masses[i] << " = " << num_gluglu_2012[i] / n_gluglu_2012[i] << endl;
      outfile << "eff nhig vbf    "<< h_masses[i] << " = " << num_vbf_2012[i] / n_vbf_2012[i] << endl;
      outfile << "eff nhig wzh    "<< h_masses[i] << " = " << num_wzh_2012[i] / n_wzh_2012[i] << endl;
      //
      outfile << "nhig gluglu "<< h_masses[i] << " = " << num_gluglu_2012[i] << endl;
      outfile << "nhig vbf    "<< h_masses[i] << " = " << num_vbf_2012[i]    << endl;
      outfile << "nhig wzh    "<< h_masses[i] << " = " << num_wzh_2012[i]    << endl;
    }

    outfile.close();
  }

  hOutputFile->Write() ;
  hOutputFile->Close() ;
  hOutputFile->Delete();

  delete data;
  //delete datacs;

  for(int i=0; i<33; i++){
    delete mc_2012_fill[i];
    delete mc_2012[i];
  }

  for(int i=0; i<11; i++){  
    delete mc_gluglu_2012_fill[i];
    delete mc_vbf_2012_fill[i];
    delete mc_wzh_2012_fill[i];
    delete mc_gluglu_2012[i];
    delete mc_vbf_2012[i];
    delete mc_wzh_2012[i];
  }

  vector<double> values;
  if (variableData == "massgg"){
    values.push_back(integralhiggs);
    values.push_back(integralbkg);
    values.push_back(integralhiggs/sqrt(integralbkg));
  }

  return values;
}








