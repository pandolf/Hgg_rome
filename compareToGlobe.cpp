#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"

#include <iostream>
#include <cmath>





int main() {

  TFile* file_globe = TFile::Open("/afs/cern.ch/work/p/pandolf/CMSSW_5_3_9_patch1_globe_synch/src/HiggsAnalysis/HiggsTo2photons/h2gglobe/AnalysisScripts/histograms_CMS-HGG_data_2.root");
  TTree* tree_globe = (TTree*)file_globe->Get("Data");

  TFile* file_rome = TFile::Open("/afs/cern.ch/work/p/pandolf/CMSSW_5_2_5/src/Analysis/Higgs/redntp_data_prova_synch.root");
  TTree* tree_rome = (TTree*)file_rome->Get("AnaTree");


  int event_globe;
  tree_globe->SetBranchAddress("event", &event_globe);
  float ptPhot1_globe; 
  tree_globe->SetBranchAddress("ph1_pt", &ptPhot1_globe);
  float ptPhot2_globe; 
  tree_globe->SetBranchAddress("ph2_pt", &ptPhot2_globe);
  float etaPhot1_globe; 
  tree_globe->SetBranchAddress("ph1_eta", &etaPhot1_globe);
  float etaPhot2_globe; 
  tree_globe->SetBranchAddress("ph2_eta", &etaPhot2_globe);
  float ePhot1_globe; 
  tree_globe->SetBranchAddress("ph1_e", &ePhot1_globe);
  float ePhot2_globe; 
  tree_globe->SetBranchAddress("ph2_e", &ePhot2_globe);
  float mgg_globe; 
  tree_globe->SetBranchAddress("PhotonsMass", &mgg_globe);
  float ptJet1_globe; 
  tree_globe->SetBranchAddress("j1_pt", &ptJet1_globe);
  float ptJet2_globe; 
  tree_globe->SetBranchAddress("j2_pt", &ptJet2_globe);
  float ptJet3_globe; 
  tree_globe->SetBranchAddress("j3_pt", &ptJet3_globe);
  float ptJet4_globe; 
  tree_globe->SetBranchAddress("j4_pt", &ptJet4_globe);
  float etaJet1_globe; 
  tree_globe->SetBranchAddress("j1_eta", &etaJet1_globe);
  float etaJet2_globe; 
  tree_globe->SetBranchAddress("j2_eta", &etaJet2_globe);
  float etaJet3_globe; 
  tree_globe->SetBranchAddress("j3_eta", &etaJet3_globe);
  float etaJet4_globe; 
  tree_globe->SetBranchAddress("j4_eta", &etaJet4_globe);
  float csvJet1_globe; 
  tree_globe->SetBranchAddress("j1_algoPF1_csvBtag", &csvJet1_globe);
  float csvJet2_globe; 
  tree_globe->SetBranchAddress("j2_algoPF1_csvBtag", &csvJet2_globe);
  float csvJet3_globe; 
  tree_globe->SetBranchAddress("j3_algoPF1_csvBtag", &csvJet3_globe);
  float csvJet4_globe; 
  tree_globe->SetBranchAddress("j4_algoPF1_csvBtag", &csvJet4_globe);
  float betastarJet1_globe; 
  tree_globe->SetBranchAddress("j1_betaStarClassic", &betastarJet1_globe);
  float betastarJet2_globe; 
  tree_globe->SetBranchAddress("j2_betaStarClassic", &betastarJet2_globe);
  float betastarJet3_globe; 
  tree_globe->SetBranchAddress("j3_betaStarClassic", &betastarJet3_globe);
  float betastarJet4_globe; 
  tree_globe->SetBranchAddress("j4_betaStarClassic", &betastarJet4_globe);
  float rmscandJet1_globe; 
  tree_globe->SetBranchAddress("j1_rmscand", &rmscandJet1_globe);
  float rmscandJet2_globe; 
  tree_globe->SetBranchAddress("j2_rmscand", &rmscandJet2_globe);
  float rmscandJet3_globe; 
  tree_globe->SetBranchAddress("j3_rmscand", &rmscandJet3_globe);
  float rmscandJet4_globe; 
  tree_globe->SetBranchAddress("j4_rmscand", &rmscandJet4_globe);


  int event_rome;
  tree_rome->SetBranchAddress("event", &event_rome);
  float ptPhot1_rome; 
  tree_rome->SetBranchAddress("ptphot1", &ptPhot1_rome);
  float ptPhot2_rome; 
  tree_rome->SetBranchAddress("ptphot2", &ptPhot2_rome);
  float etaPhot1_rome; 
  tree_rome->SetBranchAddress("etaphot1", &etaPhot1_rome);
  float etaPhot2_rome; 
  tree_rome->SetBranchAddress("etaphot2", &etaPhot2_rome);
  float ePhot1_rome; 
  tree_rome->SetBranchAddress("ephot1", &ePhot1_rome);
  float ePhot2_rome; 
  tree_rome->SetBranchAddress("ephot2", &ePhot2_rome);
  float mgg_rome; 
  tree_rome->SetBranchAddress("massggnewvtx", &mgg_rome);
  int njets_rome;
  tree_rome->SetBranchAddress("njets", &njets_rome);
  float ptJet_rome[10]; 
  tree_rome->SetBranchAddress("ptcorrjet", ptJet_rome);
  float etaJet_rome[10]; 
  tree_rome->SetBranchAddress("etajet", etaJet_rome);
  float betastarJet_rome[10]; 
  tree_rome->SetBranchAddress("betastarjet", betastarJet_rome);
  float rmscandJet_rome[10]; 
  tree_rome->SetBranchAddress("rmsjet", rmscandJet_rome);
  float csvJet_rome[10]; 
  tree_rome->SetBranchAddress("btagcsvjet", csvJet_rome);



  TFile* outfile = TFile::Open("comparison.root", "recreate");
  outfile->cd();

  TH1D* h1_ptPhot1_diff = new TH1D("ptPhot1_diff", "", 100, -0.5, 0.5);
  TH1D* h1_ptPhot2_diff = new TH1D("ptPhot2_diff", "", 100, -0.5, 0.5);

  TH1D* h1_etaPhot1_diff = new TH1D("etaPhot1_diff", "", 100, -0.1, 0.1);
  TH1D* h1_etaPhot2_diff = new TH1D("etaPhot2_diff", "", 100, -0.1, 0.1);

  TH1D* h1_ePhot1_diff = new TH1D("ePhot1_diff", "", 100, -1., 1.);
  TH1D* h1_ePhot2_diff = new TH1D("ePhot2_diff", "", 100, -1., 1.);

  TH1D* h1_mgg_diff = new TH1D("mgg_diff", "", 100, -1., 1.);

  TH1D* h1_ptJet1_diff = new TH1D("ptJet1_diff", "", 100, -3., 3.);
  TH1D* h1_ptJet2_diff = new TH1D("ptJet2_diff", "", 100, -3., 3.);
  TH1D* h1_ptJet3_diff = new TH1D("ptJet3_diff", "", 100, -3., 3.);
  TH1D* h1_ptJet4_diff = new TH1D("ptJet4_diff", "", 100, -3., 3.);

  TH1D* h1_etaJet1_diff = new TH1D("etaJet1_diff", "", 100, -0.1, 0.1);
  TH1D* h1_etaJet2_diff = new TH1D("etaJet2_diff", "", 100, -0.1, 0.1);
  TH1D* h1_etaJet3_diff = new TH1D("etaJet3_diff", "", 100, -0.1, 0.1);
  TH1D* h1_etaJet4_diff = new TH1D("etaJet4_diff", "", 100, -0.1, 0.1);

  TH1D* h1_betastarJet1_diff = new TH1D("betastarJet1_diff", "", 100, -3., 3.);
  TH1D* h1_betastarJet2_diff = new TH1D("betastarJet2_diff", "", 100, -3., 3.);
  TH1D* h1_betastarJet3_diff = new TH1D("betastarJet3_diff", "", 100, -3., 3.);
  TH1D* h1_betastarJet4_diff = new TH1D("betastarJet4_diff", "", 100, -3., 3.);

  TH1D* h1_csvJet1_diff = new TH1D("csvJet1_diff", "", 100, -3., 3.);
  TH1D* h1_csvJet2_diff = new TH1D("csvJet2_diff", "", 100, -3., 3.);
  TH1D* h1_csvJet3_diff = new TH1D("csvJet3_diff", "", 100, -3., 3.);
  TH1D* h1_csvJet4_diff = new TH1D("csvJet4_diff", "", 100, -3., 3.);



  int nentries_globe = tree_globe->GetEntries();
  int nentries_rome = tree_rome->GetEntries();

  for( unsigned ientry_globe=0; ientry_globe<nentries_globe; ++ientry_globe ) {

    tree_globe->GetEntry(ientry_globe);

    bool foundInRome = false;

    for( unsigned ientry_rome=0; ientry_rome<nentries_rome; ++ientry_rome ) {

      tree_rome->GetEntry(ientry_rome);

      if( event_rome!=event_globe ) continue;

      bool printEvent = false;
      foundInRome = true; 

      h1_ptPhot1_diff->Fill( ptPhot1_globe-ptPhot1_rome );
      if( fabs(ptPhot1_globe-ptPhot1_rome)>0.5 ) printEvent=true;
      h1_ptPhot2_diff->Fill( ptPhot2_globe-ptPhot2_rome );
      if( fabs(ptPhot2_globe-ptPhot2_rome)>0.5 ) printEvent=true;

      h1_etaPhot1_diff->Fill( etaPhot1_globe-etaPhot1_rome );
      h1_etaPhot2_diff->Fill( etaPhot2_globe-etaPhot2_rome );

      h1_ePhot1_diff->Fill( ePhot1_globe-ePhot1_rome );
      h1_ePhot2_diff->Fill( ePhot2_globe-ePhot2_rome );

      h1_mgg_diff->Fill( mgg_globe-mgg_rome );
      if( fabs(mgg_globe-mgg_rome)>0.5 ) printEvent=true;


      float ptJet1_rome;
      float ptJet2_rome;
      float ptJet3_rome;
      float ptJet4_rome;

      float etaJet1_rome;
      float etaJet2_rome;
      float etaJet3_rome;
      float etaJet4_rome;

      float csvJet1_rome;
      float csvJet2_rome;
      float csvJet3_rome;
      float csvJet4_rome;

      float betastarJet1_rome;
      float betastarJet2_rome;
      float betastarJet3_rome;
      float betastarJet4_rome;

      if( njets_rome>0 ) {
        ptJet1_rome = ptJet_rome[0];
        etaJet1_rome = etaJet_rome[0];
        csvJet1_rome = csvJet_rome[0];
        betastarJet1_rome = betastarJet_rome[0];
      } else {
        ptJet1_rome = -999.;
        etaJet1_rome = -999.;
        csvJet1_rome = -999.;
        betastarJet1_rome = -999.;
      }

      if( njets_rome>1 ) {
        ptJet2_rome = ptJet_rome[1];
        etaJet2_rome = etaJet_rome[1];
        csvJet2_rome = csvJet_rome[1];
        betastarJet2_rome = betastarJet_rome[1];
      } else {
        ptJet2_rome = -999.;
        etaJet2_rome = -999.;
        csvJet2_rome = -999.;
        betastarJet2_rome = -999.;
      }
      if( njets_rome>2 ) {
        ptJet3_rome = ptJet_rome[2];
        etaJet3_rome = etaJet_rome[2];
        csvJet3_rome = csvJet_rome[2];
        betastarJet3_rome = betastarJet_rome[2];
      } else {
        ptJet3_rome = -999.;
        etaJet3_rome = -999.;
        csvJet3_rome = -999.;
        betastarJet3_rome = -999.;
      }
      if( njets_rome>3 ) {
        ptJet4_rome = ptJet_rome[3];
        etaJet4_rome = etaJet_rome[3];
        csvJet4_rome = csvJet_rome[3];
        betastarJet4_rome = betastarJet_rome[3];
      } else {
        ptJet4_rome = -999.;
        etaJet4_rome = -999.;
        csvJet4_rome = -999.;
        betastarJet4_rome = -999.;
      }


//std::cout << std::endl << "event: " << event_globe << std::endl;
//std::cout << "ptjet1: globe: " << ptJet1_globe << " rome: " << ptJet1_rome << std::endl;
//std::cout << "ptjet2: globe: " << ptJet2_globe << " rome: " << ptJet2_rome << std::endl;
//std::cout << "ptjet3: globe: " << ptJet3_globe << " rome: " << ptJet3_rome << std::endl;
//std::cout << "ptjet4: globe: " << ptJet4_globe << " rome: " << ptJet4_rome << std::endl;

      h1_ptJet1_diff->Fill( ptJet1_globe-ptJet1_rome );
      if( fabs(ptJet1_globe-ptJet1_rome)>0.5 ) printEvent=true;
      h1_ptJet2_diff->Fill( ptJet2_globe-ptJet2_rome );
      if( fabs(ptJet2_globe-ptJet2_rome)>0.5 ) printEvent=true;
      h1_ptJet3_diff->Fill( ptJet3_globe-ptJet3_rome );
      if( fabs(ptJet3_globe-ptJet3_rome)>0.5 ) printEvent=true;
      h1_ptJet4_diff->Fill( ptJet4_globe-ptJet4_rome );
      if( fabs(ptJet4_globe-ptJet4_rome)>0.5 ) printEvent=true;

      h1_etaJet1_diff->Fill( etaJet1_globe-etaJet1_rome );
      h1_etaJet2_diff->Fill( etaJet2_globe-etaJet2_rome );
      h1_etaJet3_diff->Fill( etaJet3_globe-etaJet3_rome );
      h1_etaJet4_diff->Fill( etaJet4_globe-etaJet4_rome );

      if( betastarJet1_rome==-1 ) betastarJet1_rome=0.;
      if( betastarJet2_rome==-1 ) betastarJet2_rome=0.;
      if( betastarJet3_rome==-1 ) betastarJet3_rome=0.;
      if( betastarJet4_rome==-1 ) betastarJet4_rome=0.;

      h1_betastarJet1_diff->Fill( betastarJet1_globe-betastarJet1_rome );
      if( fabs(betastarJet4_globe-betastarJet4_rome)>0.02 ) printEvent=true;
      h1_betastarJet2_diff->Fill( betastarJet2_globe-betastarJet2_rome );
      if( fabs(betastarJet4_globe-betastarJet4_rome)>0.02 ) printEvent=true;
      h1_betastarJet3_diff->Fill( betastarJet3_globe-betastarJet3_rome );
      if( fabs(betastarJet4_globe-betastarJet4_rome)>0.02 ) printEvent=true;
      h1_betastarJet4_diff->Fill( betastarJet4_globe-betastarJet4_rome );
      if( fabs(betastarJet4_globe-betastarJet4_rome)>0.02 ) printEvent=true;

      h1_csvJet1_diff->Fill( csvJet1_globe-csvJet1_rome );
      h1_csvJet2_diff->Fill( csvJet2_globe-csvJet2_rome );
      h1_csvJet3_diff->Fill( csvJet3_globe-csvJet3_rome );
      h1_csvJet4_diff->Fill( csvJet4_globe-csvJet4_rome );

      if( printEvent ) std::cout << std::endl << "event: " << event_globe << " is suspicious." << std::endl;

      break;

    } // for rome

    if( !foundInRome ) std::cout << "Event " << event_globe << " not found in rome." << std::endl;

  } // for globe


  outfile->cd();

  h1_ptPhot1_diff->Write();
  h1_ptPhot2_diff->Write();

  h1_etaPhot1_diff->Write();
  h1_etaPhot2_diff->Write();

  h1_ePhot1_diff->Write();
  h1_ePhot2_diff->Write();

  h1_mgg_diff->Write();

  h1_ptJet1_diff->Write();
  h1_ptJet2_diff->Write();
  h1_ptJet3_diff->Write();
  h1_ptJet4_diff->Write();

  h1_etaJet1_diff->Write();
  h1_etaJet2_diff->Write();
  h1_etaJet3_diff->Write();
  h1_etaJet4_diff->Write();

  h1_betastarJet1_diff->Write();
  h1_betastarJet2_diff->Write();
  h1_betastarJet3_diff->Write();
  h1_betastarJet4_diff->Write();

  h1_csvJet1_diff->Write();
  h1_csvJet2_diff->Write();
  h1_csvJet3_diff->Write();
  h1_csvJet4_diff->Write();


  outfile->Close();

  return 0;

}
  
