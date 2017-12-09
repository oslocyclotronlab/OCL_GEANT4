#define AnalyseSims_cxx
#include "AnalyseSims.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void AnalyseSims::Loop()
{
//*****************************************************************************//
// Example class for analysing OSCAR-data simulated with GEANT4. As long as the
// NTuple for storing data is not changed, you can use this script as a starting point. 
// As an example, the script now plots the electron energy deposited in the detectors
// and the folded energy together. 
//   
//   In a ROOT session, you can do:
//      root> .L AnalyseSims.C
//      root> AnalyseSims t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch

   fChain->SetBranchStatus("*",1);  // enable all branches, practical since we need all branches in this analysis

   if (fChain == 0) return; //aborts if it cannot find the tree/chain of trees

   TH1D *h1 = new TH1D("h1","Simulated electron energy deposition",2000,0,10);
   TH1D *h2 = new TH1D("h2","Simulated and folded energy deposition",2000,0,10); 

   Double_t smoothingFactor = 0.00019;

   TRandom *r3 = new TRandom3(0.);

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      if(EdepInCrystal1  > 0.) h1->Fill(EdepInCrystal1);
      if(EdepInCrystal2  > 0.) h1->Fill(EdepInCrystal2);
      if(EdepInCrystal3  > 0.) h1->Fill(EdepInCrystal3);
      if(EdepInCrystal4  > 0.) h1->Fill(EdepInCrystal4);
      if(EdepInCrystal5  > 0.) h1->Fill(EdepInCrystal5);
      if(EdepInCrystal6  > 0.) h1->Fill(EdepInCrystal6);
      if(EdepInCrystal7  > 0.) h1->Fill(EdepInCrystal7);
      if(EdepInCrystal8  > 0.) h1->Fill(EdepInCrystal8);
      if(EdepInCrystal9  > 0.) h1->Fill(EdepInCrystal9);
      if(EdepInCrystal10 > 0.) h1->Fill(EdepInCrystal10);
      if(EdepInCrystal11 > 0.) h1->Fill(EdepInCrystal11);
      if(EdepInCrystal12 > 0.) h1->Fill(EdepInCrystal12);
      if(EdepInCrystal13 > 0.) h1->Fill(EdepInCrystal13);
      if(EdepInCrystal14 > 0.) h1->Fill(EdepInCrystal14);
      if(EdepInCrystal15 > 0.) h1->Fill(EdepInCrystal15);
      if(EdepInCrystal16 > 0.) h1->Fill(EdepInCrystal16);
      if(EdepInCrystal17 > 0.) h1->Fill(EdepInCrystal17);
      if(EdepInCrystal18 > 0.) h1->Fill(EdepInCrystal18);
      if(EdepInCrystal19 > 0.) h1->Fill(EdepInCrystal19);
      if(EdepInCrystal20 > 0.) h1->Fill(EdepInCrystal20);
      if(EdepInCrystal21 > 0.) h1->Fill(EdepInCrystal21);
      if(EdepInCrystal22 > 0.) h1->Fill(EdepInCrystal22);
      if(EdepInCrystal23 > 0.) h1->Fill(EdepInCrystal23);
      if(EdepInCrystal24 > 0.) h1->Fill(EdepInCrystal24);
      if(EdepInCrystal25 > 0.) h1->Fill(EdepInCrystal25);
      if(EdepInCrystal26 > 0.) h1->Fill(EdepInCrystal26);
      if(EdepInCrystal27 > 0.) h1->Fill(EdepInCrystal27);
      if(EdepInCrystal28 > 0.) h1->Fill(EdepInCrystal28);
      if(EdepInCrystal29 > 0.) h1->Fill(EdepInCrystal29);
      if(EdepInCrystal30 > 0.) h1->Fill(EdepInCrystal30);

      if(EdepInCrystal1  > 0.) h2->Fill(r3->Gaus(EdepInCrystal1 , sqrt(smoothingFactor*EdepInCrystal1 )));
      if(EdepInCrystal2  > 0.) h2->Fill(r3->Gaus(EdepInCrystal2 , sqrt(smoothingFactor*EdepInCrystal2 )));
      if(EdepInCrystal3  > 0.) h2->Fill(r3->Gaus(EdepInCrystal3 , sqrt(smoothingFactor*EdepInCrystal3 )));
      if(EdepInCrystal4  > 0.) h2->Fill(r3->Gaus(EdepInCrystal4 , sqrt(smoothingFactor*EdepInCrystal4 )));
      if(EdepInCrystal5  > 0.) h2->Fill(r3->Gaus(EdepInCrystal5 , sqrt(smoothingFactor*EdepInCrystal5 )));
      if(EdepInCrystal6  > 0.) h2->Fill(r3->Gaus(EdepInCrystal6 , sqrt(smoothingFactor*EdepInCrystal6 )));
      if(EdepInCrystal7  > 0.) h2->Fill(r3->Gaus(EdepInCrystal7 , sqrt(smoothingFactor*EdepInCrystal7 )));
      if(EdepInCrystal8  > 0.) h2->Fill(r3->Gaus(EdepInCrystal8 , sqrt(smoothingFactor*EdepInCrystal8 )));
      if(EdepInCrystal9  > 0.) h2->Fill(r3->Gaus(EdepInCrystal9 , sqrt(smoothingFactor*EdepInCrystal9 )));
      if(EdepInCrystal10 > 0.) h2->Fill(r3->Gaus(EdepInCrystal10, sqrt(smoothingFactor*EdepInCrystal10)));
      if(EdepInCrystal11 > 0.) h2->Fill(r3->Gaus(EdepInCrystal11, sqrt(smoothingFactor*EdepInCrystal11)));
      if(EdepInCrystal12 > 0.) h2->Fill(r3->Gaus(EdepInCrystal12, sqrt(smoothingFactor*EdepInCrystal12)));
      if(EdepInCrystal13 > 0.) h2->Fill(r3->Gaus(EdepInCrystal13, sqrt(smoothingFactor*EdepInCrystal13)));
      if(EdepInCrystal14 > 0.) h2->Fill(r3->Gaus(EdepInCrystal14, sqrt(smoothingFactor*EdepInCrystal14)));
      if(EdepInCrystal15 > 0.) h2->Fill(r3->Gaus(EdepInCrystal15, sqrt(smoothingFactor*EdepInCrystal15)));
      if(EdepInCrystal16 > 0.) h2->Fill(r3->Gaus(EdepInCrystal16, sqrt(smoothingFactor*EdepInCrystal16)));
      if(EdepInCrystal17 > 0.) h2->Fill(r3->Gaus(EdepInCrystal17, sqrt(smoothingFactor*EdepInCrystal17)));
      if(EdepInCrystal18 > 0.) h2->Fill(r3->Gaus(EdepInCrystal18, sqrt(smoothingFactor*EdepInCrystal18)));
      if(EdepInCrystal19 > 0.) h2->Fill(r3->Gaus(EdepInCrystal19, sqrt(smoothingFactor*EdepInCrystal19)));
      if(EdepInCrystal20 > 0.) h2->Fill(r3->Gaus(EdepInCrystal20, sqrt(smoothingFactor*EdepInCrystal20)));
      if(EdepInCrystal21 > 0.) h2->Fill(r3->Gaus(EdepInCrystal21, sqrt(smoothingFactor*EdepInCrystal21)));
      if(EdepInCrystal22 > 0.) h2->Fill(r3->Gaus(EdepInCrystal22, sqrt(smoothingFactor*EdepInCrystal22)));
      if(EdepInCrystal23 > 0.) h2->Fill(r3->Gaus(EdepInCrystal23, sqrt(smoothingFactor*EdepInCrystal23)));
      if(EdepInCrystal24 > 0.) h2->Fill(r3->Gaus(EdepInCrystal24, sqrt(smoothingFactor*EdepInCrystal24)));
      if(EdepInCrystal25 > 0.) h2->Fill(r3->Gaus(EdepInCrystal25, sqrt(smoothingFactor*EdepInCrystal25)));
      if(EdepInCrystal26 > 0.) h2->Fill(r3->Gaus(EdepInCrystal26, sqrt(smoothingFactor*EdepInCrystal26)));
      if(EdepInCrystal27 > 0.) h2->Fill(r3->Gaus(EdepInCrystal27, sqrt(smoothingFactor*EdepInCrystal27)));
      if(EdepInCrystal28 > 0.) h2->Fill(r3->Gaus(EdepInCrystal28, sqrt(smoothingFactor*EdepInCrystal28)));
      if(EdepInCrystal29 > 0.) h2->Fill(r3->Gaus(EdepInCrystal29, sqrt(smoothingFactor*EdepInCrystal29)));
      if(EdepInCrystal30 > 0.) h2->Fill(r3->Gaus(EdepInCrystal30, sqrt(smoothingFactor*EdepInCrystal30)));

      // if (Cut(ientry) < 0) continue;
   }
   h1->SetLineColor(kBlack);
   h1->SetLineStyle(4);
   h1->Draw("hist");
   h2->Draw("same hist");
}
