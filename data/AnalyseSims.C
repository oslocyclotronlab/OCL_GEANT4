#define AnalyseSims_cxx
#include "AnalyseSims.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <fstream>

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

   TH1D *h1 = new TH1D("h1","Simulated electron energy deposition",4200,0,21);
   TH1D *h2 = new TH1D("h2","Simulated and folded energy deposition",4200,0,21);

   Double_t cSmooth[] = {2.03936976e-04, 6.82322078e-23,  3.76053110e-05}; // make sure cal is in same units (MeVor keV)

   // generate a random number
   TRandom *r3 = new TRandom3(0.); // initiated with random seed

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

                                        // use the random numb generator to sample from a gaussian with the parameters below
      if(EdepInCrystal1  > 0.) h2->Fill(r3->Gaus(EdepInCrystal1 , sqrt(cSmooth[0] + cSmooth[1]*EdepInCrystal1  + cSmooth[2]*EdepInCrystal1 *EdepInCrystal1 )));
      if(EdepInCrystal2  > 0.) h2->Fill(r3->Gaus(EdepInCrystal2 , sqrt(cSmooth[0] + cSmooth[1]*EdepInCrystal2  + cSmooth[2]*EdepInCrystal2 *EdepInCrystal2 )));
      if(EdepInCrystal3  > 0.) h2->Fill(r3->Gaus(EdepInCrystal3 , sqrt(cSmooth[0] + cSmooth[1]*EdepInCrystal3  + cSmooth[2]*EdepInCrystal3 *EdepInCrystal3 )));
      if(EdepInCrystal4  > 0.) h2->Fill(r3->Gaus(EdepInCrystal4 , sqrt(cSmooth[0] + cSmooth[1]*EdepInCrystal4  + cSmooth[2]*EdepInCrystal4 *EdepInCrystal4 )));
      if(EdepInCrystal5  > 0.) h2->Fill(r3->Gaus(EdepInCrystal5 , sqrt(cSmooth[0] + cSmooth[1]*EdepInCrystal5  + cSmooth[2]*EdepInCrystal5 *EdepInCrystal5 )));
      if(EdepInCrystal6  > 0.) h2->Fill(r3->Gaus(EdepInCrystal6 , sqrt(cSmooth[0] + cSmooth[1]*EdepInCrystal6  + cSmooth[2]*EdepInCrystal6 *EdepInCrystal6 )));
      if(EdepInCrystal7  > 0.) h2->Fill(r3->Gaus(EdepInCrystal7 , sqrt(cSmooth[0] + cSmooth[1]*EdepInCrystal7  + cSmooth[2]*EdepInCrystal7 *EdepInCrystal7 )));
      if(EdepInCrystal8  > 0.) h2->Fill(r3->Gaus(EdepInCrystal8 , sqrt(cSmooth[0] + cSmooth[1]*EdepInCrystal8  + cSmooth[2]*EdepInCrystal8 *EdepInCrystal8 )));
      if(EdepInCrystal9  > 0.) h2->Fill(r3->Gaus(EdepInCrystal9 , sqrt(cSmooth[0] + cSmooth[1]*EdepInCrystal9  + cSmooth[2]*EdepInCrystal9 *EdepInCrystal9 )));
      if(EdepInCrystal10 > 0.) h2->Fill(r3->Gaus(EdepInCrystal10, sqrt(cSmooth[0] + cSmooth[1]*EdepInCrystal10 + cSmooth[2]*EdepInCrystal10*EdepInCrystal10)));
      if(EdepInCrystal11 > 0.) h2->Fill(r3->Gaus(EdepInCrystal11, sqrt(cSmooth[0] + cSmooth[1]*EdepInCrystal11 + cSmooth[2]*EdepInCrystal11*EdepInCrystal11)));
      if(EdepInCrystal12 > 0.) h2->Fill(r3->Gaus(EdepInCrystal12, sqrt(cSmooth[0] + cSmooth[1]*EdepInCrystal12 + cSmooth[2]*EdepInCrystal12*EdepInCrystal12)));
      if(EdepInCrystal13 > 0.) h2->Fill(r3->Gaus(EdepInCrystal13, sqrt(cSmooth[0] + cSmooth[1]*EdepInCrystal13 + cSmooth[2]*EdepInCrystal13*EdepInCrystal13)));
      if(EdepInCrystal14 > 0.) h2->Fill(r3->Gaus(EdepInCrystal14, sqrt(cSmooth[0] + cSmooth[1]*EdepInCrystal14 + cSmooth[2]*EdepInCrystal14*EdepInCrystal14)));
      if(EdepInCrystal15 > 0.) h2->Fill(r3->Gaus(EdepInCrystal15, sqrt(cSmooth[0] + cSmooth[1]*EdepInCrystal15 + cSmooth[2]*EdepInCrystal15*EdepInCrystal15)));
      if(EdepInCrystal16 > 0.) h2->Fill(r3->Gaus(EdepInCrystal16, sqrt(cSmooth[0] + cSmooth[1]*EdepInCrystal16 + cSmooth[2]*EdepInCrystal16*EdepInCrystal16)));
      if(EdepInCrystal17 > 0.) h2->Fill(r3->Gaus(EdepInCrystal17, sqrt(cSmooth[0] + cSmooth[1]*EdepInCrystal17 + cSmooth[2]*EdepInCrystal17*EdepInCrystal17)));
      if(EdepInCrystal18 > 0.) h2->Fill(r3->Gaus(EdepInCrystal18, sqrt(cSmooth[0] + cSmooth[1]*EdepInCrystal18 + cSmooth[2]*EdepInCrystal18*EdepInCrystal18)));
      if(EdepInCrystal19 > 0.) h2->Fill(r3->Gaus(EdepInCrystal19, sqrt(cSmooth[0] + cSmooth[1]*EdepInCrystal19 + cSmooth[2]*EdepInCrystal19*EdepInCrystal19)));
      if(EdepInCrystal20 > 0.) h2->Fill(r3->Gaus(EdepInCrystal20, sqrt(cSmooth[0] + cSmooth[1]*EdepInCrystal20 + cSmooth[2]*EdepInCrystal20*EdepInCrystal20)));
      if(EdepInCrystal21 > 0.) h2->Fill(r3->Gaus(EdepInCrystal21, sqrt(cSmooth[0] + cSmooth[1]*EdepInCrystal21 + cSmooth[2]*EdepInCrystal21*EdepInCrystal21)));
      if(EdepInCrystal22 > 0.) h2->Fill(r3->Gaus(EdepInCrystal22, sqrt(cSmooth[0] + cSmooth[1]*EdepInCrystal22 + cSmooth[2]*EdepInCrystal22*EdepInCrystal22)));
      if(EdepInCrystal23 > 0.) h2->Fill(r3->Gaus(EdepInCrystal23, sqrt(cSmooth[0] + cSmooth[1]*EdepInCrystal23 + cSmooth[2]*EdepInCrystal23*EdepInCrystal23)));
      if(EdepInCrystal24 > 0.) h2->Fill(r3->Gaus(EdepInCrystal24, sqrt(cSmooth[0] + cSmooth[1]*EdepInCrystal24 + cSmooth[2]*EdepInCrystal24*EdepInCrystal24)));
      if(EdepInCrystal25 > 0.) h2->Fill(r3->Gaus(EdepInCrystal25, sqrt(cSmooth[0] + cSmooth[1]*EdepInCrystal25 + cSmooth[2]*EdepInCrystal25*EdepInCrystal25)));
      if(EdepInCrystal26 > 0.) h2->Fill(r3->Gaus(EdepInCrystal26, sqrt(cSmooth[0] + cSmooth[1]*EdepInCrystal26 + cSmooth[2]*EdepInCrystal26*EdepInCrystal26)));
      if(EdepInCrystal27 > 0.) h2->Fill(r3->Gaus(EdepInCrystal27, sqrt(cSmooth[0] + cSmooth[1]*EdepInCrystal27 + cSmooth[2]*EdepInCrystal27*EdepInCrystal27)));
      if(EdepInCrystal28 > 0.) h2->Fill(r3->Gaus(EdepInCrystal28, sqrt(cSmooth[0] + cSmooth[1]*EdepInCrystal28 + cSmooth[2]*EdepInCrystal28*EdepInCrystal28)));
      if(EdepInCrystal29 > 0.) h2->Fill(r3->Gaus(EdepInCrystal29, sqrt(cSmooth[0] + cSmooth[1]*EdepInCrystal29 + cSmooth[2]*EdepInCrystal29*EdepInCrystal29)));
      if(EdepInCrystal30 > 0.) h2->Fill(r3->Gaus(EdepInCrystal30, sqrt(cSmooth[0] + cSmooth[1]*EdepInCrystal30 + cSmooth[2]*EdepInCrystal30*EdepInCrystal30)));

      // if (Cut(ientry) < 0) continue;
   }
   h1->SetLineColor(kBlack);
   h1->SetLineStyle(4);
   h1->Draw("hist");
   h2->Draw("same hist");

   ////////////
   cout << "Writing histograms to file" << endl;

   ofstream myfile;
   myfile.open ("h1.txt");
   myfile << "# E(MeV) counts\n";
   Int_t n = h1->GetNbinsX();
   for (Int_t i=1; i<=n; i++) {
      myfile << h1->GetBinLowEdge(i)+h1->GetBinWidth(i)/2 << "\t"
             << h1->GetBinContent(i)
             << "\n";
   }
   myfile.close();

   myfile.open ("h2.txt");
   myfile << "# E(MeV) counts\n";
   n = h2->GetNbinsX();
   for (Int_t i=1; i<=n; i++) {
      myfile << h2->GetBinLowEdge(i)+h2->GetBinWidth(i)/2 << "\t"
             << h2->GetBinContent(i)
             << "\n";
   }
   myfile.close();

}
