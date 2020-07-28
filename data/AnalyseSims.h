//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat Dec  9 14:34:09 2017 by ROOT version 6.11/01
// from TTree OSCARhits/Energy deposited in the OSCAR detectors
// found on file: ../data/Edep.root
//////////////////////////////////////////////////////////

#ifndef AnalyseSims_h
#define AnalyseSims_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class AnalyseSims {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Double_t        EdepInCrystal1;
   Double_t        EdepInCrystal2;
   Double_t        EdepInCrystal3;
   Double_t        EdepInCrystal4;
   Double_t        EdepInCrystal5;
   Double_t        EdepInCrystal6;
   Double_t        EdepInCrystal7;
   Double_t        EdepInCrystal8;
   Double_t        EdepInCrystal9;
   Double_t        EdepInCrystal10;
   Double_t        EdepInCrystal11;
   Double_t        EdepInCrystal12;
   Double_t        EdepInCrystal13;
   Double_t        EdepInCrystal14;
   Double_t        EdepInCrystal15;
   Double_t        EdepInCrystal16;
   Double_t        EdepInCrystal17;
   Double_t        EdepInCrystal18;
   Double_t        EdepInCrystal19;
   Double_t        EdepInCrystal20;
   Double_t        EdepInCrystal21;
   Double_t        EdepInCrystal22;
   Double_t        EdepInCrystal23;
   Double_t        EdepInCrystal24;
   Double_t        EdepInCrystal25;
   Double_t        EdepInCrystal26;
   Double_t        EdepInCrystal27;
   Double_t        EdepInCrystal28;
   Double_t        EdepInCrystal29;
   Double_t        EdepInCrystal30;

   // List of branches
   TBranch        *b_EdepInCrystal1;   //!
   TBranch        *b_EdepInCrystal2;   //!
   TBranch        *b_EdepInCrystal3;   //!
   TBranch        *b_EdepInCrystal4;   //!
   TBranch        *b_EdepInCrystal5;   //!
   TBranch        *b_EdepInCrystal6;   //!
   TBranch        *b_EdepInCrystal7;   //!
   TBranch        *b_EdepInCrystal8;   //!
   TBranch        *b_EdepInCrystal9;   //!
   TBranch        *b_EdepInCrystal10;   //!
   TBranch        *b_EdepInCrystal11;   //!
   TBranch        *b_EdepInCrystal12;   //!
   TBranch        *b_EdepInCrystal13;   //!
   TBranch        *b_EdepInCrystal14;   //!
   TBranch        *b_EdepInCrystal15;   //!
   TBranch        *b_EdepInCrystal16;   //!
   TBranch        *b_EdepInCrystal17;   //!
   TBranch        *b_EdepInCrystal18;   //!
   TBranch        *b_EdepInCrystal19;   //!
   TBranch        *b_EdepInCrystal20;   //!
   TBranch        *b_EdepInCrystal21;   //!
   TBranch        *b_EdepInCrystal22;   //!
   TBranch        *b_EdepInCrystal23;   //!
   TBranch        *b_EdepInCrystal24;   //!
   TBranch        *b_EdepInCrystal25;   //!
   TBranch        *b_EdepInCrystal26;   //!
   TBranch        *b_EdepInCrystal27;   //!
   TBranch        *b_EdepInCrystal28;   //!
   TBranch        *b_EdepInCrystal29;   //!
   TBranch        *b_EdepInCrystal30;   //!

   AnalyseSims(TTree *tree=0);
   virtual ~AnalyseSims();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef AnalyseSims_cxx
AnalyseSims::AnalyseSims(TTree *tree) : fChain(0)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../data/Edep.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("../data/Edep.root");
      }
      f->GetObject("OSCARhits",tree);

   }
   Init(tree);
}

AnalyseSims::~AnalyseSims()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t AnalyseSims::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t AnalyseSims::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void AnalyseSims::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("EdepInCrystal1", &EdepInCrystal1, &b_EdepInCrystal1);
   fChain->SetBranchAddress("EdepInCrystal2", &EdepInCrystal2, &b_EdepInCrystal2);
   fChain->SetBranchAddress("EdepInCrystal3", &EdepInCrystal3, &b_EdepInCrystal3);
   fChain->SetBranchAddress("EdepInCrystal4", &EdepInCrystal4, &b_EdepInCrystal4);
   fChain->SetBranchAddress("EdepInCrystal5", &EdepInCrystal5, &b_EdepInCrystal5);
   fChain->SetBranchAddress("EdepInCrystal6", &EdepInCrystal6, &b_EdepInCrystal6);
   fChain->SetBranchAddress("EdepInCrystal7", &EdepInCrystal7, &b_EdepInCrystal7);
   fChain->SetBranchAddress("EdepInCrystal8", &EdepInCrystal8, &b_EdepInCrystal8);
   fChain->SetBranchAddress("EdepInCrystal9", &EdepInCrystal9, &b_EdepInCrystal9);
   fChain->SetBranchAddress("EdepInCrystal10", &EdepInCrystal10, &b_EdepInCrystal10);
   fChain->SetBranchAddress("EdepInCrystal11", &EdepInCrystal11, &b_EdepInCrystal11);
   fChain->SetBranchAddress("EdepInCrystal12", &EdepInCrystal12, &b_EdepInCrystal12);
   fChain->SetBranchAddress("EdepInCrystal13", &EdepInCrystal13, &b_EdepInCrystal13);
   fChain->SetBranchAddress("EdepInCrystal14", &EdepInCrystal14, &b_EdepInCrystal14);
   fChain->SetBranchAddress("EdepInCrystal15", &EdepInCrystal15, &b_EdepInCrystal15);
   fChain->SetBranchAddress("EdepInCrystal16", &EdepInCrystal16, &b_EdepInCrystal16);
   fChain->SetBranchAddress("EdepInCrystal17", &EdepInCrystal17, &b_EdepInCrystal17);
   fChain->SetBranchAddress("EdepInCrystal18", &EdepInCrystal18, &b_EdepInCrystal18);
   fChain->SetBranchAddress("EdepInCrystal19", &EdepInCrystal19, &b_EdepInCrystal19);
   fChain->SetBranchAddress("EdepInCrystal20", &EdepInCrystal20, &b_EdepInCrystal20);
   fChain->SetBranchAddress("EdepInCrystal21", &EdepInCrystal21, &b_EdepInCrystal21);
   fChain->SetBranchAddress("EdepInCrystal22", &EdepInCrystal22, &b_EdepInCrystal22);
   fChain->SetBranchAddress("EdepInCrystal23", &EdepInCrystal23, &b_EdepInCrystal23);
   fChain->SetBranchAddress("EdepInCrystal24", &EdepInCrystal24, &b_EdepInCrystal24);
   fChain->SetBranchAddress("EdepInCrystal25", &EdepInCrystal25, &b_EdepInCrystal25);
   fChain->SetBranchAddress("EdepInCrystal26", &EdepInCrystal26, &b_EdepInCrystal26);
   fChain->SetBranchAddress("EdepInCrystal27", &EdepInCrystal27, &b_EdepInCrystal27);
   fChain->SetBranchAddress("EdepInCrystal28", &EdepInCrystal28, &b_EdepInCrystal28);
   fChain->SetBranchAddress("EdepInCrystal29", &EdepInCrystal29, &b_EdepInCrystal29);
   fChain->SetBranchAddress("EdepInCrystal30", &EdepInCrystal30, &b_EdepInCrystal30);
   Notify();
}

Bool_t AnalyseSims::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void AnalyseSims::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t AnalyseSims::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef AnalyseSims_cxx
