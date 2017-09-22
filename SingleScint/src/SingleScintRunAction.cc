///////////////////////////////////////////
//
// Oct/2015  Fabio -> RunAction.cc
//
///////////////////////////////////////////

#include "SingleScintRunAction.hh"
#include "SingleScintAnalysis.hh"

#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include "G4Run.hh"
#include "G4ios.hh"

#include <iomanip>

SingleScintRunAction::SingleScintRunAction() : G4UserRunAction()
{
	  // Create analysis manager
	  // The choice of analysis technology is done via selecting of a namespace
	  // in B4Analysis.hh
	  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
	  G4cout << "Using " << analysisManager->GetType() << G4endl;

	  // Create directories
	  //analysisManager->SetHistoDirectoryName("histograms");
	  //analysisManager->SetNtupleDirectoryName("ntuple");
	  analysisManager->SetVerboseLevel(1);
	  analysisManager->SetFirstHistoId(1);

	  // Book histograms, ntuple
	  //

	  // Creating histograms

	  //  G4int CreateH1(const G4String& name, const G4String& title,
	  //                 G4int nbins, G4double xmin, G4double xmax,
	  //                 const G4String& unitName = "none",
	  //                 const G4String& fcnName = "none",
	  //                 const G4String& binSchemeName = "linear");

	  xmin = 0; // in keV
	  xmax = 12e3; // in keV
	  binsize = 2; // in keV
	  nbins= (int)((xmax-xmin)/binsize);
	  analysisManager->CreateH1("Histo1","Edep in Crystal", nbins, xmin*keV, xmax*keV);

	  xmin = 0; //
	  xmax = 12e3; //
	  binsize = 2; //
	  nbins= (int)(xmax-xmin)/binsize;
          analysisManager->CreateH1("Histo2","Absorbed Photons", nbins, xmin, xmax);

      // Here we need some units!
      xmin = 0; // in ns
	  xmax = 500; // in ns
	  binsize = 2; // in ns
	  nbins= (int)(xmax-xmin)/binsize;
          analysisManager->CreateH1("Histo3","Time of Absorption", nbins, xmin, xmax*ns);

	  // Creating ntuple
	  //
	  analysisManager->CreateNtuple("B4", "Edep and TrackL...");
	  analysisManager->CreateNtupleDColumn("Edep");
	  analysisManager->CreateNtupleDColumn("nAbsPhotons");
	  analysisManager->CreateNtupleDColumn("absTime");
	  analysisManager->FinishNtuple();

}

SingleScintRunAction::~SingleScintRunAction()
{
	delete G4AnalysisManager::Instance();
}

void SingleScintRunAction::BeginOfRunAction(const G4Run*)
{
	  //inform the runManager to save random number seed
	  //G4RunManager::GetRunManager()->SetRandomNumberStore(true);

	  // Get analysis manager
	  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

	  // Open an output file
	  //
	  G4String fileName = "../data/Edep.root";
	  analysisManager->OpenFile(fileName);
}

void SingleScintRunAction::EndOfRunAction(const G4Run*)
{
	  // print histogram statistics
	  //
	  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
	  if ( analysisManager->GetH1(1) ) {
	    G4cout << G4endl << " ----> print histograms statistic ";
	    if(isMaster) {
	      G4cout << "for the entire run " << G4endl << G4endl;
	    }
	    else {
	      G4cout << "for the local thread " << G4endl << G4endl;
	    }

	    G4cout << G4endl << " EAbs : mean = "
	       << G4BestUnit(analysisManager->GetH1(1)->mean(), "Energy")
	       << " rms = "
	       << G4BestUnit(analysisManager->GetH1(1)->rms(),  "Energy") << G4endl;

	  }

	  // save histograms & ntuple
	  //
	  analysisManager->Write();
	  analysisManager->CloseFile();

}


