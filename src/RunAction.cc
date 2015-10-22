///////////////////////////////////////////
//
// Oct/2013  E. Nacher -> RunAction.cc
//
///////////////////////////////////////////

#include "RunAction.hh"
#include "Analysis.hh"

#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include "G4Run.hh"
#include "G4ios.hh"

#include <iomanip>

RunAction::RunAction() : G4UserRunAction()
{
	  // Create analysis manager
	  // The choice of analysis technology is done via selectin of a namespace
	  // in B4Analysis.hh
	  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
	  G4cout << "Using " << analysisManager->GetType() << G4endl;

	  // Create directories
	  //analysisManager->SetHistoDirectoryName("histograms");
	  //analysisManager->SetNtupleDirectoryName("ntuple");
	  analysisManager->SetVerboseLevel(1);
	  analysisManager->SetFirstHistoId(1);
	  G4int histID; // Histogram ID staring with SetFirstHistoId(xxx)

	  // Book histograms, ntuple
	  //

	  // Creating histograms
//	  G4String Histoname;
//	  G4String Histotitle
	  G4double xmin = 0; // in keV
	  G4double xmax = 12e3; // in keV
	  G4int binsize = 2; // in keV
	  G4int nbins= (int)(xmax-xmin)/binsize;
	  analysisManager->CreateH1("Histo1","Edep in Crystal", nbins, xmin, xmax*keV);
	  //G4bool SetH1XAxisTitle(G4int id, const G4String& title);

	  analysisManager->SetH1XAxisTitle(histID=1, "Edep (MeV)");

	  // Creating ntuple
	  //
	  analysisManager->CreateNtuple("B4", "Edep and TrackL...");
	  analysisManager->CreateNtupleDColumn("Edep");
}

RunAction::~RunAction()
{
	delete G4AnalysisManager::Instance();
}

void RunAction::BeginOfRunAction(const G4Run*)
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

void RunAction::EndOfRunAction(const G4Run*)
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



