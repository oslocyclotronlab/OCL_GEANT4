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

	  // Book ntuple

	  // Creating ntuple
	  //
	  analysisManager->CreateNtuple("OSCARhits", "Energy deposited in the OSCAR detectors");
	  analysisManager->CreateNtupleDColumn("EdepInCrystal1");
	  analysisManager->CreateNtupleDColumn("EdepInCrystal2");
	  analysisManager->CreateNtupleDColumn("EdepInCrystal3");
	  analysisManager->CreateNtupleDColumn("EdepInCrystal4");
	  analysisManager->CreateNtupleDColumn("EdepInCrystal5");
	  analysisManager->CreateNtupleDColumn("EdepInCrystal6");
	  analysisManager->CreateNtupleDColumn("EdepInCrystal7");
	  analysisManager->CreateNtupleDColumn("EdepInCrystal8");
	  analysisManager->CreateNtupleDColumn("EdepInCrystal9");
	  analysisManager->CreateNtupleDColumn("EdepInCrystal10");
	  analysisManager->CreateNtupleDColumn("EdepInCrystal11");
	  analysisManager->CreateNtupleDColumn("EdepInCrystal12");
	  analysisManager->CreateNtupleDColumn("EdepInCrystal13");
	  analysisManager->CreateNtupleDColumn("EdepInCrystal14");
	  analysisManager->CreateNtupleDColumn("EdepInCrystal15");
	  analysisManager->CreateNtupleDColumn("EdepInCrystal16");
	  analysisManager->CreateNtupleDColumn("EdepInCrystal17");
	  analysisManager->CreateNtupleDColumn("EdepInCrystal18");
	  analysisManager->CreateNtupleDColumn("EdepInCrystal19");
	  analysisManager->CreateNtupleDColumn("EdepInCrystal20");
	  analysisManager->CreateNtupleDColumn("EdepInCrystal21");
	  analysisManager->CreateNtupleDColumn("EdepInCrystal22");
	  analysisManager->CreateNtupleDColumn("EdepInCrystal23");
	  analysisManager->CreateNtupleDColumn("EdepInCrystal24");
	  analysisManager->CreateNtupleDColumn("EdepInCrystal25");
	  analysisManager->CreateNtupleDColumn("EdepInCrystal26");
	  analysisManager->CreateNtupleDColumn("EdepInCrystal27");
	  analysisManager->CreateNtupleDColumn("EdepInCrystal28");
	  analysisManager->CreateNtupleDColumn("EdepInCrystal29");
	  analysisManager->CreateNtupleDColumn("EdepInCrystal30");
	  analysisManager->CreateNtupleDColumn("EdepInCrystal31");
	  analysisManager->CreateNtupleDColumn("EdepInCrystal32");

	  analysisManager->CreateNtupleDColumn("EdepInCrystal32", "EdepInCrystal32");

	  // analysisManager->CreateNtupleDColumn("EdepInCrystal2","EdepInCrystal1");

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


