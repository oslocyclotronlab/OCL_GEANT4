///////////////////////////////////////////
//
// Oct/2015  Fabio -> RunAction.cc
//
///////////////////////////////////////////

#include "OCLRunAction.hh"
#include "OCLAnalysis.hh"

#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include "G4Run.hh"
#include "G4ios.hh"

#include <iomanip>

OCLRunAction::OCLRunAction() : G4UserRunAction()
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

	  analysisManager->FinishNtuple();


}

OCLRunAction::~OCLRunAction()
{
	delete G4AnalysisManager::Instance();
}

void OCLRunAction::BeginOfRunAction(const G4Run*)
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

void OCLRunAction::EndOfRunAction(const G4Run*)
{
	  // save ntuples
	  //
	  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

	  analysisManager->Write();
	  analysisManager->CloseFile();

}


