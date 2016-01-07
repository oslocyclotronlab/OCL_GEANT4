#include "SingleScintEventAction.hh"
#include "SingleScintAnalysis.hh"
//#include "Randomize.hh" // do we really need this?
#include <iomanip>

#include "SingleScintRunAction.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4UnitsTable.hh"

#include "globals.hh"

//#include "iomanip"

using namespace std;
using namespace CLHEP;

SingleScintEventAction::SingleScintEventAction(SingleScintRunAction* run)
:G4UserEventAction(), EdepInCrystal(0.), nAbsPhotons(0.), absTime(0.)
{}

SingleScintEventAction::~SingleScintEventAction()
{}

void SingleScintEventAction::BeginOfEventAction(const G4Event*)
{
	// initialisation per event
	EdepInCrystal = 0.;
	nAbsPhotons = 0;
	absTime = 0;
}

void SingleScintEventAction::EndOfEventAction(const G4Event* evt)
{

	  // Accumulate statistics
	  //

	  // get analysis manager
	  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

	  // fill histograms
	  analysisManager->FillH1(1, EdepInCrystal);
	  analysisManager->FillH1(2, nAbsPhotons);
	  analysisManager->FillH1(3, absTime);

	  // fill ntuple
	  analysisManager->FillNtupleDColumn(0, EdepInCrystal);
	  analysisManager->FillNtupleDColumn(1, nAbsPhotons);
	  analysisManager->FillNtupleDColumn(2, absTime);
	  analysisManager->AddNtupleRow();

  // Print per event (modulo n)
  //
  G4int eventID = 1 + evt->GetEventID();
  //G4int printModulo = G4RunManager::GetRunManager()->GetPrintProgress();
  G4int printModulo = 100;
  if ( ( printModulo > 0 ) && ( eventID % printModulo == 0 ) )
  {
    G4cout << "---> End of event: " << eventID << G4endl;
  }
}
