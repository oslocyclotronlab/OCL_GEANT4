#include "OCLEventAction.hh"
#include "OCLAnalysis.hh"
//#include "Randomize.hh" // do we really need this?
#include <iomanip>

#include "OCLRunAction.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4UnitsTable.hh"

#include "globals.hh"

//#include "iomanip"

using namespace std;
using namespace CLHEP;

OCLEventAction::OCLEventAction(OCLRunAction* run)
:G4UserEventAction(), EdepInCrystalTest{0.7,9.,9.,9.,9.,9.,9.,9.,9.,9.,9.,9.,9.,9.,9.,9.,9.,9.,9.,9.,9.,9.,9.,9.,9.,9.,9.,9.,9.,9.,9.,9.}, nAbsPhotons(0.), absTime(0.)
{}

OCLEventAction::~OCLEventAction()
{}

void OCLEventAction::BeginOfEventAction(const G4Event*)
{
	// initialisation per event
	nAbsPhotons = 0;
	absTime = 0;

	for(G4int k=0; k<32; k++) EdepInCrystalTest[k]=0.;
}

void OCLEventAction::EndOfEventAction(const G4Event* evt)
{

	  // Accumulate statistics
	  //

	  // get analysis manager
	  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

// fill ntuple
	  //Ignoring the copy numbers belonging to the beam entry and exit:
	  for(G4int k=1; k<31; k++)analysisManager->FillNtupleDColumn(k-1,EdepInCrystalTest[k]);

	  analysisManager->AddNtupleRow();


  // Print per event (modulo n)
  //
  G4int eventID = 1 + evt->GetEventID();
  //G4int printModulo = G4RunManager::GetRunManager()->GetPrintProgress();
  G4int printModulo = 1000;
  if ( ( printModulo > 0 ) && ( eventID % printModulo == 0 ) )
  {
    G4cout << "---> End of event: " << eventID << G4endl;
  }
}
