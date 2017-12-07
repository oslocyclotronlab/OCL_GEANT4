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
:G4UserEventAction(), EdepInCrystalTest{0.7,9.,9.,9.,9.,9.,9.,9.,9.,9.,9.,9.,9.,9.,9.,9.,9.,9.,9.,9.,9.,9.,9.,9.,9.,9.,9.,9.,9.,9.,9.,9.}, EdepInCrystal20(0.),EdepInCrystal27(0.),EdepInCrystal1(0.),EdepInCrystal2(0.),EdepInCrystal3(0.),EdepInCrystal4(0.),EdepInCrystal5(0.),EdepInCrystal6(0.),EdepInCrystal7(0.),EdepInCrystal8(0.),EdepInCrystal9(0.),EdepInCrystal10(0.),EdepInCrystal11(0.),EdepInCrystal12(0.),EdepInCrystal13(0.),EdepInCrystal14(0.),EdepInCrystal15(0.),EdepInCrystal16(0.),EdepInCrystal17(0.),EdepInCrystal18(0.),EdepInCrystal19(0.),EdepInCrystal21(0.),EdepInCrystal22(0.),EdepInCrystal23(0.),EdepInCrystal24(0.),EdepInCrystal25(0.),EdepInCrystal26(0.),EdepInCrystal28(0.),EdepInCrystal29(0.),EdepInCrystal30(0.),EdepInCrystal31(0.),EdepInCrystal32(0.), nAbsPhotons(0.), absTime(0.), foldedEdep(0.)
{}

SingleScintEventAction::~SingleScintEventAction()
{}

void SingleScintEventAction::BeginOfEventAction(const G4Event*)
{
	// initialisation per event
	EdepInCrystal1 =0.7;
	EdepInCrystal2 =0.;
	EdepInCrystal3 =0.;
	EdepInCrystal4 =0.;
	EdepInCrystal5 =0.;
	EdepInCrystal6 =0.;
	EdepInCrystal7 =0.;
	EdepInCrystal8 =0.;
	EdepInCrystal9 =0.;
	EdepInCrystal10=0.;
	EdepInCrystal11=0.;
	EdepInCrystal12=0.;
	EdepInCrystal13=0.;
	EdepInCrystal14=0.;
	EdepInCrystal15=0.;
	EdepInCrystal16=0.;
	EdepInCrystal17=0.;
	EdepInCrystal18=0.;
	EdepInCrystal19=0.;
	EdepInCrystal20=0.;
	EdepInCrystal21=0.;
	EdepInCrystal22=0.;
	EdepInCrystal23=0.;
	EdepInCrystal24=0.;
	EdepInCrystal25=0.;
	EdepInCrystal26=0.;
	EdepInCrystal27=0.;
	EdepInCrystal28=0.;
	EdepInCrystal29=0.;
	EdepInCrystal30=0.;
	EdepInCrystal31=0.;
	EdepInCrystal32=0.;
	nAbsPhotons = 0;
	absTime = 0;
	foldedEdep=0.;

	for(G4int k=0; k<32; k++) EdepInCrystalTest[k]=0.;
}

void SingleScintEventAction::EndOfEventAction(const G4Event* evt)
{

	  // Accumulate statistics
	  //

	  // get analysis manager
	  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

// fill ntuple

	  // analysisManager->FillNtupleDColumn(0 ,EdepInCrystal1) ;
	  // analysisManager->FillNtupleDColumn(1 ,EdepInCrystal2) ;
	  // analysisManager->FillNtupleDColumn(2 ,EdepInCrystal3) ;
	  // analysisManager->FillNtupleDColumn(3 ,EdepInCrystal4) ;
	  // analysisManager->FillNtupleDColumn(4 ,EdepInCrystal5) ;
	  // analysisManager->FillNtupleDColumn(5 ,EdepInCrystal6) ;
	  // analysisManager->FillNtupleDColumn(6 ,EdepInCrystal7) ;
	  // analysisManager->FillNtupleDColumn(7 ,EdepInCrystal8) ;
	  // analysisManager->FillNtupleDColumn(8 ,EdepInCrystal9) ;
	  // analysisManager->FillNtupleDColumn(9 ,EdepInCrystal10);
	  // analysisManager->FillNtupleDColumn(10,EdepInCrystal11);
	  // analysisManager->FillNtupleDColumn(11,EdepInCrystal12);
	  // analysisManager->FillNtupleDColumn(12,EdepInCrystal13);
	  // analysisManager->FillNtupleDColumn(13,EdepInCrystal14);
	  // analysisManager->FillNtupleDColumn(14,EdepInCrystal15);
	  // analysisManager->FillNtupleDColumn(15,EdepInCrystal16);
	  // analysisManager->FillNtupleDColumn(16,EdepInCrystal17);
	  // analysisManager->FillNtupleDColumn(17,EdepInCrystal18);
	  // analysisManager->FillNtupleDColumn(18,EdepInCrystal19);
	  // analysisManager->FillNtupleDColumn(19,EdepInCrystal20);
	  // analysisManager->FillNtupleDColumn(20,EdepInCrystal21);
	  // analysisManager->FillNtupleDColumn(21,EdepInCrystal22);
	  // analysisManager->FillNtupleDColumn(22,EdepInCrystal23);
	  // analysisManager->FillNtupleDColumn(23,EdepInCrystal24);
	  // analysisManager->FillNtupleDColumn(24,EdepInCrystal25);
	  // analysisManager->FillNtupleDColumn(25,EdepInCrystal26);
	  // analysisManager->FillNtupleDColumn(26,EdepInCrystal27);
	  // analysisManager->FillNtupleDColumn(27,EdepInCrystal28);
	  // analysisManager->FillNtupleDColumn(28,EdepInCrystal29);
	  // analysisManager->FillNtupleDColumn(29,EdepInCrystal30);
	  // analysisManager->FillNtupleDColumn(30,EdepInCrystal31);
	  // analysisManager->FillNtupleDColumn(32,EdepInCrystal32);

	  for(G4int k=0; k<32; k++)analysisManager->FillNtupleDColumn(k,EdepInCrystalTest[k]);

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
