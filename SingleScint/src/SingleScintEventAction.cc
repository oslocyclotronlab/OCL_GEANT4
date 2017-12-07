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
:G4UserEventAction(), EdepInCrystal20(0.),EdepInCrystal27(0.),EdepInCrystal1(0.),EdepInCrystal2(0.),EdepInCrystal3(0.),EdepInCrystal4(0.),EdepInCrystal5(0.),EdepInCrystal6(0.),EdepInCrystal7(0.),EdepInCrystal8(0.),EdepInCrystal9(0.),EdepInCrystal10(0.),EdepInCrystal11(0.),EdepInCrystal12(0.),EdepInCrystal13(0.),EdepInCrystal14(0.),EdepInCrystal15(0.),EdepInCrystal16(0.),EdepInCrystal17(0.),EdepInCrystal18(0.),EdepInCrystal19(0.),EdepInCrystal21(0.),EdepInCrystal22(0.),EdepInCrystal23(0.),EdepInCrystal24(0.),EdepInCrystal25(0.),EdepInCrystal26(0.),EdepInCrystal28(0.),EdepInCrystal29(0.),EdepInCrystal30(0.),EdepInCrystal31(0.),EdepInCrystal32(0.), nAbsPhotons(0.), absTime(0.), foldedEdep(0.)
{}

SingleScintEventAction::~SingleScintEventAction()
{}

void SingleScintEventAction::BeginOfEventAction(const G4Event*)
{
	// initialisation per event
	EdepInCrystal1 =0.;
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
}

void SingleScintEventAction::EndOfEventAction(const G4Event* evt)
{

	  // Accumulate statistics
	  //

	  // get analysis manager
	  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

	  // fill histograms
	
	  // analysisManager->FillH1(1 , EdepInCrystal1);
	  // analysisManager->FillH1(2 , EdepInCrystal2);
	  // analysisManager->FillH1(3 , EdepInCrystal3);
	  // analysisManager->FillH1(4 , EdepInCrystal4);
	  // analysisManager->FillH1(5 , EdepInCrystal5);
	  // analysisManager->FillH1(6 , EdepInCrystal6);
	  // analysisManager->FillH1(7 , EdepInCrystal7);
	  // analysisManager->FillH1(8 , EdepInCrystal8);
	  // analysisManager->FillH1(9 , EdepInCrystal9);
	  // analysisManager->FillH1(10,EdepInCrystal10);
	  // analysisManager->FillH1(11,EdepInCrystal11);
	  // analysisManager->FillH1(12,EdepInCrystal12);
	  // analysisManager->FillH1(13,EdepInCrystal13);
	  // analysisManager->FillH1(14,EdepInCrystal14);
	  // analysisManager->FillH1(15,EdepInCrystal15);
	  // analysisManager->FillH1(16,EdepInCrystal16);
	  // analysisManager->FillH1(17,EdepInCrystal17);
	  // analysisManager->FillH1(18,EdepInCrystal18);
	  // analysisManager->FillH1(19,EdepInCrystal19);
	  // analysisManager->FillH1(20,EdepInCrystal20);
	  // analysisManager->FillH1(21,EdepInCrystal21);
	  // analysisManager->FillH1(22,EdepInCrystal22);
	  // analysisManager->FillH1(23,EdepInCrystal23);
	  // analysisManager->FillH1(24,EdepInCrystal24);
	  // analysisManager->FillH1(25,EdepInCrystal25);
	  // analysisManager->FillH1(26,EdepInCrystal26);
	  // analysisManager->FillH1(27,EdepInCrystal27);
	  // analysisManager->FillH1(28,EdepInCrystal28);
	  // analysisManager->FillH1(29,EdepInCrystal29);
	  // analysisManager->FillH1(30,EdepInCrystal30);
	  // analysisManager->FillH1(31,EdepInCrystal31);
	  // analysisManager->FillH1(32,EdepInCrystal32);

   //    G4double foldFactor = 0.0185;
	  // analysisManager->FillH1(33, (G4RandGauss::shoot(EdepInCrystal1 ,foldFactor*sqrt(EdepInCrystal1 ))) );
	  // analysisManager->FillH1(34, (G4RandGauss::shoot(EdepInCrystal2 ,foldFactor*sqrt(EdepInCrystal2 ))) );
	  // analysisManager->FillH1(35, (G4RandGauss::shoot(EdepInCrystal3 ,foldFactor*sqrt(EdepInCrystal3 ))) );
	  // analysisManager->FillH1(36, (G4RandGauss::shoot(EdepInCrystal4 ,foldFactor*sqrt(EdepInCrystal4 ))) );
	  // analysisManager->FillH1(37, (G4RandGauss::shoot(EdepInCrystal5 ,foldFactor*sqrt(EdepInCrystal5 ))) );
	  // analysisManager->FillH1(38, (G4RandGauss::shoot(EdepInCrystal6 ,foldFactor*sqrt(EdepInCrystal6 ))) );
	  // analysisManager->FillH1(39, (G4RandGauss::shoot(EdepInCrystal7 ,foldFactor*sqrt(EdepInCrystal7 ))) );
	  // analysisManager->FillH1(40, (G4RandGauss::shoot(EdepInCrystal8 ,foldFactor*sqrt(EdepInCrystal8 ))) );
	  // analysisManager->FillH1(41, (G4RandGauss::shoot(EdepInCrystal9 ,foldFactor*sqrt(EdepInCrystal9 ))) );
	  // analysisManager->FillH1(42, (G4RandGauss::shoot(EdepInCrystal10,foldFactor*sqrt(EdepInCrystal10))) );
	  // analysisManager->FillH1(43, (G4RandGauss::shoot(EdepInCrystal11,foldFactor*sqrt(EdepInCrystal11))) );
	  // analysisManager->FillH1(44, (G4RandGauss::shoot(EdepInCrystal12,foldFactor*sqrt(EdepInCrystal12))) );
	  // analysisManager->FillH1(45, (G4RandGauss::shoot(EdepInCrystal13,foldFactor*sqrt(EdepInCrystal13))) );
	  // analysisManager->FillH1(46, (G4RandGauss::shoot(EdepInCrystal14,foldFactor*sqrt(EdepInCrystal14))) );
	  // analysisManager->FillH1(47, (G4RandGauss::shoot(EdepInCrystal15,foldFactor*sqrt(EdepInCrystal15))) );
	  // analysisManager->FillH1(48, (G4RandGauss::shoot(EdepInCrystal16,foldFactor*sqrt(EdepInCrystal16))) );
	  // analysisManager->FillH1(49, (G4RandGauss::shoot(EdepInCrystal17,foldFactor*sqrt(EdepInCrystal17))) );
	  // analysisManager->FillH1(50, (G4RandGauss::shoot(EdepInCrystal18,foldFactor*sqrt(EdepInCrystal18))) );
	  // analysisManager->FillH1(51, (G4RandGauss::shoot(EdepInCrystal19,foldFactor*sqrt(EdepInCrystal19))) );
	  // analysisManager->FillH1(52, (G4RandGauss::shoot(EdepInCrystal20,foldFactor*sqrt(EdepInCrystal20))) );
	  // analysisManager->FillH1(53, (G4RandGauss::shoot(EdepInCrystal21,foldFactor*sqrt(EdepInCrystal21))) );
	  // analysisManager->FillH1(54, (G4RandGauss::shoot(EdepInCrystal22,foldFactor*sqrt(EdepInCrystal22))) );
	  // analysisManager->FillH1(55, (G4RandGauss::shoot(EdepInCrystal23,foldFactor*sqrt(EdepInCrystal23))) );
	  // analysisManager->FillH1(56, (G4RandGauss::shoot(EdepInCrystal24,foldFactor*sqrt(EdepInCrystal24))) );
	  // analysisManager->FillH1(57, (G4RandGauss::shoot(EdepInCrystal25,foldFactor*sqrt(EdepInCrystal25))) );
	  // analysisManager->FillH1(58, (G4RandGauss::shoot(EdepInCrystal26,foldFactor*sqrt(EdepInCrystal26))) );
	  // analysisManager->FillH1(59, (G4RandGauss::shoot(EdepInCrystal27,foldFactor*sqrt(EdepInCrystal27))) );
	  // analysisManager->FillH1(60, (G4RandGauss::shoot(EdepInCrystal28,foldFactor*sqrt(EdepInCrystal28))) );
	  // analysisManager->FillH1(61, (G4RandGauss::shoot(EdepInCrystal29,foldFactor*sqrt(EdepInCrystal29))) );
	  // analysisManager->FillH1(62, (G4RandGauss::shoot(EdepInCrystal30,foldFactor*sqrt(EdepInCrystal30))) );
	  // analysisManager->FillH1(63, (G4RandGauss::shoot(EdepInCrystal31,foldFactor*sqrt(EdepInCrystal31))) );
	  // analysisManager->FillH1(64, (G4RandGauss::shoot(EdepInCrystal32,foldFactor*sqrt(EdepInCrystal32))) );

// fill ntuple

	  analysisManager->FillNtupleDColumn(0 ,EdepInCrystal1) ;
	  analysisManager->FillNtupleDColumn(1 ,EdepInCrystal2) ;
	  analysisManager->FillNtupleDColumn(2 ,EdepInCrystal3) ;
	  analysisManager->FillNtupleDColumn(3 ,EdepInCrystal4) ;
	  analysisManager->FillNtupleDColumn(4 ,EdepInCrystal5) ;
	  analysisManager->FillNtupleDColumn(5 ,EdepInCrystal6) ;
	  analysisManager->FillNtupleDColumn(6 ,EdepInCrystal7) ;
	  analysisManager->FillNtupleDColumn(7 ,EdepInCrystal8) ;
	  analysisManager->FillNtupleDColumn(8 ,EdepInCrystal9) ;
	  analysisManager->FillNtupleDColumn(9 ,EdepInCrystal10);
	  analysisManager->FillNtupleDColumn(10,EdepInCrystal11);
	  analysisManager->FillNtupleDColumn(11,EdepInCrystal12);
	  analysisManager->FillNtupleDColumn(12,EdepInCrystal13);
	  analysisManager->FillNtupleDColumn(13,EdepInCrystal14);
	  analysisManager->FillNtupleDColumn(14,EdepInCrystal15);
	  analysisManager->FillNtupleDColumn(15,EdepInCrystal16);
	  analysisManager->FillNtupleDColumn(16,EdepInCrystal17);
	  analysisManager->FillNtupleDColumn(17,EdepInCrystal18);
	  analysisManager->FillNtupleDColumn(18,EdepInCrystal19);
	  analysisManager->FillNtupleDColumn(19,EdepInCrystal20);
	  analysisManager->FillNtupleDColumn(20,EdepInCrystal21);
	  analysisManager->FillNtupleDColumn(21,EdepInCrystal22);
	  analysisManager->FillNtupleDColumn(22,EdepInCrystal23);
	  analysisManager->FillNtupleDColumn(23,EdepInCrystal24);
	  analysisManager->FillNtupleDColumn(24,EdepInCrystal25);
	  analysisManager->FillNtupleDColumn(25,EdepInCrystal26);
	  analysisManager->FillNtupleDColumn(26,EdepInCrystal27);
	  analysisManager->FillNtupleDColumn(27,EdepInCrystal28);
	  analysisManager->FillNtupleDColumn(28,EdepInCrystal29);
	  analysisManager->FillNtupleDColumn(29,EdepInCrystal30);
	  analysisManager->FillNtupleDColumn(30,EdepInCrystal31);
	  analysisManager->FillNtupleDColumn(31,EdepInCrystal32);

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
