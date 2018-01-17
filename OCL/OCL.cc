

#include "OCLDetectorConstruction.hh"
#include "OCLPhysicsList.hh"
#include "OCLPrimaryGeneratorAction.hh"
#include "OCLRunAction.hh"
#include "OCLEventAction.hh"
#include "OCLSteppingAction.hh"

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "G4UImanager.hh"

#include "G4PhysListFactory.hh"
#include "QGSP_BIC_HP.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

#include "Randomize.hh"
#include "time.h"
#include "globals.hh"

#include "G4ios.hh"
#include "fstream"
#include "iomanip"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv)
{
  // Choose the Random engine
  //
  G4Random::setTheEngine(new CLHEP::RanecuEngine);
//  G4long seed = time(0);
//  G4Random::setTheSeed(seed);
  
  // Construct the default run manager
  //
// #ifdef G4MULTITHREADED
//   G4MTRunManager* runManager = new G4MTRunManager;
// #else
//   G4RunManager* runManager = new G4RunManager;
// #endif
  G4RunManager* runManager = new G4RunManager; // Hotfix as long as we have not implemented a Multithreaded version

  // Set mandatory initialization classes
  // Detector construction
  runManager->SetUserInitialization(new OCLDetectorConstruction());

  // Physics list
  G4PhysListFactory *physListFactory = new G4PhysListFactory();
  G4VUserPhysicsList *physicsList =
            physListFactory->GetReferencePhysList("QGSP_BIC_HP");
  physicsList->SetVerboseLevel(1);
  runManager->SetUserInitialization(physicsList);

  
  // runManager->SetUserInitialization(new OCLPhysicsList);

	// set aditional user action classes
    OCLRunAction* run = new OCLRunAction;
	runManager->SetUserAction(run);
	
	OCLEventAction* event = new OCLEventAction(run);
	runManager->SetUserAction(event);
	
	OCLSteppingAction* step = new OCLSteppingAction(event);
	runManager->SetUserAction(step);
	
	// set mandatory user action class
	runManager->SetUserAction(new OCLPrimaryGeneratorAction);



  // Initialize G4 kernel
  //
  runManager->Initialize();
  
#ifdef G4VIS_USE
  // Initialize visualization
  G4VisManager* visManager = new G4VisExecutive;
  // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
  // G4VisManager* visManager = new G4VisExecutive("Quiet");
  visManager->Initialize();
#endif

  // Get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  if (argc!=1) {
    // batch mode
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command+fileName);
  }
  else {
    // interactive mode : define UI session
#ifdef G4UI_USE
    G4UIExecutive* ui = new G4UIExecutive(argc, argv);
#ifdef G4VIS_USE
    UImanager->ApplyCommand("/control/execute init_vis.mac"); 
#else
    UImanager->ApplyCommand("/control/execute init.mac"); 
#endif
    ui->SessionStart();
    delete ui;
#endif
  }

  // Job termination
  // Free the store: user actions, physics_list and detector_description are
  // owned and deleted by the run manager, so they should not be deleted 
  // in the main() program !
  
#ifdef G4VIS_USE
  delete visManager;
#endif
  delete runManager;

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
