#include "OCLDetectorConstruction.hh"
#include "OCLParallelWorldTargetsOnWheel.hh"
#include "OCLParallelWorldSiRi.hh"
#include "OCLParallelWorldFrameOuter.hh"
#include "OCLParallelWorldTargetChamber.hh"

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
#include "G4ParallelWorldPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

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

  // Detect interactive mode (if no arguments) and define UI session
  //
  G4UIExecutive* ui = 0;
  if ( argc == 1 ) {
    ui = new G4UIExecutive(argc, argv);
  }

  // Construct the default run manager
  //
  // #ifdef G4MULTITHREADED
  //   G4MTRunManager* runManager = new G4MTRunManager;
  // #else
  //   G4RunManager* runManager = new G4RunManager;
  // #endif
  G4RunManager* runManager = new G4RunManager; // Hotfix as long as we have not implemented a Multithreaded version


  // mass world
  //
  OCLDetectorConstruction* massWorld = new OCLDetectorConstruction;

  // // parallel world
  // //
  // // note that the order of initialization is important now:
  // //  "If more than one parallel worlds are defined,
  // //  later-defined world comes on top of others."
  G4String paraWorldName1 = "ParallelWorld Target Chamber";
  massWorld->RegisterParallelWorld(new OCLParallelWorldTargetChamber(paraWorldName1));
  G4String paraWorldName2 = "ParallelWorld Frame Outer";
  massWorld->RegisterParallelWorld(new OCLParallelWorldFrameOuter(paraWorldName2));
  G4String paraWorldName3 = "ParallelWorld SiRi";
  massWorld->RegisterParallelWorld(new OCLParallelWorldSiRi(paraWorldName3));
  G4String paraWorldName4 = "ParallelWorld Targets on Wheel";
  massWorld->RegisterParallelWorld(new OCLParallelWorldTargetsOnWheel(paraWorldName4));

  // // Set mandatory initialization classes
  // // Detector construction
  runManager->SetUserInitialization(massWorld);

  // Physics list
  G4PhysListFactory factory;
  G4VModularPhysicsList* physicsList = 0;
  G4String physName = "QGSP_BIC_HP";
  // reference PhysicsList via its name
  physicsList = factory.GetReferencePhysList(physName);
  physicsList->RegisterPhysics(new G4ParallelWorldPhysics(paraWorldName1, true));
  physicsList->RegisterPhysics(new G4ParallelWorldPhysics(paraWorldName2, true));
  physicsList->RegisterPhysics(new G4ParallelWorldPhysics(paraWorldName3, true));
  physicsList->RegisterPhysics(new G4ParallelWorldPhysics(paraWorldName4, true));
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
  // runManager->Initialize();

  // Initialize visualization
  //
  G4VisManager* visManager = new G4VisExecutive;
  // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
  // G4VisManager* visManager = new G4VisExecutive("Quiet");
  visManager->Initialize();

  // Get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  // Process macro or start UI session
  //
  if ( ! ui ) {
    // batch mode
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command+fileName);
  }
  else {
    // interactive mode
    UImanager->ApplyCommand("/control/execute init_vis.mac");
    ui->SessionStart();
    delete ui;
  }


  // Job termination
  // Free the store: user actions, physics_list and detector_description are
  // owned and deleted by the run manager, so they should not be deleted
  // in the main() program !

  delete visManager;
  delete runManager;

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
