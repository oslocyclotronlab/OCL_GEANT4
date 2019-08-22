#include "OCLParallelWorldFrameOuter.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Box.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4GDMLParser.hh"
#include "G4GenericMessenger.hh"


OCLParallelWorldFrameOuter::OCLParallelWorldFrameOuter(G4String worldName)
: G4VUserParallelWorld(worldName),
  fMessenger(nullptr),
  fuseThisParallelWorld(true)
{
  DefineCommands();
}

OCLParallelWorldFrameOuter::~OCLParallelWorldFrameOuter()
{
  delete fMessenger;
}

void OCLParallelWorldFrameOuter::Construct()
{
  G4cout << "Use Parallel World CADFrameOuter: " << fuseThisParallelWorld << G4endl;

  if (fuseThisParallelWorld) {
    // Get mass world
    G4VPhysicalVolume* ghostWorld = GetWorld();
    G4LogicalVolume* worldLogical = ghostWorld->GetLogicalVolume();

    // place volumes in the parallel world here.
    // read from CAD
    G4GDMLParser parser;
    // bool flag for reading is for validation -> gives lots of unnecessary warnings
    parser.Read("../OCL/Mesh-Models/frame_outer.gdml", false);
    // remember: each "Wold" volume read with the CAD parser
    // must have a unique name (corresponding to the one used in the gdml file)
    G4LogicalVolume* CADWorldLog = parser.GetVolume("World_frame_outer");

    // set the CAD world volume to
    CADWorldLog->SetMaterial(nullptr);

    // Colors
    G4VisAttributes* VisAtt = new G4VisAttributes(G4Colour(0.95,0.95,0.95, 0.1));
    CADWorldLog->SetVisAttributes(VisAtt);
    for (int i=0; i < CADWorldLog->GetNoDaughters(); i++){
      CADWorldLog->GetDaughter(i)->GetLogicalVolume ()->SetVisAttributes(VisAtt);
    }


    new G4PVPlacement(0, G4ThreeVector(), CADWorldLog,
                      "ParallelWorld Frame Outer", worldLogical, 0, 0);
    }
  else {} // Do nothing
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void OCLParallelWorldFrameOuter::DefineCommands()
{
  // Define command directory using generic messenger class
  fMessenger = new G4GenericMessenger(this,
                                      "/OCL/det/",
                                      "Detector control");

  // define the command
  auto& useCmd
    = fMessenger->DeclareProperty("useCADFrameOuter",
                                  fuseThisParallelWorld,
                                  "Set usage of CAD Frame");
  useCmd.SetStates(G4State_PreInit);
  useCmd.SetParameterName("bool", true);
  useCmd.SetDefaultValue("true");
  useCmd.SetGuidance("Due to the useage of the _HP physics list, we have to set this command before initialization");
}
