#include "OCLParallelWorldTargetChamber.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Box.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4GDMLParser.hh"
// #include "G4GDMLReadSetup.hh"
#include "G4PhysicalVolumeStore.hh"

OCLParallelWorldTargetChamber::OCLParallelWorldTargetChamber(G4String worldName)
:G4VUserParallelWorld(worldName)
{;}

OCLParallelWorldTargetChamber::~OCLParallelWorldTargetChamber()
{;}

void OCLParallelWorldTargetChamber::Construct()
{
  G4bool useThisParallelWorld = true;

  if (useThisParallelWorld) {
    // Get mass world
    G4VPhysicalVolume* ghostWorld = GetWorld();
    G4LogicalVolume* worldLogical = ghostWorld->GetLogicalVolume();

    // place volumes in the parallel world here.
    // read from CAD
    G4GDMLParser parser;
    // bool flag for reading is for validation -> gives lots of unnecessary warnings
    parser.Read("../OCL/Mesh-Models/target_chamber.gdml", false);
    // remember: each "Wold" volume read with the CAD parser
    // must have a unique name (corresponding to the one used in the gdml file)
    G4LogicalVolume* CADWorldLog = parser.GetVolume("World_target_chamber");

    // set the CAD world volume to
    CADWorldLog->SetMaterial(nullptr);

    // Colors
    G4VisAttributes* VisAtt = new G4VisAttributes(G4Colour(0.85,0.91,0.97, 0.5));

    CADWorldLog->SetVisAttributes(VisAtt);
    for (int i=0; i < CADWorldLog->GetNoDaughters(); i++){
      G4LogicalVolume* daughterLog = CADWorldLog->GetDaughter(i)->GetLogicalVolume();

      // colors for the plexi glass
      if (daughterLog->GetName() == "Target_Ball_Plexi_Lower"
          || daughterLog->GetName() == "Target_Ball_Plexi_Upper"){
        daughterLog->SetVisAttributes(G4VisAttributes(G4Colour(0.80,0.90,0.90, 0.1)));
      }
      else { daughterLog->SetVisAttributes(VisAtt);}
    }

    new G4PVPlacement(0, G4ThreeVector(), CADWorldLog,
                      "ParallelWorld Target Chamber", worldLogical, 0, 0);
    }
  else {} // Do nothing
}
