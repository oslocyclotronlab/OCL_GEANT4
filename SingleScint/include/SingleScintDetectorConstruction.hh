
#ifndef SingleScintDetectorConstruction_h
#define SingleScintDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

/// Detector construction class to define materials and geometry.
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4PVPlacement;
class G4Box;

class SingleScintDetectorConstruction : public G4VUserDetectorConstruction
{
  public:

    SingleScintDetectorConstruction();
    ~SingleScintDetectorConstruction();

  public:
    virtual G4VPhysicalVolume* Construct();
    //virtual void ConstructSDandField();

  private:
    
    G4Box* solidWorld;
    G4LogicalVolume* WorldLog;
    G4PVPlacement* WorldPhys;
    
};

#endif

