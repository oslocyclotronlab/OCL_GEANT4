
#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"

class G4VPhysicalVolume;

class TestDetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    TestDetectorConstruction();
    ~TestDetectorConstruction();

  public:
     G4VPhysicalVolume* Construct();

};

#endif

