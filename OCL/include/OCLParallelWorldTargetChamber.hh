#ifndef OCLParallelWorldTargetChamber_h
#define OCLParallelWorldTargetChamber_h 1

#include "globals.hh"
#include "G4VUserParallelWorld.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"

class OCLParallelWorldTargetChamber : public G4VUserParallelWorld
{
  public:
  OCLParallelWorldTargetChamber(G4String worldName);
  virtual ~OCLParallelWorldTargetChamber();

  public:
  virtual void Construct();
  // virtual void ConstructSD();
};

#endif
