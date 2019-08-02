#ifndef OCLParallelWorldTargetWheel_h
#define OCLParallelWorldTargetWheel_h 1

#include "globals.hh"
#include "G4VUserParallelWorld.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"

class OCLParallelWorldTargetWheel : public G4VUserParallelWorld
{
  public:
  OCLParallelWorldTargetWheel(G4String worldName);
  virtual ~OCLParallelWorldTargetWheel();

  public:
  virtual void Construct();
  // virtual void ConstructSD();
};

#endif
