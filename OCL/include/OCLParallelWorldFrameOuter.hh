#ifndef OCLParallelWorldFrameOuter_h
#define OCLParallelWorldFrameOuter_h 1

#include "globals.hh"
#include "G4VUserParallelWorld.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

class OCLParallelWorldFrameOuter : public G4VUserParallelWorld
{
  public:
  OCLParallelWorldFrameOuter(G4String worldName);
  virtual ~OCLParallelWorldFrameOuter();

  public:
  virtual void Construct();
  // virtual void ConstructSD();
};

#endif
