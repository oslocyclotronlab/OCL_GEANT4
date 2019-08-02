#ifndef OCLParallelWorldSiRi_h
#define OCLParallelWorldSiRi_h 1

#include "globals.hh"
#include "G4VUserParallelWorld.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"

class OCLParallelWorldSiRi : public G4VUserParallelWorld
{
  public:
  OCLParallelWorldSiRi(G4String worldName);
  virtual ~OCLParallelWorldSiRi();

  public:
  virtual void Construct();
  // virtual void ConstructSD();
};

#endif
