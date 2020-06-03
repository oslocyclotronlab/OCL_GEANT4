#ifndef OCLParallelWorldTargetsOnWheel_h
#define OCLParallelWorldTargetsOnWheel_h 1

#include "globals.hh"
#include "G4VUserParallelWorld.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"

class G4GenericMessenger;

class OCLParallelWorldTargetsOnWheel : public G4VUserParallelWorld
{
  public:
  OCLParallelWorldTargetsOnWheel(G4String worldName);
  virtual ~OCLParallelWorldTargetsOnWheel();

  public:
  virtual void Construct();
  // virtual void ConstructSD();

  void DefineCommands();
  G4bool GetUse() { return fuseThisParallelWorld; }

private:
  G4GenericMessenger* fMessenger;
  G4bool fuseThisParallelWorld;
};

#endif
