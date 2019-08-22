#ifndef OCLParallelWorldFrameOuter_h
#define OCLParallelWorldFrameOuter_h 1

#include "globals.hh"
#include "G4VUserParallelWorld.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

class G4GenericMessenger;

class OCLParallelWorldFrameOuter : public G4VUserParallelWorld
{
  public:
  OCLParallelWorldFrameOuter(G4String worldName);
  virtual ~OCLParallelWorldFrameOuter();

  public:
  virtual void Construct();

  void DefineCommands();
  G4bool GetUse() { return fuseThisParallelWorld; }

  private:
  G4GenericMessenger* fMessenger;
  G4bool fuseThisParallelWorld;
};

#endif
