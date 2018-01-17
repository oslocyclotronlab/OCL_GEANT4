
#ifndef OCLEventAction_h
#define OCLEventAction_h 1

#include <fstream>

#include "G4UserEventAction.hh"
#include "globals.hh"

using namespace std;


class G4Event;
class OCLRunAction;

class OCLEventAction : public G4UserEventAction
{
  public:
	OCLEventAction(OCLRunAction*);
    ~OCLEventAction();

  public:
    void BeginOfEventAction(const G4Event*);
    void EndOfEventAction(const G4Event*);


	void accumulateEnergy(G4double);

	G4int        nAbsPhotons;
	//G4double     totEnergyDep;
	G4double     absTime;
	//G4double     totEnergyDepCathod;


	G4double  EdepInCrystalTest[32];
   
  private:
	OCLRunAction*       runAction;

};

#endif
