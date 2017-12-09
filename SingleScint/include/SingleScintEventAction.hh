
#ifndef SingleScintEventAction_h
#define SingleScintEventAction_h 1

#include <fstream>

#include "G4UserEventAction.hh"
#include "globals.hh"

using namespace std;


class G4Event;
class SingleScintRunAction;

class SingleScintEventAction : public G4UserEventAction
{
  public:
	SingleScintEventAction(SingleScintRunAction*);
    ~SingleScintEventAction();

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
	SingleScintRunAction*       runAction;

};

#endif
