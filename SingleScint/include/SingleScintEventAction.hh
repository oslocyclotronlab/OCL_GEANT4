
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
	G4double  EdepInCrystal1 ;
	G4double  EdepInCrystal2 ;
	G4double  EdepInCrystal3 ;
	G4double  EdepInCrystal4 ;
	G4double  EdepInCrystal5 ;
	G4double  EdepInCrystal6 ;
	G4double  EdepInCrystal7 ;
	G4double  EdepInCrystal8 ;
	G4double  EdepInCrystal9 ;
	G4double  EdepInCrystal10;
	G4double  EdepInCrystal11;
	G4double  EdepInCrystal12;
	G4double  EdepInCrystal13;
	G4double  EdepInCrystal14;
	G4double  EdepInCrystal15;
	G4double  EdepInCrystal16;
	G4double  EdepInCrystal17;
	G4double  EdepInCrystal18;
	G4double  EdepInCrystal19;
	G4double  EdepInCrystal20;
	G4double  EdepInCrystal21;
	G4double  EdepInCrystal22;
	G4double  EdepInCrystal23;
	G4double  EdepInCrystal24;
	G4double  EdepInCrystal25;
	G4double  EdepInCrystal26;
    G4double  EdepInCrystal27;
	G4double  EdepInCrystal28;
	G4double  EdepInCrystal29;
	G4double  EdepInCrystal30;
	G4double  EdepInCrystal31;
	G4double  EdepInCrystal32;

	G4double foldedEdep;
   
  private:
	SingleScintRunAction*       runAction;

};

#endif
