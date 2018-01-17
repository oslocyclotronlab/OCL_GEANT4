
#ifndef OCLRunAction_h
#define OCLRunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

#include <fstream>

using namespace std;

class G4Run;

class OCLRunAction : public G4UserRunAction
{
  public:
    OCLRunAction();
    ~OCLRunAction();

  public:
    void BeginOfRunAction(const G4Run*);
    void EndOfRunAction(const G4Run*);

	G4double xmin;
	G4double xmax;
	G4double binsize;
	G4int nbins;
		
};

#endif

