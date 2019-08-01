
#ifndef OCLRunAction_h
#define OCLRunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

#include <fstream>

using namespace std;

class G4Run;
class OCLRunMessenger;
class G4String;

class OCLRunAction : public G4UserRunAction
{
  public:
    OCLRunAction();
    ~OCLRunAction();

    void SetOutName(G4String nameChoice);
    G4String GetOutName() const { return fOutName; }

  public:
    void BeginOfRunAction(const G4Run*);
    void EndOfRunAction(const G4Run*);

    OCLRunMessenger* fRunMessenger;

    G4String fOutName;

  	G4double xmin;
  	G4double xmax;
  	G4double binsize;
  	G4int nbins;

};

#endif

