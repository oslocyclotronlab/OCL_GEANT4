#ifndef OCLMaterials_h
#define OCLMaterials_h 1

#include "globals.hh"
#include "G4SystemOfUnits.hh"

class G4Material;

class OCLMaterials
{
  public:
    //OCLMaterials();
    ~OCLMaterials();
    static OCLMaterials* GetInstance();

    G4Material* GetMaterial(const G4String&);


  private:
    OCLMaterials();
    void CreateMaterials();

    static OCLMaterials* instance;

};

#endif



