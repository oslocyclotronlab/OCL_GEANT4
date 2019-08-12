#include "OCLMaterials.hh"

#include "G4UnitsTable.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4NistManager.hh"

OCLMaterials::OCLMaterials() {
  CreateMaterials();
}

OCLMaterials::~OCLMaterials() {}

OCLMaterials* OCLMaterials::instance = 0;

OCLMaterials* OCLMaterials::GetInstance() {
  if (instance == 0) {
      instance = new OCLMaterials();
    }
  return instance;
}


void OCLMaterials::CreateMaterials() {
  G4String name, symbol;
  G4double a, z;
  G4double pressure, density, temperature, fractionmass;
  G4int nel, ncomponents, natoms;
  G4Material* ma;

  // load NIST material database manager
  G4NistManager * nist = G4NistManager::Instance();

  //---------------------------------
  //Define simple materials
  G4Element* Br = new G4Element("Bromium",    "Br",   z=35.,  a=79.904*g/mole);
  G4Element* Ce = new G4Element("Cerium",     "Cl",   z=58.,  a=140.116*g/mole);
  G4Element* La = new G4Element("Lanthanum",  "La",   z=57.,  a=138.90547*g/mole);
  G4Element* C = new G4Element(name="Carbon", symbol="C", z=6., a = 12.01*g/mole);

  // add more elements from NIST database
  G4Element* Ar  = nist->FindOrBuildElement("Ar");
  G4Element* Cs = nist->FindOrBuildElement("Cs");
  G4Element* K  = nist->FindOrBuildElement("K");
  G4Element* Mg = nist->FindOrBuildElement("Mg");
  G4Element* N  = nist->FindOrBuildElement("N");
  G4Element* O  = nist->FindOrBuildElement("O");
  G4Element* Sb = nist->FindOrBuildElement("Sb");
  G4Element* H = nist->FindOrBuildElement("H");
  G4Element* He = nist->FindOrBuildElement("He");
  G4Element* Ca = nist->FindOrBuildElement("C");

  //
  // define materials from elements.
  //

  G4Material* Al2O3      = nist->FindOrBuildMaterial("G4_ALUMINUM_OXIDE");
  G4Material* Aluminium  = nist->FindOrBuildMaterial("G4_Al");
  G4Material* B2O3       = nist->FindOrBuildMaterial("G4_BORON_OXIDE");
  G4Material* Copper = nist->FindOrBuildMaterial("G4_Cu");
  G4Material* K2O        = nist->FindOrBuildMaterial("G4_POTASSIUM_OXIDE");
  G4Material* lead = nist->FindOrBuildMaterial("G4_Pb");
  G4Material* Na2O       = nist->FindOrBuildMaterial("G4_SODIUM_MONOXIDE");
  G4Material* PlexiGlass = nist->FindOrBuildMaterial("G4_PLEXIGLASS");
  G4Material* Pyrex       = nist->FindOrBuildMaterial("G4_Pyrex_Glass");
  G4Material* Silicon    = nist->FindOrBuildMaterial("G4_Si");
  G4Material* SiO2       = nist->FindOrBuildMaterial("G4_SILICON_DIOXIDE");

  G4Material* Nitrogen = new G4Material(name="N2", density = 1.25053*mg/cm3, ncomponents=1);
  Nitrogen->AddElement(N, 2);

  G4Material* Oxygen = new G4Material(name="O2", density = 1.4289*mg/cm3, ncomponents=1);
  Oxygen->AddElement(O, 2);

  G4Material* Argon = new G4Material(name="Argon", density = 1.7836*mg/cm3, ncomponents=1);
  Argon->AddElement(Ar, 1);


  //
  // composites
  //

  //---------------------------------
  // AIR
  //1 hPa = 1 mbar
  //1 atm = 1013.25 hPa
  // density = P / R.T where R=286.9 J/(Kg.K) [or is it 287.058]
  // standard (haha - whose?) : T = 20C, P = 1013.25 kPa [1013.25 mbar]

  // Dry Air (average composition with Ar), STP
  density = 1.2928*mg/cm3 ;       // STP
  G4Material* Air = new G4Material(name="Air", density, ncomponents=3);
  Air->AddMaterial( Nitrogen, fractionmass = 0.7557 );
  Air->AddMaterial( Oxygen,   fractionmass = 0.2315 );
  Air->AddMaterial( Argon,    fractionmass = 0.0128 );

  // Galactic VACUUM
  density     = CLHEP::universe_mean_density;
  pressure    = 3.e-18*pascal;
  temperature = 2.73 * kelvin;
  G4Material* Vacuum_galactic
    = new G4Material("Galactic", 1, 0.01*g/mole,density,kStateGas,temperature,pressure);

  // ''Normal'' laboratory vacuum
  pressure = 1.e-7*100.*pascal;  // 1e-6 mbar
  density  = 1.2928*mg/cm3 * pressure / (1.01325 * 100.e3*pascal);
  temperature = (273.15+25.) * kelvin;
  G4Material* Vacuum
    = new G4Material(name="Vacuum", density, ncomponents=1,
            kStateGas, temperature, pressure);
  Vacuum->AddMaterial(Air, fractionmass=1.);

  //LaBr3
  G4Material* LaBr3 = new G4Material("LaBr3", density = 5.07*g/cm3, ncomponents=2);
  LaBr3->AddElement(La, natoms=1);
  LaBr3->AddElement(Br, natoms=3);

  //CeBr3
  G4Material* CeBr3 = new G4Material("CeBr3", density = 5.07*g/cm3, ncomponents=2);
  CeBr3->AddElement(Ce, natoms=1);
  CeBr3->AddElement(Br, natoms=3);

  //LaBr3_Ce
  //with 5% dopping, see technical note "BrilLanCe Scintillators Performance Summary"
  // -- didn't fint the doping there any longer, however, adopted the numbers (and doping "method")
  //from http://dx.doi.org/10.1063/1.4810848 now.
  // Potentially it should 5% of the molecules, and not of the weight. However, as CeBr3 has
  // almost the same weight as LaBr3, this shoudln't make a big difference.

  G4Material* LaBr3_Ce = new G4Material("LaBr3_Ce", density = 5.08*g/cm3, ncomponents=2);
  LaBr3_Ce->AddMaterial(LaBr3,  fractionmass=95.*perCent);
  LaBr3_Ce->AddMaterial(CeBr3,   fractionmass=5.*perCent);


  // PMT-materials
  // Borosilicate
  G4Material* Borosilicate = new G4Material("Borosilicate glass",
                                            density= 2.23*g/cm3, ncomponents=5);
  Borosilicate->AddMaterial(SiO2,   fractionmass=80.6 * perCent);
  Borosilicate->AddMaterial(B2O3,  fractionmass=13.0 * perCent);
  Borosilicate->AddMaterial(Na2O,  fractionmass=2.   * perCent); // 1/2 of wt% for (Na20+K20)
  Borosilicate->AddMaterial(K2O,   fractionmass=2.   * perCent); // 1/2 of wt% for (Na20+K20)
  Borosilicate->AddMaterial(Al2O3, fractionmass=2.31  * perCent);

  // Bialkali
  // (Bialkali KCsSb,  Density=?, Thickness=?)?
  G4Material* Bialkali = new G4Material("Bialkali", density= 2.*g/cm3, ncomponents=3);
  Bialkali->AddElement(K,  natoms=2);
  Bialkali->AddElement(Cs, natoms=1);
  Bialkali->AddElement(Sb, natoms=1);

  // MgO reflector
  G4Material* MgO = new G4Material("MgO", density = 3.6*g/cm3, ncomponents=2);
  MgO->AddElement(Mg, natoms=1);
  MgO->AddElement(O, natoms=1);

  //NIFF
  G4Material* isobutane = new G4Material("isoC4H10",density = 2.67*mg/cm3, nel=2);
  isobutane->AddElement(C,4);
  isobutane->AddElement(H,10);

  /*
     iso-butane with conditions specific for NIFF module at OCL
     use conditions: T = 293 K, p = 5 mbar (= 500 Pa)
     To make live easier, we will consider ideal gas law to calculate the density.
  */
  G4double iso_pressure = 500*hep_pascal; //5.e-3*bar;
  G4double iso_Temperature = 293*kelvin;
  G4double iso_molarMass = 58.12*g/mole;
  // G4double gas_cte = CLHEP::k_Boltzmann * CLHEP::Avogadro;
  G4double gas_cte = 8.3144598 * joule / kelvin / mole;
  // G4cout << gas_cte << G4endl;
  G4double iso_density = iso_pressure*iso_molarMass/gas_cte/iso_Temperature; // in g/m3
  // G4cout << iso_density << G4endl;
  G4Material* isobutane_ppac = new G4Material("isoC4H10_PPAC", iso_density,
                                  ncomponents = 1, kStateGas, iso_Temperature, iso_pressure);
  isobutane_ppac->AddMaterial(isobutane, fractionmass = 1.);

  G4Material* Mylar = new G4Material("Mylar", density = 1.39*g/cm3, nel=3);
  Mylar->AddElement(O,2);
  Mylar->AddElement(C,5);
  Mylar->AddElement(H,4);



  // Print all the materials defined.
  // G4cout << G4endl << "The materials defined are : " << G4endl << G4endl;
  // G4cout << *(G4Material::GetMaterialTable()) << G4endl;

}


G4Material* OCLMaterials::GetMaterial(const G4String& name)
{
   // const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();
  G4Material* ma = G4Material::GetMaterial(name);
  if (ma)
  {
    // G4cout << "Material is selected: " << ma->GetName() << G4endl;
    return ma;  // Proceed
  }
  else
  {
    G4cerr << "Material " << name << " does not exist" << G4endl;
    return ma;  // Handle null-pointer error
  }

}


