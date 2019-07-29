#include <iostream>
#include <fstream>
#include <string>

using namespace std;

/***************** Written by Gry in december 2015 *******************
****************** changes made in january 2016
This file generates a macros for runnin GEANT4 simulations with
N gamma rays of a given gamma ray energy Egamma isotropically distributed
between the theta angles thetamin and thetamax.

To generate a set of macros run the bash script generatemacros.sh like this
source generatemacros.sh

To run all the simulations use the script called manyruns.sh like this
source manyruns.sh
******************/
int PrintSimMacs(){
Int_t N = 1e6; //Number of gammas to be simulated
// Double_t mintheta = 160.; //smallest theta angle to be used for the gamma distribution [deg]
// Double_t maxtheta = 180.; //largest theta angle
Double_t Egammalow = 400.; //Lowest gamma energy to be simulated
Double_t Egammastep = 200.; //Steps of gamma energy
int nSteps = 50; // Number of Steps to simulate
// Double_t Egammahigh = 3000.; //Highest gamma energy to be simulated
Double_t Egamma = 0.;
ofstream macfile;
string fout; // basename

for(int k=0;k<nSteps;k++)
 {
fout = "sim" + to_string(k) + ".mac";
macfile.open(fout);

Egamma = Egammalow + k*Egammastep;

cout << "Writing file: " << fout << endl;
cout << "Gamma energy: " << Egamma << endl;
// macfile<<endl;
macfile << "## Particle type, position, energy... \n\
/gps/particle gamma\n\
/gps/number 1\n\
\n\
## Particle source distribution\n\
/gps/pos/type Plane\n\
/gps/pos/shape Ellipse\n\
/gps/pos/centre 0. 0. 0. mm\n\
/gps/pos/halfx 0.75 mm\n\
/gps/pos/halfy 1.25 mm" << endl;
macfile<<endl;
macfile<<"/gps/energy "<<Egamma<<" keV"<<endl;
macfile<<endl;
macfile<<"/gps/ang/type iso"<<endl;
// macfile<<"/gps/ang/mintheta "<< mintheta <<" degree"<<endl;
// macfile<<"/gps/ang/maxtheta "<< maxtheta <<" degree"<<endl;
macfile<<endl;
// macfile<<"/process/inactivate Cerenkov"<<endl;
macfile<<endl;
macfile<<"## Number of events to run"<<endl;
macfile <<"/run/beamOn "<<N<<endl;

macfile.close();
 }
 return 0;
}
