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
int PrintSimMacs(int k=0){
Int_t N = 5000000; //Number of gammas to be simulated
Double_t mintheta = 160.; //smallest theta angle to be used for the gamma distribution [deg]
Double_t maxtheta = 180.; //largest theta angle
Double_t Egammalow = 500.; //Lowest gamma energy to be simulated
Double_t Egammastep = 1000.; //Steps of gamma energy 
Double_t Egammahigh = 17500.; //Highest gamma energy to be simulated
Double_t Egamma = 0.; 
ofstream macfile;
macfile.open("macfile.txt");

Egamma = Egammalow + k*Egammastep;

macfile<<endl;
macfile<<"## Particle type, position, energy..."<<endl;

macfile<<"/gps/particle gamma"<<endl;
macfile<<"/gps/number 1"<<endl;
macfile<<"# Source position now hardcoded"<<endl;
macfile<<endl;
macfile<<"/gps/energy "<<Egamma<<" keV"<<endl;
macfile<<endl;
//macfile<<"/gps/direction 0 0 1"<<endl;
macfile<<"/gps/ang/type iso"<<endl;
macfile<<"/gps/ang/mintheta "<< mintheta <<" degree"<<endl;
macfile<<"/gps/ang/maxtheta "<< maxtheta <<" degree"<<endl;
macfile<<endl;
macfile<<"/process/inactivate Cerenkov"<<endl;
macfile<<endl;
macfile<<"## Number of events to run"<<endl;
macfile <<"/run/beamOn "<<N<<endl;

macfile.close();
    
}
