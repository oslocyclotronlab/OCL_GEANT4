// Retreives the number of counts in the peaks and background of the spectra 
// and writes this to file
// additional files to be read can be inserted via the /fnames.insert/ function

#include <map>
#include <sstream>
#include <stdexcept>
#include <string>

#include <TROOT.h>
// #include <TSystem.h>
#include <iostream>
#include <cstdio>
#include "TFile.h"

#include "AnalyseSims.C"
#ifdef MAKECINT
#pragma link C++ defined_in “AnalyseSims.h”;
#pragma link C++ defined_in “AnalyseSims.C”;
#endif

#include <vector>
using namespace std;

////////////////////////////
// Parameters
double threshold = 0.3; // Detector threshold in MeV
////////////////////////////


typedef std::map<std::string, std::string> ConfigInfo;

ConfigInfo GetConfigValues(std::string fname){
  std::ifstream fileStream(fname);
      ConfigInfo configValues;
      std::string line;
          while (std::getline(fileStream, line))
          {
              std::istringstream is_line(line);
              std::string key;
              if (std::getline(is_line, key, ' '/*delimiter*/))
              {
                  std::string value;
                  if (key[0] == '#')
                      continue;

                  if (std::getline(is_line, value))
                  {
                      configValues[key] = value;
                      // cout << value << endl;
                  }
              }
          }
    return configValues;
}

double GetValueFromConfigs(ConfigInfo configValues, string searchKey, int checkPos, string checkValue){
  std::string token;
  std::string keyVal;

  // find the correct key
  for (auto E: configValues) 
  { 
      if (E.first==searchKey)
      {
        keyVal = E.second;
      } 
  }; 

  // split the key to vector with all strings
  std::istringstream iss(keyVal);
  std::vector<std::string> keyVal_split(std::istream_iterator<std::string>{iss},
                                   std::istream_iterator<std::string>());

  // // Iterate and print values of vector
  // for(auto i : results) {
  //     std::cout << i << '\n';
  //     }

  if(checkPos<0){
    return std::stod(keyVal_split[0]);
  }
  else if(keyVal_split[checkPos]==checkValue){
    return std::stod(keyVal_split[0]);
  }
  else{
    cout << "Didn't match checkValue/units" << endl;
    exit(1);
  }
}


double GetBinFromXval(TH1D* h1, double xValueofBin)
{
  TAxis *xaxis = h1->GetXaxis();
  TAxis *yaxis = h1->GetYaxis();
  Int_t binx = xaxis->FindBin(xValueofBin);
  double binvalue = h1->GetBinContent(binx);
  return binvalue;
}


double GetBinWoBG(TH1D* h1,double xValueofBin ,Int_t nChannelsBG)
{
  TAxis *xaxis = h1->GetXaxis();
  TAxis *yaxis = h1->GetYaxis();
  Int_t binx = xaxis->FindBin(xValueofBin);
  double binvalue = h1->GetBinContent(binx);


  double cntsInPeak = binvalue;
  cntsInPeak += h1->GetBinContent(binx+1) + h1->GetBinContent(binx-1);

  // get BG
  Int_t min = binx-nChannelsBG;
  Int_t max = binx + nChannelsBG;
  cout << "Subtract BG from " <<  h1->GetXaxis()->GetBinCenter(min)
       << "to " <<  h1->GetXaxis()->GetBinCenter(max) << endl;
  double sum = 0;
  for(int ch=min;ch<max+1;ch++ ){
    sum += h1->GetBinContent(ch);
  }
  // a little messy -- this accounts also for thaat counts can come in the bin below/above
  double cntbg = (sum - cntsInPeak) /((2.*nChannelsBG)-3);
  double cntWobg = cntsInPeak - cntbg*3.;
  return cntWobg;

}





/////////////////////////////////////////
/////////////////////////////////////////

void GetPeaks(){


vector<double> Eg;
vector<double> nEvents;
vector<double> nCounts;
vector<double> cntFE ;
vector<double> cntSE ;
vector<double> cntDE ;
vector<double> cnt511;
vector<double> cntCompt;
Int_t nChannelsBG;// # of channels, symmetric, to choose to average out the BG around the peaks
int xstart; // dummy for startbin
int xstop; // dummy for stop bin
double sum; // dummy for a sum

std::vector<string> fnames;

string fout_name = "Peaks.dat";

//////////////////////

// get list of files matching certain criterium

string fname_tmp;
system("xterm -e 'find sim*.mac > tmp.txt'");
ifstream tmpfile("tmp.txt", ios::in);
while(tmpfile>>fname_tmp)
{
    fnames.push_back(fname_tmp); //adding data in to the vector
}

// for(auto n : fnames) {
//   cout << n << endl;
// }

// additional files to be analyzed
// fnames.insert(fnames.end(), { "myfile1.mac", "myfile2.mac" });
// fnames.insert(fnames.end(), { "sim10.mac"});

for(auto fname : fnames)
{
  // get config value from file
  // std::string fname = "sim2.mac";
  cout << "Working on " << fname << endl;
  auto configValues=GetConfigValues(fname);

  // get gamma-ray energy from config values
  int checkPos;
  std::string searchKey = "/gps/energy";
  checkPos = 1;
  string checkValue="keV";
  double EgFE = GetValueFromConfigs(configValues, searchKey, checkPos, checkValue);
  cout << "Full energy peak at:" << EgFE << " " << checkValue << endl;
  EgFE /=1000.; //  keV to MeV
  Eg.push_back(EgFE);

  // read in root file
  auto fname_rootStr = fname + ".root";
  const char *fname_root = fname_rootStr.c_str();
  TFile *file_root = new TFile(fname_root, "READ");

  TTree *OSCARhits;
  file_root->GetObject("OSCARhits",OSCARhits);

  // execute the macro that combines all detectors
  AnalyseSims t;
  t.Init(OSCARhits);
  t.Loop();

  // get the histogram from AnalyseSims

  TH1D *h1 = (TH1D*)gDirectory->GetList()->FindObject("h1");
  h1->Draw();
  //

  // cntFE = GetBinFromXval(h1, EgFE);
  cntFE.push_back(GetBinWoBG(h1,EgFE,nChannelsBG));

  // get peaks on top of bg -- 
  nChannelsBG = 4;
  cntSE.push_back(GetBinWoBG(h1,EgFE-0.511,nChannelsBG));
  cntDE.push_back(GetBinWoBG(h1,EgFE-2.*0.511,nChannelsBG));
  cnt511.push_back(GetBinWoBG(h1,0.511,nChannelsBG));
  if(EgFE<1.22){
    cntSE.back()=0;
    cntDE.back()=0;
    cnt511.back()=0;
  }

  nCounts.push_back(h1->GetEntries());

  // assume that compton is all but the peaks; staring at threshold
  xstart = h1->GetXaxis()->FindBin(threshold);
  sum = 0;
  for(int i=xstart;i<h1->GetSize();i++)
    sum += h1->GetBinContent(i);
  cntCompt.push_back(sum - ( cntFE.back() + cntSE.back() + cntDE.back() + cnt511.back()));


  // find number of Events simulatd in total
  searchKey = "/run/beamOn";
  checkPos = -1;
  checkValue="dummy//dummy";
  nEvents.push_back(GetValueFromConfigs(configValues, searchKey, checkPos, checkValue));
}

// print all of them
// cout << cntFE <<" " << cntSE << " " << cntDE << " " << cnt511 << " " << cntCompt << endl;
// cout << cntFE + cntSE + cntDE + cnt511 + cntCompt << endl; // consistency check
// cout << nEvents << " " << nCounts << endl;


// for(auto n : cntFE) {
//   cout << n << endl;
// }

// Write vectors to file
std::ofstream fout(fout_name);
fout << "# Eg[MeV?]" << "\t" 
<< "nEvents" << "\t" 
<< "nCounts" << "\t"  
<< "cntFE" << "\t"
<< "cntSE" << "\t"
<< "cntDE" << "\t"
<< "cnt511" << "\t" <<
"cntCompt=Rest" << endl;
for(vector<double>::size_type i = 0; i != cntFE.size(); i++) {
    fout 
    << Eg[i] << "\t"
    << nEvents[i] << "\t" 
    << nCounts[i] << "\t" 
    << cntFE[i] << "\t" 
    << cntSE[i] << "\t" 
    << cntDE[i] << "\t" 
    << cnt511[i] << "\t" 
    << cntCompt[i] << endl;
}
fout.close();

}