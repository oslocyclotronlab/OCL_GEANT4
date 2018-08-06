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
#include <iomanip>
#include <cstdio>
#include "TFile.h"

#include "AnalyseSims.C"
#ifdef MAKECINT
#pragma link C++ defined_in “AnalyseSims.h”;
#pragma link C++ defined_in “AnalyseSims.C”;
#endif

#include "th22mama_hist.C"

#include <vector>
using namespace std;

////////////////////////////
// Parameters
// double threshold = 300.; // Detector threshold in MeV
double threshold = 0.; // Detector threshold in MeV
// double smoothingFactor = 100.; // smoothing for the mama output 
// smoothing now externally: works badly with root, creates spikes!
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
  int nchTot = 2*nChannelsBG +1;
  int nchBGonly = nchTot -3;
  double cntbg = (sum - cntsInPeak) / nchBGonly; // bg per channel
  double cntWobg = cntsInPeak - cntbg*3.;
  return cntWobg;

}


TH1D* CreateMamaSpectrum(TH1D* h, double EgFE, int nChannelsBG=3)
{
  TH1D* h1 = (TH1D*) h->Clone(); // create output spectrum as a clone
  h1->SetTitle("hmama");
  h1->SetName("hmama");

  vector<double> peaks;
  if(EgFE>=1022.)
  {
  peaks.insert(peaks.end(), {EgFE, EgFE-511., EgFE-2.*511., 511.});
  }
  else
  {
  peaks.insert(peaks.end(), {EgFE});
  }

  TAxis *xaxis = h1->GetXaxis();
  TAxis *yaxis = h1->GetYaxis();


  for(auto Eg : peaks)
  {
    Int_t  binx     = xaxis->FindBin(Eg);
    double binvalue = h1->GetBinContent(binx);

    double cntsInPeak = binvalue;
    cntsInPeak += h1->GetBinContent(binx+1) + h1->GetBinContent(binx-1);

    // get BG
    Int_t min = binx - nChannelsBG;
    Int_t max = binx + nChannelsBG;
    // cout << "Subtract BG from " <<  h1->GetXaxis()->GetBinCenter(min)
    // << "to " <<  h1->GetXaxis()->GetBinCenter(max) << endl;
    double sum = 0;
    for(int ch=min;ch<max+1;ch++ ){
      sum += h1->GetBinContent(ch);
    }
    // a little messy -- this accounts also for thaat counts can come in the bin below/above
    int nchTot = 2*nChannelsBG +1;
    int nchBGonly = nchTot -3;
    double cntbg = (sum - cntsInPeak) / nchBGonly; // bg per channel

    // set all counts equal to bg
    for(int ch=min;ch<max+1;ch++ ){
      h1->SetBinContent(ch,cntbg);
    }
  }

  return h1;
}

// // generate a random number
// TRandom *r3 = new TRandom3(0.); // initiated with random seed

// TH1D* SmoothSpectrum(TH1D* h, double smoothingFactor){

//   TH1D* h1 = (TH1D*) h->Clone(); // create output spectrum as a clone
//   int xstart = 0;
//   double E;
//   // for(int i=xstart;i<h1->GetSize();i++){
//   //   h1->SetBinContent(i,0);
//   // }
//   for(int i=xstart;i<h1->GetSize();i++){
//     E = h1->GetBinCenter(i); 
//     h1->Fill(r3->Gaus(E , sqrt(smoothingFactor*E )));
//     if(i<50)
//       {cout << r3->Gaus(E , sqrt(smoothingFactor*E ))<< endl;}
//   }
//   return h1;
// }


template <typename T>
std::string to_string_with_precision(const T a_value, const int n = 2)
{
    std::ostringstream out;
    out << std::setprecision(n) << a_value;
    return out.str();
}



TH1D* ScaleX(TH1D* h, double scaleFactor){
  // Scale histogram X axis 
  // double ScaleX = 1.e3; //

  int nBinsx = h->GetXaxis()->GetNbins();
  double xmin = h->GetXaxis()->GetXmin();
  double xmax = h->GetXaxis()->GetXmax();
  cout << "h before scaling"<< endl;;
  cout << "xmin: " << xmin << "; xmax:" << xmax << endl;
  xmax *= scaleFactor;
  h->GetXaxis()->Set(nBinsx,xmin,xmax);
  cout << "h after scaling" << endl;;
  cout << "xmin: " << xmin << "; xmax:" << xmax << endl;
    return h;
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
vector<double> cntRest;
vector<double> cntRestNoThres;
Int_t nChannelsBG;// # of channels, symmetric, to choose to average out the BG around the peaks
int xstart; // dummy for startbin
int xstop; // dummy for stop bin
double sum; // dummy for a sum

TH1D* hmama; // output spectrum for mama
TH1D* hmama_raw; // before smoothing


std::vector<string> fnames;

string fout_name = "Peaks.dat";
system("xterm -e 'mkdir mama_spectra'"); // create dir for mama spectra, if not already existend
string outdir = "mama_spectra";

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
// fnames.insert(fnames.end(), { "sim29.mac"});

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
  // EgFE /=1000.; //  keV to MeV
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
  h1=ScaleX(h1, 1e3);
  h1->Draw();
  //

  // cntFE = GetBinFromXval(h1, EgFE);
  cntFE.push_back(GetBinWoBG(h1,EgFE,nChannelsBG));

  // get peaks on top of bg -- 
  nChannelsBG = 3;
  cntSE.push_back(GetBinWoBG(h1,EgFE-511.,nChannelsBG));
  cntDE.push_back(GetBinWoBG(h1,EgFE-2.*511.,nChannelsBG));
  cnt511.push_back(GetBinWoBG(h1,511.,nChannelsBG));
  if(EgFE<1022.){
    cntSE.back()=0;
    cntDE.back()=0;
    cnt511.back()=0;
  }
  nCounts.push_back(h1->GetEntries());

  // assume that compton is all but the peaks; staring at threshold
  xstart = h1->GetXaxis()->FindBin(threshold);
  sum = 0;
  for(int i=xstart;i<h1->GetSize();i++){
    sum += h1->GetBinContent(i);
  }
  cntRest.push_back(sum - ( cntFE.back() + cntSE.back() + cntDE.back() + cnt511.back()));

  // find number of Events simulatd in total
  searchKey = "/run/beamOn";
  checkPos = -1;
  checkValue="dummy//dummy";
  nEvents.push_back(GetValueFromConfigs(configValues, searchKey, checkPos, checkValue));

  // Create spectra of "compton" BG, here= everything that was not peaks
  hmama = CreateMamaSpectrum(h1, EgFE, nChannelsBG);
  // hmama->Smooth(smoothingFactor);
  hmama->Draw();

  cntRestNoThres.push_back(hmama->GetEffectiveEntries());

  // string Egpp = to_string_with_precision(EgFE/1000, 0);
  string Egpp = to_string(int(EgFE));

  auto fname_mamaStr = outdir + "/" + "cmp" + Egpp;
  const char *fname_mama = fname_mamaStr.c_str();
  th22mama_hist(hmama, fname_mama);

  // cout << nCounts.back() - cntFE.back() -  cntSE.back() - cntDE.back() - cnt511.back() << "\t" << cntRestNoThres.back() << endl;

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
<< "cnt511" << "\t"
<< "cntRest" << "\t" <<
"cntRestNoThres" << endl;
for(vector<double>::size_type i = 0; i != cntFE.size(); i++) {
    fout 
    << Eg[i] << "\t"
    << nEvents[i] << "\t" 
    << nCounts[i] << "\t" 
    << cntFE[i] << "\t" 
    << cntSE[i] << "\t" 
    << cntDE[i] << "\t" 
    << cnt511[i] << "\t" 
    << cntRest[i] << "\t" 
    << cntRestNoThres[i] << endl;
}
fout.close();

}