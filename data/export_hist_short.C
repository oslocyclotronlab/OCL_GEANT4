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

#include "th1_to_mama.C"
#include "th22mama.c"

#include <vector>
using namespace std;


/////////////////////////////////////////
/////////////////////////////////////////

void export_hist_short(){

TH1D* hmama; // output spectrum for mama
TH1D* hmama_raw; // before smoothing


std::vector<string> fnames;

string fout_name = "Peaks.dat";
system("xterm -e 'mkdir mama_spectra'"); // create dir for mama spectra, if not already existend
string outdir = "mama_spectra";

//////////////////////

// get list of files matching certain criterium

string fname_tmp;
system("xterm -e 'find root_files/*.root > tmp.txt'");
ifstream tmpfile("tmp.txt", ios::in);
while(tmpfile>>fname_tmp)
{
    fnames.push_back(fname_tmp); //adding data in to the vector
}

for(auto fname : fnames)
{
  cout << "Working on " << fname << endl;

  // read in root file
  auto fname_rootStr = fname; //+ ".root";
  const char *fname_root = fname_rootStr.c_str();
  TFile *file_root = new TFile(fname_root, "READ");

  TTree *OSCARhits;
  file_root->GetObject("OSCARhits",OSCARhits);

  // execute the macro that combines all detectors
  AnalyseSims t(OSCARhits);
  t.Loop();

  // get the histogram from AnalyseSims

  // TH1D *h1 = (TH1D*)gDirectory->GetList()->FindObject("h1");
  auto fname_mamaStr = outdir + "/" + fname + ".m";
  // const char *fname_mama = fname_mamaStr.c_str();
  // th1_to_mama(h1, fname_mama);

  TH2D *h1_all = (TH2D*)gDirectory->GetList()->FindObject("h1_all");
  fname_mamaStr = outdir + "/" + fname + "_all.m";
  const char *fname_mama_all = fname_mamaStr.c_str();
  th22mama(h1_all, fname_mama_all);
  }
}
