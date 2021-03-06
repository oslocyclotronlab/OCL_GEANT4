/*
 Usage:
    .x th1_to_mama.C+ ( histogram_name, "output-filename")
 or
    .x th1_to_mama.C+ ( histogram_name, "output-filename", "comment")
*/

#include <TH2.h>
#include <TH1.h>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <time.h>

void th1_to_mama(TH1* m, const char* filename, const char* comment="none")
{
    char tmp[128];
    time_t now = time(0);
    strftime(tmp, sizeof(tmp), "%Y-%m-%d %T", localtime(&now)); // not the 23-Mar-07 16:02:34 format

    TAxis *xax = m->GetXaxis();

    const int nx = xax->GetNbins();
    int n_max_Singles_Mama = 8191; // Different for TH1 and TH2!
    if(nx>n_max_Singles_Mama) {
            cout << "Spectrum has too many bins for mama; consider rebinning" << endl;
            // exit(1);
    }
    std::cout << "histogram is " << nx  << "; comment='" << comment << "'" << std::endl;

    std::ofstream mama(filename);
    mama << "!FILE=Disk \n"
    "!KIND=Matrix \n"
    "!LABORATORY=Oslo Cyclotron Laboratory (OCL) \n"
    "!EXPERIMENT=siri2root \n"
    "!COMMENT=" << comment << "\n"
    "!TIME=" << tmp << "\n"
    "!CALIBRATION EkeV=6";
    const double cal[3] = { xax->GetBinLowEdge(1), xax->GetBinWidth(1), 0};
    for(int i=0; i<3; ++i) {
        snprintf(tmp, sizeof(tmp), ",%13.6E", cal[i]);
        mama << tmp;
    }
    mama << "\n!PRECISION=16\n";
    snprintf(tmp, sizeof(tmp), "!DIMENSION=1,0:%4d\n", nx-1);
    mama << tmp;
    snprintf(tmp, sizeof(tmp), "!CHANNEL=(0:%4d)\n", nx-1);
    mama << tmp;

    for(int ix=1; ix<=nx; ++ix){
        mama << m->GetBinContent(ix) << ' ';
        mama << "\n";
    }
    mama << "!IDEND=" << endl << endl;

    mama.close();
}


