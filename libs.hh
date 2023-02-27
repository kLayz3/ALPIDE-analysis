#ifndef ALPIDE_LIBS_HH
#define ALPIDE_LIBS_HH

/* -- C++ -- */
#include <bits/stdc++.h>

/* -- ROOT -- */
#include "Riostream.h"
#include "TApplication.h"
#include "TArrow.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TCut.h"
#include "TCutG.h"
#include "TClonesArray.h"
#include "TDirectory.h"
#include "TEllipse.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGaxis.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH2.h"
#include "TH2F.h"
#include "TH3F.h"
#include "THStack.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLorentzVector.h"
#include "TMarker.h"
#include "TMath.h"
#include "TMatrixD.h"
#include "TMultiDimFit.h"
#include "TPaletteAxis.h"
#include "TPRegexp.h"
#include "TPrincipal.h"
#include "TProfile.h"
#include "TRandom1.h"
#include "TRandom2.h"
#include "TRandom3.h"
#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"
#include "TText.h"
#include "TVector3.h"
#include "TVectorD.h"
#include "TVectorT.h"
#include "Rtypes.h"

#define LEN(x) (sizeof x / sizeof *x)
#define timeNow() std::chrono::high_resolution_clock::now()

#define KMAG  "\x1B[35m"
#define KRED  "\x1B[31m"
#define KBLUE  "\x1B[34m"
#define KNRM  "\x1B[0m"
#define KGRN  "\x1B[32m"
#define KCYN  "\x1B[36m"

typedef uint32_t uint;
typedef uint64_t ulong;
using std::chrono::duration_cast;
using std::chrono::seconds;
using std::chrono::milliseconds;

inline extern const std::vector<std::string> maskFile {
    "",       // there's no ALPIDE0
    "16BC8",  // corresponds to MOSAIC0.cfg --ALPIDE1
    "12TB5",  // corresponds to MOSAIC1.cfg --ALPIDE2
    "17BA5",  // corresponds to MOSAIC3.cfg --ALPIDE3
    "16TA6",  // cooresponds to MOSAIC5.cfg --ALPIDE4
    "17BB2",  // cooresponds to MOSAIC7.cfg --ALPIDE5
    "17TB1",  // cooresponds to MOSAIC8.cfg --ALPIDE6
};


inline extern const std::string clusterise_help =
"\nUsage: ./clusterise <OPT1> <OPT2> ...\n\
\n\
--file=inputName.root        ..Input file.\n\
--first-event=N              ..Start from N-th event. Default 0. \n\
--max-events=N               ..Specify maximum number of events. Default all entries.\n\
--veto=N                     ..Only consider clusters with size>N. Default N=1.\n\
--output=/PATH/TO/OUT.root   ..Specify output file name. Default same as input file with 'cl' suffix.\n\
--dets=[d1,d2,..]            ..Condition to only write events which have clusters in specified detectors.\n\
      =all				     ..Equivalent to dets=1,2,...,ALPIDE_NUM. Every event must contain a cluster in all detectors.\n\
--help                       ..Print this message to stdout. \n\
\n\
The exe will cluster all the hits and write an output root file.\n\
\n\
Branch description:\n\
>> CL_NUM      : number of clusters in the event.\n\
>> ALPIDE_ID   : detector ID for each individual cluster.\n\
>> CL_SIZE     : cluster size for each cluster.\n\
>> CL_uCOL     : mean column position of each cluster.\n\
>> CL_uROW     : mean row position of each cluster.\n\
>> CL_uCOL_SIG : uncertainty of CL_uCOL, in units of col.\n\
>> CL_uROW_SIG : uncertainty of CL_uROW, in units of row.\n\
\n\
Additional branches to show individual pixels clustered:\n\
>> _N    : total number of pixels fired. Equals to sum of all CL_SIZE.\n\
>> _COLV : pixels fired which belong to a cluster. First CL_SIZE[0] elements belong to 0th cluster,\n\
>>         next CL_SIZE[1] belong to 1st cluster, etc \n\
>> _ROWV : same but for rows.\n\
\n\
Good luck, have fun <(^.^)>\n\n";

inline extern const std::string calibrate_help =
"\nUsage: ./calibrate <OPT1> <OPT2> ...\n\
\n\
file=inputName.root          ..Input file.\n\
--first-event=N              ..Start from N-th event. Default 0. \n\
--max-events=N               ..Specify maximum number of events. Default all entries.\n\
--output=/PATH/TO/OUT.root   ..Specify output file name. Default same as input file with 'calib' suffix.\n\
--help                       ..Print this message to stdout. \n\
\n\
The exe will calibrate the detectors and write an output root file.\n\
Calibration of (dx,dy) for each ALPIDE detector with reference to the first ALPIDE is stored in the\n\
aCol, bCol, aRow, bRow and corresponding sigma branches. 'a' is the slope of the fit line and 'b' the offset.\n\
\n\
Good luck, have fun <(^.^)>\n\n";

inline extern const std::string analyse_help =
"\nUsage: ./analyse <OPT1> <OPT2> ...\n\
\n\
--file=inputName.root        --Input file.\n\
--first-event=N              ..Start from N-th event. Default 0 \n\
--max-events=N	             ..Specify maximum number of events. Default all entries.\n\
\n\
--hitmap=X                   ..Plots the hitmap of AlpideX.\n\
--hitmap		             ..Plots the hitmaps of all Alpides.\n\
--track			             ..Do the tracking.\n\
	--cal=calFile.root       ..Pass calibration file to the tracker (OPTIONAL).\n\
	--save=saveFile.root     ..Save the tracks into a rootfile (OPTIONAL).\n\
--help                       ..Print this message to stdout. \n\
\n\
This exe can plot the hitmaps of all or specific detector. Works with raw or clustered data.\n\
To do tracking pass --track --cal=<calFile>.root flags. Works only for input files with clustered data. \n\
\n\
Good luck, have fun <(^.^)>\n\n";

#endif
