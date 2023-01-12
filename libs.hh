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
--veto=N                     ..Only consider clusters with size>N. Default N=2.\n\
--output=/PATH/TO/OUT.root   ..Specify output file name. Default same as input file with 'cl' suffix.\n\
--dets=[d1,d2,..]            ..Condition to only write events which have clusters in specified detectors.\n\
      =all				     ..Equivalent to dets=1,2,...,ALPIDE_NUM. Every event must contain a cluster in all detectors.\n\
--help                       ..Print this message to stdout. \n\
\n\
The exe will cluster all the hits and write an output root file.\n\
A cluster is represented as a tuple<float,float,uint> \n\
corresponding to (meanX, meanY, N) of the cluster. N is the size of the cluster.\n\
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
Calibration of (dx, dy) for each alpide with reference to the first ALPIDE is stored in the\n\
aCol, bCol, aRow, bRow and corresponding sigma branches. 'a' is the slope of the fit line and 'b' the offset.\n\
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
	--cal=calFile.root       ..Pass calibration file to the tracker. Must have.\n\
	--save=saveFile.root     ..Save the tracks into a rootfile (OPTIONAL).\n\
--help                       ..Print this message to stdout. \n\
\n\
This exe can plot the hitmaps of all or specific detector. Works with raw or clustered data.\n\
To do tracking pass --track --cal=<calFile>.root flags. Works only for input files with clustered data. \n\n\
Good luck, have fun <(^.^)>\n\n";

#endif
