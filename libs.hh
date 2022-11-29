/* -- C++ -- */
#include <vector>
#include <string>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <cassert>
#include <exception>
#include <regex>
#include <algorithm>
#include <chrono>

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

extern const std::vector<std::string> maskFile{
    "",       // there's no ALPIDE0
    "16BC8",  // corresponds to MOSAIC0.cfg --ALPIDE1
    "12TB5",  // corresponds to MOSAIC1.cfg --ALPIDE2
    "17BA5",  // corresponds to MOSAIC3.cfg --ALPIDE3
    "16TA6",  // cooresponds to MOSAIC5.cfg --ALPIDE4
    "17BB2",  // cooresponds to MOSAIC7.cfg --ALPIDE5
    "17TB1",  // cooresponds to MOSAIC8.cfg --ALPIDE6
};

extern const std::string analyse_help =
"\nUsage: ./analyse <OPT1> <OPT2> ...\n\
file=/PATH/TO/FILE/file.root\n\
--first-event=N   Start from N-th event. Default 0 \n\
--max-events=N	  Specify maximum number of events. Default all entries.\n\
--veto=N          Remove clusters with size<=N.\n\
\n\
--raw=X,Y         Plots only raw correlations of AlpideX:AlpideY\n\
--raw=X           Plots the raw hitmap of AlpideX  \n\
\n\
The script can cluster all the hits and make analysis on the clustered hits.\n\n\
Good luck, have fun <(^.^)>\n\n";

extern const std::string clusterise_help =
"\nUsage: ./clusterise <OPT1> <OPT2> ...\n\
file=/PATH/TO/FILE/file.root\n\
--first-event=N   Start from N-th event. Default 0. \n\
--max-events=N    Specify maximum number of events. Default all entries.\n\
--veto=N          Remove clusters with size<=N.\n\
--output=/PATH/TO/FILE/file.root \n\
				  Specify output file name. Default same as input file with cl suffix.\n\
\n\
The exe will cluster all the hits in a (selected) format and write an output root file.\n\
A cluster is represented as a tuple<float,float,float,float,uint> \n\
corresponding to (meanX, meanY, sigmaX, sigmaY, N) of the cluster. N is the size of the cluster.\n\
Good luck, have fun <(^.^)>\n\n";

/* --raw=0,1         Writes output root tree as a list of branches of vectors of AlpideClustering::Point structures\n\ */
/* 				  (equivalent to std::pair<uint,uint>) corresponding to all fired pixels in the same cluster. \n\ */
/* 				  If unspecified, the cluster is represented as a tuple<float,float,float,float,uint> \n\ */
/* 				  corresponding to (meanX, meanY, sigmaX, sigmaY, N) of the cluster. N is the size of the cluster.\n\ */
/* 				  File structure can be checked in detail using the parser (to-be-done).\n\ */
