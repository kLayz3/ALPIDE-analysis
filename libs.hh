/* -- C++ -- */
#include <vector>
#include <string>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <cassert>
#include <exception>
#include <regex>
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

extern const std::string helpString =
"\nUsage: ./analysis <OPT>\n\
file=/PATH/TO/FILE/file.root\n\
--first-event=N   Start from N-th event. Default 0 \n\
--max-events=N    Fix maximum number of events. Default Nentries.\n\
--veto=N          Remove clusters with size<=N.\n\
\n\
--raw=X,Y         Plots only raw correlations of AlpideX:AlpideY\n\
--raw=X           Plots the raw hitmap of AlpideX  \n\
\n\
The script will cluster all the hits and make analysis on the clustered hits.\n\n\
Good luck, have fun <(^.^)>\n\n";
