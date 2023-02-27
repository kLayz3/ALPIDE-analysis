#include "AlpideClustering.h"
#include "CMDLineParser.h"
#include "AuxFunctions.h"
#include "libs.hh"
#include "spec.hh"

extern const std::string clusterise_help;

using namespace std;
using namespace AlpideClustering;

void SetOneBranchAddress(TTree* h101, int x, uint32_t* Col, uint* Row, uint& rowM, uint& colM) {
    assert(x<=ALPIDE_NUM && x>=1);
    h101->SetBranchAddress(TString::Format("ALPIDE%dCOLv", x), Col);
    h101->SetBranchAddress(TString::Format("ALPIDE%dROWv", x), Row);
    h101->SetBranchAddress(TString::Format("ALPIDE%dCOL", x), &colM);
    h101->SetBranchAddress(TString::Format("ALPIDE%dROW", x), &rowM);
}

void SetAllBranchAddress(TTree* h101, uint32_t (*Col)[MAX_HITS], uint (*Row)[MAX_HITS], uint* rowM, uint* colM, uint& tHi, uint& tLo) {
	if(!h101 || h101->IsZombie()) return;
	for(int x=1; x<=ALPIDE_NUM; ++x) {
		h101->SetBranchAddress(TString::Format("ALPIDE%dCOLv", x), Col[x]);
		h101->SetBranchAddress(TString::Format("ALPIDE%dROWv", x), Row[x]);
		h101->SetBranchAddress(TString::Format("ALPIDE%dCOL", x), &colM[x]);
		h101->SetBranchAddress(TString::Format("ALPIDE%dROW", x), &rowM[x]);
	}
	h101->SetBranchAddress("ALPIDE1T_HI", &tHi);
	h101->SetBranchAddress("ALPIDE1T_LO", &tLo);
}

void CoarseClusterise(const char* fileName, const char* outFile, ulong firstEvent=0, ulong maxEvents=0, int veto=0, const std::bitset<ALPIDE_NUM+1> mandatory = {0}) {	
	TFile* in = new TFile(fileName,"READ");
    if(!in || in->IsZombie()) {cerr << "Can't open rootfile with name: " << fileName << "\n"; exit(EXIT_FAILURE);}
    TTree* h101 = static_cast<TTree*>(in->Get("h101"));
	if(veto < 0) veto=0;
	
	/* Read containers */
	UInt_t tHi, tLo;
	UInt_t Col[ALPIDE_NUM+1][MAX_HITS];
	UInt_t Row[ALPIDE_NUM+1][MAX_HITS];
	UInt_t colM[ALPIDE_NUM+1];
	UInt_t rowM[ALPIDE_NUM+1];
	SetAllBranchAddress(h101, Col, Row, colM, rowM, tHi, tLo);
	
	/* MARK: Write containers */
	uint32_t cNum{0}; // number of clusters
	uint32_t AlpideID[MAX_CLUSTERS]; 
	uint32_t cSize[MAX_CLUSTERS];

	double uCol[MAX_CLUSTERS];
	double uRow[MAX_CLUSTERS];
	double uCol_SIG[MAX_CLUSTERS];
	double uRow_SIG[MAX_CLUSTERS];

	uint32_t _N{0};
	uint32_t _ColV[MAX_CLUSTERS * 10];
	uint32_t _RowV[MAX_CLUSTERS * 10];

	TFile *out = new TFile(outFile, "RECREATE");
	TTree *tree = new TTree("h101", "h101");
	
	/* MARK: outward branches */
	tree->Branch("T_HI", &tHi);
	tree->Branch("T_LO", &tLo);
	tree->Branch("CL_NUM", &cNum);
	tree->Branch("ALPIDE_ID", AlpideID, "ALPIDE_ID[CL_NUM]/i");
	tree->Branch("CL_SIZE", cSize, "CL_SIZE[CL_NUM]/i");
	tree->Branch("CL_uCOL", uCol, "CL_uCOL[CL_NUM]/D");
	tree->Branch("CL_uROW", uRow, "CL_uROW[CL_NUM]/D");	
	tree->Branch("CL_uCOL_SIG", uCol_SIG, "CL_uCOL_SIG[CL_NUM]/D");
	tree->Branch("CL_uROW_SIG", uRow_SIG, "CL_uROW_SIG[CL_NUM]/D");	
	tree->Branch("_N", &_N);
	tree->Branch("_COLV", _ColV, "_COLV[_N]/i");
	tree->Branch("_ROWV", _RowV, "_ROWV[_N]/i");

	auto t1 = timeNow();

	ulong lastEvent = SortEntries(firstEvent, maxEvents, h101);
	printf("Entries in file: %lld\n", h101->GetEntries());

    ulong evCounter(0);
    for(ulong evNum = firstEvent; evNum < lastEvent; ++evNum) {
        ++evCounter; if(evCounter%100 == 0) PrintProgress((double)evCounter/maxEvents);
		h101->GetEntry(evNum);
	
		cNum = 0; _N = 0;
		std::bitset<ALPIDE_NUM+1> b = mandatory;
		
		/* MARK: clustering, i loops over all alpides {1,2,3 ... ALPIDE_NUM} */
		for(int i=1; i<=ALPIDE_NUM; ++i) {
			if(rowM[i]==0 || rowM[i] != colM[i]) continue;

			auto clusters = ConstructClusters(Col[i], Row[i], rowM[i], veto);
			if(clusters.size() == 0) continue;

			b[i]=0;

			for(auto& cluster : clusters) {
				// Save cluster mean & sigma into the arguments //
				cSize[cNum] = FitCluster(cluster, uCol[cNum],uRow[cNum], uCol_SIG[cNum],uRow_SIG[cNum]); 
				AlpideID[cNum] = i;
				++cNum;

				for(auto [px,py] : cluster) {
					_ColV[_N] = px;
					_RowV[_N] = py;
					++_N;
				}
			}
		}

		if(cNum>0 && b.none())
			tree->Fill();
	}

	out->Write();
	out->Close();
	in->Close();
	//ReleaseMalloc(AlpideID, cSize, uCol, uRow);

    auto t2 = timeNow();
    cout << "\nTime taken: " << duration_cast<seconds>(t2-t1).count() << "s\n";
}

auto main(int argc, char* argv[]) -> int {
    if(IsCmdArg("help", argc, argv)) {cout << clusterise_help; return 0;}
	
	string pStr;
    string fileName;
	string outFile;
    ulong firstEvent(0);
    ulong maxEvents(0);
    int veto(0);
	std::bitset<ALPIDE_NUM+1> mandatory{0};
    
    if(!ParseCmdLine("file", fileName, argc, argv)) {
		cerr << "No file specified!\n";
		cout << clusterise_help; return 0;
	}
	if(!ParseCmdLine("output", outFile, argc, argv)) {
		outFile = fileName.substr(0, fileName.find('.')) + "_cl.root";
		cout << "No output file specified. Writing into file: " << outFile << endl;
	}
    
	if(ParseCmdLine("veto", pStr, argc, argv)) {
        try {
            veto = stoi(pStr);
            printf("Only taking clusters size > %d\n", veto);
        }
        catch(exception& e) {}
    }
    
    if(ParseCmdLine("first-event", pStr, argc, argv)) {
        try {
            firstEvent = stoul(pStr);
            printf("Starting from ev#: %ld\n", firstEvent);
        }
        catch(exception& e) {}
    }
    
    if(ParseCmdLine("max-events", pStr, argc, argv)) {
        try {
            maxEvents = stoul(pStr);
            printf("Max events: %ld\n", maxEvents);
        }
        catch(exception& e) {}
    }
	if(ParseCmdLine("dets", pStr, argc, argv)) {
		RemoveCharsFromString(pStr, "[({})]");
		if(pStr=="all" || pStr=="a") {
			mandatory.set(); 
			mandatory[0] = 0; // sets all to 1
		}
		else {
			auto parts = SplitStringToSet(pStr, ',');
			for(string part : parts) {
				try {
					int x = stoi(part);
					if(x>ALPIDE_NUM || x<=0) throw std::runtime_error("Out of Bounds.");
					mandatory[x] = 1;
				}
				catch(exception& e) {cout << "In " << __PRETTY_FUNCTION__ << ": " << e.what() << endl;}
			}
		}
	}

	/* ### ALPIDE_NUM is NOT defined at runtime! Once it's changed, the binary has to be recompiled. ### */

	CoarseClusterise(fileName.c_str(), outFile.c_str(), firstEvent, maxEvents, veto, mandatory);
	
	cout<<endl; return 0;
}


