#include "AlpideClustering.h"
#include "CMDLineParser.h"
#include "AuxFunctions.h"
#include "libs.hh"
#include "spec.hh"

extern const std::string clusterise_help;

using namespace std;
using namespace AlpideClustering;
using namespace AlpideAuxFunctions;

vector<uint32_t> SetReadBranchAddresses(TTree* h101, uint64_t* ts, uint32_t* nPix, uint32_t (*Chip)[MAX_HITS], uint32_t (*Col)[MAX_HITS], uint (*Row)[MAX_HITS]) {
	if(!h101 || h101->IsZombie()) throw std::runtime_error("Bad TTree pointer passed to SetAllBranchAddress.");
	/* Try to find all MOSAIC##x##CHIP branches in the TTree */ 
	vector<uint32_t> valid_mosaics;
	TObjArray* branch_names = h101->GetListOfBranches();
	for(int x=0; x<256; ++x) {
		if(branch_names->FindObject(TString::Format("MOSAIC%dCHIP", x)) == nullptr) continue;
		
		h101->SetBranchAddress(TString::Format("MOSAIC%dT_LO", x), (uint32_t*)&ts[x]);
		h101->SetBranchAddress(TString::Format("MOSAIC%dT_HI", x), (uint32_t*)((char*)&ts[x] + 4));
		h101->SetBranchAddress(TString::Format("MOSAIC%dCHIP", x), &nPix[x]); // number of pixels fired == size of CHIPv,COLv,ROWv
		h101->SetBranchAddress(TString::Format("MOSAIC%dCHIPv", x), Chip[x]); // ChipId array
		h101->SetBranchAddress(TString::Format("MOSAIC%dCOLv", x), Col[x]);   // Col array
		h101->SetBranchAddress(TString::Format("MOSAIC%dROWv", x), Row[x]);   // Row array

		valid_mosaics.push_back(x);
	}
	if(valid_mosaics.size() == 0) throw std::runtime_error("SetReadBranchAddresses: couldn't find a single valid MOSAIC branch.\n\
			Expecting branches with names: MOSAIC%dCHIP etc.");
	return valid_mosaics;
}

void clusterise(const char* fileName, const char* outFile, ulong firstEvent=0, ulong maxEvents=0, int veto=0) {	
	TFile* in = new TFile(fileName,"READ");
    if(!in || in->IsZombie()) {cerr << "Can't open rootfile with name: " << fileName << "\n"; exit(EXIT_FAILURE);}
	
	auto t1 = timeNow();

    TTree* h101 = dynamic_cast<TTree*>(in->Get("h101"));
	if(veto < 0) veto=0;
	
	/* Read containers */
	uint64_t ts[256];
	uint64_t alpide_timestamp;
	uint32_t nPix[256];
	UInt_t Chip[256][MAX_HITS];
	UInt_t Col[256][MAX_HITS];
	UInt_t Row[256][MAX_HITS];
	auto valid_boards = SetReadBranchAddresses(h101, ts, nPix, Chip, Col, Row);
	
	/* MARK: Write containers */
	uint32_t cNum = 0;              // total count of clusters in the event
	uint32_t boardId[MAX_CLUSTERS]; // array of boardId for each cluster
	uint32_t chipId[MAX_CLUSTERS];  // array of chipId  for each cluster
	uint32_t cSize[MAX_CLUSTERS];   // array of cluster sizes
	double   uCol[MAX_CLUSTERS];    // mean x cluster
	double   uRow[MAX_CLUSTERS];    // mean x cluster
	double   sCol[MAX_CLUSTERS];    // sig x cluster
	double   sRow[MAX_CLUSTERS];    // sig y cluster

	/* For checking, we're branching also the raw data that got clustered */
	uint32_t _N = 0;                      // total number of pixels fired in the event                     
	uint32_t _boardId[MAX_CLUSTERS * 10]; // boardId for individual fired pixel
	uint32_t _chipId[MAX_CLUSTERS * 10];  // chipId of individual fired pixel
	uint32_t _ColV[MAX_CLUSTERS * 10];    // col of individual fired pixel
	uint32_t _RowV[MAX_CLUSTERS * 10];    // row of individual fired pixel

	TFile *out = new TFile(outFile, "RECREATE");
	TTree *tree = new TTree("h101", "h101");
	
	/* MARK: outward branches */
	// For timestamp, we can take the first board's. All the others are asserted to be within
	// stitch window in the drasi process
	tree->Branch("T", &alpide_timestamp, "T/l");
	tree->Branch("ALPIDE_cluster_count", &cNum);
	tree->Branch("ALPIDE_boardId", boardId, "ALPIDE_boardId[ALPIDE_cluster_count]/i");
	tree->Branch("ALPIDE_chipId", chipId, "ALPIDE_chipId[ALPIDE_cluster_count]/i");
	tree->Branch("ALPIDE_cluster_size", cSize, "ALPIDE_cluster_size[ALPIDE_cluster_count]/i");
	tree->Branch("ALPIDE_cluster_col", uCol, "ALPIDE_cluster_col[ALPIDE_cluster_count]/D");
	tree->Branch("ALPIDE_cluster_row", uRow, "ALPIDE_cluster_row[ALPIDE_cluster_count]/D");	
	tree->Branch("ALPIDE_cluster_col_sig", sCol, "ALPIDE_cluster_col_sig[ALPIDE_cluster_count]/D");
	tree->Branch("ALPIDE_cluster_row_sig", sRow, "ALPIDE_cluster_row_sig[ALPIDE_cluster_count]/D");	

	tree->Branch("_N", &_N);
	tree->Branch("_BOARD_ID", &_boardId, "_BOARD_ID[_N]/i");
	tree->Branch("_CHIP_ID", &_chipId, "_CHIP_ID[_N]/i");
	tree->Branch("_COLV", _ColV, "_COLV[_N]/i");
	tree->Branch("_ROWV", _RowV, "_ROWV[_N]/i");

	ulong lastEvent = SortEntries(firstEvent, maxEvents, h101);
	printf("Entries in file: %lld\n", h101->GetEntries());

    ulong evCounter(0);
    for(ulong evNum = firstEvent; evNum < lastEvent; ++evNum) {
        ++evCounter; if(evCounter%100 == 0) PrintProgress((double)evCounter/maxEvents);
		h101->GetEntry(evNum);
		cNum = 0; _N = 0;
		
		for(auto board : valid_boards) {
			// Check if board i has data in this event //
			if(nPix[board] == 0) continue;
			if(nPix[board] > MAX_HITS) {
				cerr << "Board: " << board << " has " << nPix[board] \
					 << " hits. Can't stuff them in the MAX_HITS = " \
					 << MAX_HITS << " array. \n" \
					 << " EvNum: " << evNum << endl << std::flush;
				continue;
			}
			alpide_timestamp = ts[board];
			
			auto chip_clusters_pairs = ConstructClusters(Chip[board], Col[board], Row[board], nPix[board], veto);
			// type: vector<pair<uint32_t, vector<vector<Point>>>>
			for(auto& [chip, clusters] : chip_clusters_pairs) {
				for(auto& cluster : clusters) {
					// Save cluster mean & sigma into the arguments //
					cSize[cNum] = FitCluster(cluster, uCol[cNum], uRow[cNum], sCol[cNum], sRow[cNum]); 
					chipId[cNum] = chip;

					boardId[cNum++] = board;

					for(auto& [px,py] : cluster) {
						_boardId[_N] = board;
						_chipId[_N] = chip;
						_ColV[_N] = px;
						_RowV[_N++] = py;
					}
				}
			}
		}

		if(cNum>0) tree->Fill();
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
#if 0
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
#endif
	/* ### ALPIDE_NUM is NOT defined at runtime! Once it's changed, the binary has to be recompiled. ### */

	clusterise(fileName.c_str(), outFile.c_str(), firstEvent, maxEvents, veto);
	
	cout<<endl; return 0;
}


