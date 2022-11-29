#include "includes/AlpideClustering.cpp"
#include "includes/CMDLineParser.cc"
#include "libs.hh"
#include "spec.hh"

#define LEN(x) (sizeof x / sizeof *x)
#define timeNow() std::chrono::high_resolution_clock::now()

typedef uint32_t uint;
typedef uint64_t ulong;
extern const std::string clusterise_help;

using namespace std;
using std::chrono::duration_cast;
using std::chrono::seconds;
using namespace AlpideClustering;

/* Read containers */
UInt_t Row[ALPIDE_NUM+1][MAX_HITS];
UInt_t Col[ALPIDE_NUM+1][MAX_HITS];
UInt_t rowM[ALPIDE_NUM+1];
UInt_t colM[ALPIDE_NUM+1];

/* Write containers */
vector<float> uCol[ALPIDE_NUM+1];  // mean X of the clusters
vector<float> uRow[ALPIDE_NUM+1];  // mean Y of the clusters
vector<float> sCol[ALPIDE_NUM+1];  // sig  X of the clusters
vector<float> sRow[ALPIDE_NUM+1];  // sig  Y of the clusters
vector<uint>  cSize[ALPIDE_NUM+1]; // sizes of each cluster
/* All containers are sensor specific! */

void SetOneBranchAddress(TTree* h101, int x) {
    assert(x<=ALPIDE_NUM && x>=1);
    h101->SetBranchAddress(TString::Format("ALPIDE%dCOLv", x), Col[x]);
    h101->SetBranchAddress(TString::Format("ALPIDE%dROWv", x), Row[x]);
    h101->SetBranchAddress(TString::Format("ALPIDE%dCOL", x), &colM[x]);
    h101->SetBranchAddress(TString::Format("ALPIDE%dROW", x), &rowM[x]);
}

void SetAllBranchAddress(TTree* h101) {
	if(!h101) return;
	for(int x=1; x<=ALPIDE_NUM; ++x) {
		h101->SetBranchAddress(TString::Format("ALPIDE%dCOLv", x), Col[x]);
		h101->SetBranchAddress(TString::Format("ALPIDE%dROWv", x), Row[x]);
		h101->SetBranchAddress(TString::Format("ALPIDE%dCOL", x), &colM[x]);
		h101->SetBranchAddress(TString::Format("ALPIDE%dROW", x), &rowM[x]);
	}
}

ulong sortEntries(ulong& firstEvent, ulong& maxEvents, TTree* h101) {
	firstEvent = std::min(firstEvent, (ulong)h101->GetEntries());
    ulong n = (maxEvents==0 || firstEvent+maxEvents > h101->GetEntries()) ? (h101->GetEntries()) : (firstEvent+maxEvents);
    maxEvents = n - firstEvent;
	return n;
}

/* TODO Later: This will clusterise hits in all alpides (specified by the directive ALPIDE_NUM) 
void RawClusterise(const char* fileName, const char* outFile = nullptr, ulong firstEvent=0, ulong maxEvents=0, int veto=1) {}
*/

/* Cluster being represented as a tuple of (meanX,meanY,sigmaX,sigmaY,N) */
void CoarseClusterise(const char* fileName, const char* outFile, ulong firstEvent=0, ulong maxEvents=0, int veto=1) {
	TFile* in = new TFile(fileName,"READ");
    if(!in || in->IsZombie()) {cerr << "Can't open rootfile with name: " << fileName << "\n"; exit(EXIT_FAILURE);}
    TTree* h101 = (TTree*)in->Get("h101");
	if(veto < 0) veto=0;	
	
	TFile *out = new TFile(outFile, "RECREATE");
	TTree *tree = new TTree("h101", "h101");
	
	for(int i=1; i<=ALPIDE_NUM; ++i) {
		tree->Branch(TString::Format("A%dcSize",i), &cSize[i]);
		tree->Branch(TString::Format("A%duCOL", i), &uCol[i]);
		tree->Branch(TString::Format("A%duROW", i), &uRow[i]);	
		tree->Branch(TString::Format("A%dsCOL", i), &sCol[i]);	
		tree->Branch(TString::Format("A%dsROW", i), &sRow[i]);
	}
	
	auto t1 = timeNow();
	SetAllBranchAddress(h101);

	ulong lastEvent = sortEntries(firstEvent, maxEvents, h101);
	printf("Entries in file: %lld\n", h101->GetEntries());

    ulong evCounter(0);
    for(ulong evNum = firstEvent; evNum < lastEvent; ++evNum) {
        ++evCounter;
        if(evCounter%100 == 0) {
			printProgress((float)evCounter/maxEvents);
		}
		h101->GetEntry(evNum);
		bool hasCluster(false); /* doesn't fill the file if no cluster identified */
		
		/* i loops over all alpides {1,2,3 ... MAX_ALPIDES} */
		for(int i=1; i<=ALPIDE_NUM; ++i) {
			if(rowM[i]==0 || rowM[i] != colM[i]) continue;
			auto clusters = ConstructClusters(Col[i], Row[i], rowM[i], veto);
			if(clusters.size() == 0) continue;

			for(auto& cluster: clusters) {
				float uX,uY,sX,sY;
				uint cN = FitCluster(cluster, uX,uY,sX,sY);
				uRow[i].push_back(uX); uCol[i].push_back(uY);
				sRow[i].push_back(sX); sCol[i].push_back(sY);	
				cSize[i].push_back(cN);
			}
			hasCluster = true;
		}
		if(hasCluster) tree->Fill();
		
		/* clear the vectors for the next event */
		for_each(uCol, uCol+ALPIDE_NUM, [](auto& vec){vec.clear();});
		for_each(sCol, sCol+ALPIDE_NUM, [](auto& vec){vec.clear();});
		for_each(uRow, uRow+ALPIDE_NUM, [](auto& vec){vec.clear();});
		for_each(sRow, sRow+ALPIDE_NUM, [](auto& vec){vec.clear();});
		for_each(cSize, cSize+ALPIDE_NUM, [](auto& vec){vec.clear();});
	}

	tree->Write();
    auto t2 = timeNow();
    cout << "\nTime taken: " << duration_cast<seconds>(t2-t1).count() << "s\n";
}

auto main(int argc, char* argv[]) -> int {
    if(IsHelpArg(argc, argv)) {cout << clusterise_help; return 0;}
	
	string pStr;
    string fileName;
	string outFile;
    ulong firstEvent(0);
    ulong maxEvents(0);
    int veto(1);
    
    if(!ParseCmdLine("file", fileName, argc, argv)) {
		cerr << "No file specified!\n";
		cout << clusterise_help; return 0;
	}
	if(!ParseCmdLine("output", outFile, argc, argv)) {
		outFile = fileName.substr(0, fileName.find('.')) + "_coarse_cl.root";
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
/* void CoarseClusterise(const char* fileName, const char* outFile, ulong firstEvent=0, ulong maxEvents=0, int veto=1) { */
	CoarseClusterise(fileName.c_str(), outFile.c_str(), firstEvent, maxEvents, veto);
	cout<<endl; return 0;
}


