#include "includes/AlpideClustering.cpp"
#include "includes/CMDLineParser.cc"
#include "includes/AuxFunctions.cc"
#include "libs.hh"
#include "spec.hh"

#define LEN(x) (sizeof x / sizeof *x)
#define timeNow() std::chrono::high_resolution_clock::now()
typedef uint32_t uint;
typedef uint64_t ulong;
extern const std::string analyse_help; 

using namespace std;
using namespace AlpideClustering;
using std::chrono::duration_cast;
using std::chrono::seconds;

void SetAllBranchAddress(TTree* h101, uint& cNum, uint* AlpideID, uint* cSize, float* uCol, float* uRow) {
	if(!h101 || h101->IsZombie()) return;
	h101->SetBranchAddress("CL_NUM", &cNum);
	h101->SetBranchAddress("ALPIDE_ID", AlpideID);
	h101->SetBranchAddress("CL_uCOL", uCol);
	h101->SetBranchAddress("CL_uROW", uRow);
} 

ulong SortEntries(ulong& firstEvent, ulong& maxEvents, TTree* h101) {
	firstEvent = std::min(firstEvent, (ulong)h101->GetEntries());
    ulong n = (maxEvents==0 || firstEvent+maxEvents > h101->GetEntries()) ? (h101->GetEntries()) : (firstEvent+maxEvents);
    maxEvents = n - firstEvent;
	return n;
}

void RawHitMap(const char* fileName, int x, ulong firstEvent=0, ulong maxEvents=0) {
    assert(x>=1 && x<=ALPIDE_NUM);
    TApplication* app = new TApplication("myApp", 0, 0);
    TFile* in = new TFile(fileName,"READ");
    if(!in || in->IsZombie()) {cerr << "Can't open rootfile with name: " << fileName << "\n"; exit(EXIT_FAILURE);}
    TTree* h101 = static_cast<TTree*>(in->Get("h101"));
    auto t1 = timeNow();

	/* Read containers */
	uint cNum;
	uint AlpideID[MAX_CLUSTERS];
	uint cSize[MAX_CLUSTERS];
	float uCol[MAX_CLUSTERS];
	float uRow[MAX_CLUSTERS];

    SetAllBranchAddress(h101, cNum, AlpideID, cSize, uCol, uRow);

    TH2D* hRawHit = new TH2D(TString::Format("Raw Hitmap ALPIDE%d", x), TString::Format("Raw Hitmap ALPIDE%d", x), 1024,0,1024,512,0,512);  

	ulong lastEvent = SortEntries(firstEvent, maxEvents, h101);
	printf("Entries in file: %lld\n", h101->GetEntries());
    ulong evCounter{0};
    
    for(ulong evNum = firstEvent; evNum < lastEvent; ++evNum) {
        ++evCounter; if(evCounter%100 == 0) PrintProgress((float)evCounter/maxEvents);		
		h101->GetEntry(evNum);

		for(uint c=0; c<cNum; ++c) {
			if(AlpideID[c] != x) continue;
			hRawHit->Fill(uCol[c], uRow[c], (double)cSize[c]);
		}
    }
    
	hRawHit->Draw("colz");

	auto t2 = timeNow();
    cout << "\nTime taken: " << duration_cast<seconds>(t2-t1).count() << "s\n";
    app->Run();

	in->Close();
}


auto main(int argc, char* argv[]) -> int {
    if(IsHelpArg(argc, argv)) {cout << analyse_help; return 0;}
    
    string pStr;
    string fileName("");
    ulong firstEvent(0);
    ulong maxEvents(0);
    
	if(!ParseCmdLine("file", fileName, argc, argv)) {
		cerr << "No file specified!\n";
		cout << clusterise_help; return 0;
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
	if(ParseCmdLine("raw", pStr, argc, argv)) {
		int x;
		try {
			x = stoi(pStr);
			RawHitMap(fileName.c_str(), x, firstEvent, maxEvents);
		}
		catch(exception& e) {}
	}

    cout<<endl; return 0;
}


