#include "includes/AlpideClustering.cpp"
#include "libs.hh"
#include "spec.hh"

#define LEN(x) (sizeof x / sizeof *x)

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
vector<vector<Point>> clusters[ALPIDE_NUM+1];
UInt_t clustM[ALPIDE_NUM+1];
vector<std::tuple<float,float,float,float,uint>> fitVector[ALPIDE_NUM+1];

bool IsHelpArg(int, char**);
bool ParseCmdLine(const char*, std::string&, int, char**);
void printProgress(double);

void SetOneBranchAddress(TTree* h101, int x) {
    assert(x<=ALPIDE_NUM && x>=1); 
    h101->SetBranchAddress(Form("ALPIDE%dCOLv", x), Col[x]);
    h101->SetBranchAddress(Form("ALPIDE%dROWv", x), Row[x]);
    h101->SetBranchAddress(Form("ALPIDE%dCOL", x), &colM[x]);
    h101->SetBranchAddress(Form("ALPIDE%dROW", x), &rowM[x]);
}

void SetAllBranchAddress(TTree* h101) {
	if(!h101) return;
	for(int x=1; x<=ALPIDE_NUM; ++x) {
		h101->SetBranchAddress(Form("ALPIDE%dCOLv", x), Col[x]);
		h101->SetBranchAddress(Form("ALPIDE%dROWv", x), Row[x]);
		h101->SetBranchAddress(Form("ALPIDE%dCOL", x), &colM[x]);
		h101->SetBranchAddress(Form("ALPIDE%dROW", x), &rowM[x]);
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
	
	TString outputFileName;
	if(!outFile || LEN(outFile)==0) outputFileName += strtok(fileName, ".").str() + "_coarse_cl.root";
	else outputFileName += outFile;
	TFile *out = new TFile(outputFileName, "RECREATE");
	TTree *tree = new TTree("h101", "h101");
	
	for(int i=1; i<=ALPIDE_NUM; ++i) {
		tree->Branch(Form("ALPIDE%dClust", i), &fitVector[i]);
		tree->Branch(Form("ALPIDE%dM", i), &clustM[i];	
	}
	
	auto t1 = timeNow();
	SetAllBranchAddress(h101);

	ulong lastEvent = sortEntries(firstEvent, maxEvents, h101);
	printf("Entries in file: %d\n", h101->GetEntries());
    ulong evCounter(0);
	
    for(ulong evNum = firstEvent; evNum < lastEvent; ++evNum) {
        ++evCounter;
        if(evCounter%100 == 0) {
            cout << "\rEv#: " << evCounter << " Completed: " << Form("%.1f",(float)evCounter/maxEvents*100) <<"%"<< flush;
        }
		h101->GetEntry(evNum);
		bool hasCluster(false);
		
		/* Part I : Construct clusters from raw data */
		for(int i(1); i<=ALPIDE_NUM; ++i) {
			if(rowM[i]>1 && rowM[i]==colM[i]) {
				clusters[i] = ConstructClusters(Col[i], Row[i], rowM[i], veto);
				clustM[i] = clusters[i].size();
				if(clustM.size() == 0) continue;
			
				/* Part II : Fit and feed */ 
				for(auto& cluster: clusters[i]) {
					fitVector.push_back(FitCluster(cluster));
				}
				hasCluster = true;
			}
		}
		if(hasCluster) tree->Fill();
		
		/* clear the fit vector for the next event */
		for_each(fitVector, fitVector+ALPIDE_NUM, [](auto& vec){vec.clear();});
	}

	tree->Write();
    auto t2 = timeNow();
    cout << "\nTime taken: " << duration_cast<seconds>(t2-t1).count() << "s\n";
}

auto main(int argc, char* argv[]) -> int {
    if(IsHelpArg(argc, argv)) {cout << clusterise_help; return 0;}
	
	string pStr;
    string fileName("");
	string outFile("");
    ulong firstEvent(0);
    ulong maxEvents(0);
    int veto(0);
    
    if(!ParseCmdLine("file", fileName, argc, argv)) {
		cerr << "No file specified!\n";
		cout << clusterise_help; return 0;
	}
	if(!ParseCmdLine("output", outFile, argc, argv)) {
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
            printf("Starting from ev#: %d\n", firstEvent);
        }
        catch(exception& e) {}
    }
    
    if(ParseCmdLine("max-events", pStr, argc, argv)) {
        try {
            maxEvents = stoul(pStr);
            printf("Max events: %d\n", maxEvents);
        }
        catch(exception& e) {}
    }

	cout<<endl; return 0;
}

/* aux functions */
bool IsHelpArg(int argc, char** argv) {
    for(int i(1); i<argc; ++i) {
        if(!strcmp(argv[i], "help") || !strcmp(argv[i], "--help")) return 1;
    }   
    return 0;
}

bool ParseCmdLine(const char* line, string& parsed, int argc, char** argv) {
    cmatch m;
    std::regex r("^(--|)([^=]+)[=](.+)$");
    for(int i(1); i<argc; ++i) {
        if(regex_match(argv[i], m, r) && !strcmp(m[2].str().c_str(), line)) {
            parsed = m[3].str(); return 1;
        }   
    }   
    return 0;
}   

void printProgress(double percentage) {
    int val = (int)(percentage*100);
    int lpad =(int)(percentage*BARW);
    int rpad = BARW-lpad;
    printf("\r%3d%% [%.*s%*s]",val,lpad,BAR,rpad,"");
    fflush(stdout);
}

