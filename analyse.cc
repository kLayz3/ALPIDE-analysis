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

/* Read containers */
uint cNum;
uint AlpideID[MAX_CLUSTERS];
uint cSize[MAX_CLUSTERS];
float uCol[MAX_CLUSTERS];
float uRow[MAX_CLUSTERS];

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
	
void CorrelationAll(const char* fileName, const char* outFile, ulong firstEvent=0, ulong maxEvents=0) {
	/* File better contain info from all the possible alpides */ 
	TFile* in = new TFile(fileName,"READ");
	if(!in || in->IsZombie()) {cerr << "Can't open rootfile with name: " << fileName << "\n"; exit(EXIT_FAILURE);}
	TTree* h101 = static_cast<TTree*>(in->Get("h101"));

	TFile* out = new TFile(outFile, "RECREATE");
    TTree* h102 = new TTree("h102","h102");
	
	auto t1 = timeNow();
	
	TH2D* HitMapR[ALPIDE_NUM-1];
	TH2D* HitMapC[ALPIDE_NUM-1];
	TH1D* HitMapPR[10][1024]; 
	TH1D* HitMapPC[10][512];
	TGraph* colCal[ALPIDE_NUM+1];
	TGraph* colCalRes[ALPIDE_NUM+1];
	TGraph* rowCal[ALPIDE_NUM+1];
	TGraph* rowCalRes[ALPIDE_NUM+1];
	TF1* rowLine;
	TF1* colLine;
	
	const double posZ[ALPIDE_NUM+1] = {0.,2.5, 5.0, 7.5, 10.0, 12.5};
	const double colC = 29.0/10000.0;
	const double rowC = 27.0/10000.0;

	SetAllBranchAddress(h101, cNum, AlpideID, cSize, uCol, uRow);

	ulong lastEvent = SortEntries(firstEvent, maxEvents, h101);
	printf("Entries in file: %lld\n", h101->GetEntries());
	ulong evCounter{0};

	for(ulong evNum = firstEvent; evNum < lastEvent; ++evNum) {
		++evCounter; if(evCounter%100 == 0) PrintProgress((float)evCounter/maxEvents);		
		h101->GetEntry(evNum);

		for(int i=0; i<cNum; ++i) {
			if(AlpideID[i] != 1) break; //only correlate vs positions in 1st det
			for(int j=i+1; j<cNum; ++j) {
				int id = AlpideID[j];
				if(id == 1) continue; //only correlate vs. other detectors
				HitMapC[id]->Fill(uCol[i], uCol[j]);
				HitMapR[id]->Fill(uRow[i], uRow[j]);
			}
		}
	}
	auto t2 = timeNow();
    cout << "\nTime taken for looping: " << duration_cast<seconds>(t2-t1).count() << "s\n";
	cout << "Doing analysis now, hol'up.. \n";
	for(int det=2; det<=ALPIDE_NUM; ++det) {
		HitMapC[det]->Write();
		HitMapR[det]->Write();
		/* COLUMNS */
		cout << "ALPIDE" << det << " projections.\n";
		for(int col=0; col<1024; ++col) {
			PrintProgress((float)col/1024);
			HitMapPC[det][col] = (TH1D*)HitMapC[det]->ProjectionX(TString::Format("column %d %d", det, col), col, col+1);
			int temp = HitMapPC[det][col]->GetMaximumBin();
			colCal[det]->SetPoint(col,col,temp);
			HitMapPC[det][col]->Delete();
		}
        colCal[det]->Fit(colLine);
        double gradCol = colLine->GetParameter(0);
        double intCol = colLine->GetParameter(1);
		for(int col = 0;col<1024;col++) {
			PrintProgress((float)col/1024);
			HitMapPC[det][col] = (TH1D*)HitMapC[det]->ProjectionX(TString::Format("column %d %d",det, col), col, col+1);
			int temp = HitMapPC[det][col]->GetMaximumBin();
			colCalRes[det]->SetPoint(col,col,(gradCol*col + intCol - temp));
			HitMapPC[det][col]->Delete();
		}
		/* ROWS */
		for(int row=0; row<512; ++row) {
			PrintProgress((float)row/512);	
			HitMapPR[det][row] = (TH1D*)HitMapC[det]->ProjectionX(TString::Format("row %d %d", det, row), row, row+1);
			int temp = HitMapPR[det][row]->GetMaximumBin();
			rowCal[det]->SetPoint(row,row,temp);
			HitMapPR[det][row]->Delete();
		}
		rowCal[det]->Fit(rowLine);
		double gradRow = rowLine->GetParameter(0);
		double intRow = rowLine->GetParameter(1);
		for(int row=0;row<512;row++) {
			PrintProgress((float)row/512);	
			HitMapPR[det][row] = (TH1D*)HitMapR[det]->ProjectionX(TString::Format("column %d %d",det, row), row, row+1);
			int temp = HitMapPR[det][row]->GetMaximumBin();
			rowCalRes[det]->SetPoint(row,row,(gradRow*row + intRow - temp));  
			HitMapPR[det][row]->Delete();
		} 
	}

	for(int det=2; det<=ALPIDE_NUM; ++det) {
		colCal[det]->Write();
		rowCal[det]->Write();
		colCalRes[det]->Write();
		rowCalRes[det]->Write();
	}
	out->Write();	
	out->Close();
	in->Close();
}


auto main(int argc, char* argv[]) -> int {
    if(IsCmdArg("help", argc, argv)) {cout << analyse_help; return 0;}
    
    string pStr;
    string fileName;
	string outFile;
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
	if(!ParseCmdLine("output", outFile, argc, argv)) {
		outFile = fileName.substr(0, fileName.find('.')) + "_ANALYSIS.root";
		cout << "No output file specified. Writing into file: " << outFile << endl;
	}
	if(IsCmdArg("cal", argc, argv)) {
		CorrelationAll(fileName.c_str(), outFile.c_str(), firstEvent, maxEvents);	
	}
	
    cout<<endl; return 0;
}


