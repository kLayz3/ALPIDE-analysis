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
	
void CalibrateAll(const char* fileName, const char* outFile, ulong firstEvent=0, ulong maxEvents=0) {
	/* File better contain info from all the possible alpides */ 
	TFile* in = new TFile(fileName,"READ");
	if(!in || in->IsZombie()) {cerr << "Can't open rootfile with name: " << fileName << "\n"; exit(EXIT_FAILURE);}
    TApplication* app = new TApplication("myApp", 0, 0);
	TTree* h101 = static_cast<TTree*>(in->Get("h101"));

	TFile* out = new TFile(outFile, "RECREATE");
    TTree* h102 = new TTree("h102","h102");
	
	auto t1 = timeNow();
	TH2D* HitMapR[ALPIDE_NUM+1];
	TH2D* HitMapC[ALPIDE_NUM+1];
	TH1D* HitMapPR[ALPIDE_NUM+1][1024];
	TH1D* HitMapPC[ALPIDE_NUM+1][512];
	TGraph* colCal[ALPIDE_NUM+1];
	TGraph* colCalRes[ALPIDE_NUM+1];
	TGraph* rowCal[ALPIDE_NUM+1];
	TGraph* rowCalRes[ALPIDE_NUM+1];
	TF1* rowLine;
	TF1* colLine;
	for(int i=1; i<=ALPIDE_NUM; i++) {
		HitMapR[i] = new TH2D(TString::Format("HitMapR%d",i),TString::Format("HitMapR%d",i),512,0,512,512,0,512);
		HitMapC[i] = new TH2D(TString::Format("HitMapC%d",i),TString::Format("HitMapC%d",i),1024,0,1024,1024,0,1024);
		colCal[i] = new TGraph();
		rowCal[i] = new TGraph();
		colCal[i]->SetNameTitle(TString::Format("colCal%i",i));
		rowCal[i]->SetNameTitle(TString::Format("rowCal%i",i));
		colCal[i]->SetTitle(TString::Format("colCal%i",i));
		rowCal[i]->SetTitle(TString::Format("rowCal%i",i));
		colCalRes[i] = new TGraph();
		rowCalRes[i] = new TGraph();
		colCalRes[i]->SetNameTitle(TString::Format("colCalRes%i",i));
		rowCalRes[i]->SetNameTitle(TString::Format("rowCalRes%i",i));
		colCalRes[i]->SetTitle(TString::Format("colCalRes%i",i));
		rowCalRes[i]->SetTitle(TString::Format("rowCalRes%i",i));
	}
	rowLine = new TF1("rowLine","[0]*x+[1]",0,512); 
	colLine = new TF1("colLine","[0]*x+[1]",0,1024); 

	/* const double posZ[ALPIDE_NUM+1] = {0.0, 0., 2.5, 5.0, 7.5, 10.0, 12.5}; */
	/* const double colC = 29.0/10000.0; */
	/* const double rowC = 27.0/10000.0; */

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
				if(id == 1) continue; //don't correlate 1st vs. 1st
				auto hitIntensity = cSize[i] * cSize[j];
				/* HitMapC[id]->Fill(uCol[i], uCol[j], hitIntensity); */
				/* HitMapR[id]->Fill(uRow[i], uRow[j], hitIntensity); */

				HitMapC[id]->Fill(uCol[i], uCol[j]);
				HitMapR[id]->Fill(uRow[i], uRow[j]);
			}
		}
	}

	auto t2 = timeNow();
    cout << "\nTime taken for looping: " << duration_cast<seconds>(t2-t1).count() << "s\n";
	cout << "Doing analysis now, hol'up.. \n";
	double slopeCol[ALPIDE_NUM];
	double slopeColFine[ALPIDE_NUM+1];
	double slopeColFineSig[ALPIDE_NUM+1];
	double slopeRow[ALPIDE_NUM+1];
	double slopeRowFine[ALPIDE_NUM+1];
	double slopeRowFineSig[ALPIDE_NUM+1];

	double offsetCol[ALPIDE_NUM];
	double offsetColFine[ALPIDE_NUM+1];
	double offsetColFineSig[ALPIDE_NUM+1];
	double offsetRow[ALPIDE_NUM+1];
	double offsetRowFine[ALPIDE_NUM+1];
	double offsetRowFineSig[ALPIDE_NUM+1];
	for(int det=2; det<=ALPIDE_NUM; ++det) {
		HitMapC[det]->Write();
		HitMapR[det]->Write();
		
		// -------------------	
		// ----- COLUMNS -----
		// -------------------	
		cout << "\n\n----- ALPIDE" << det << " projections -----\n";
		for(int col=0; col<1024; ++col) {
			HitMapPC[det][col] = (TH1D*)HitMapC[det]->ProjectionX(TString::Format("det%d; col%d",det,col), col, col+1);
			int maxBin = HitMapPC[det][col]->GetMaximumBin();
			colCal[det]->AddPoint(col,maxBin);
			HitMapPC[det][col]->Delete();
		}
		cout << "\nCol: after initial calibration: \n";
        colCal[det]->Fit(colLine);
        double gradCol = colLine->GetParameter(0);
        double intCol = colLine->GetParameter(1);
		
		slopeCol[det] = gradCol;
		offsetCol[det] = intCol;

		colCal[det]->Delete();
		colCal[det] = new TGraph();
		/* Do a finer calibration ... */
		for(int col=0; col<1024; ++col) {
			HitMapPC[det][col] = (TH1D*)HitMapC[det]->ProjectionX(TString::Format("det%d; col%d",det,col), col, col+1);
			int maxBin = HitMapPC[det][col]->GetMaximumBin();
			if(abs(gradCol*col + intCol - maxBin) < 20) colCal[det]->AddPoint(col, maxBin);	
			HitMapPC[det][col]->Delete();
		}
		cout << "\nCol: after fine calibration: \n";
        colCal[det]->Fit(colLine);
        gradCol = colLine->GetParameter(0);
        intCol = colLine->GetParameter(1);

		slopeColFine[det] = gradCol;
		slopeColFineSig[det] = colLine->GetParError(0);
		offsetColFine[det] = intCol;
		offsetColFineSig[det] = colLine->GetParError(1);
		
		for(int col=0; col<1024; ++col) {
			HitMapPC[det][col] = (TH1D*)HitMapC[det]->ProjectionX(TString::Format("det%d; col%d",det,col), col, col+1);
			int maxBin = HitMapPC[det][col]->GetMaximumBin();
			colCalRes[det]->AddPoint(col, (gradCol*col + intCol - maxBin));
			HitMapPC[det][col]->Delete();
		}
		
		// ----------------
		// ----- ROWS -----
		// ----------------
		for(int row=0; row<512; ++row) {
			HitMapPR[det][row] = (TH1D*)HitMapC[det]->ProjectionX(TString::Format("det%d; row%d",det,row), row, row+1);
			int maxBin = HitMapPR[det][row]->GetMaximumBin();
			rowCal[det]->AddPoint(row,maxBin);
			HitMapPR[det][row]->Delete();
		}
		cout << "\nRow: after initial calibration: \n";
		rowCal[det]->Fit(rowLine);
		double gradRow = rowLine->GetParameter(0);
		double intRow = rowLine->GetParameter(1);

		slopeRow[det] = gradRow;
		offsetRow[det] = intRow;

		rowCal[det]->Delete();
		rowCal[det] = new TGraph();
		/* Do a finer calibration ... */
		for(int row=0; row<1024; ++row) {
			HitMapPR[det][row] = (TH1D*)HitMapC[det]->ProjectionX(TString::Format("det%d; row%d",det,row), row, row+1);
			int maxBin = HitMapPR[det][row]->GetMaximumBin();
			if(abs(gradRow*row + intRow - maxBin) < 20) rowCal[det]->AddPoint(row, maxBin);
			HitMapPR[det][row]->Delete();
		}
		cout << "\nRow: after fine calibration: \n";
        rowCal[det]->Fit(rowLine);
        gradRow = rowLine->GetParameter(0);
        intRow = rowLine->GetParameter(1);

		slopeRowFine[det] = gradRow;
		slopeRowFineSig[det] = rowLine->GetParError(0);
		offsetRowFine[det] = intRow;
		offsetRowFineSig[det] = rowLine->GetParError(1);

		for(int row=0;row<512;row++) {
			HitMapPR[det][row] = (TH1D*)HitMapR[det]->ProjectionX(TString::Format("det%d; row%d",det,row), row, row+1);
			int maxBin = HitMapPR[det][row]->GetMaximumBin();
			rowCalRes[det]->AddPoint(row,(gradRow*row + intRow - maxBin));  
			HitMapPR[det][row]->Delete();
		}
	}
	cout << std::fixed << std::setprecision(5);
	TCanvas* c[ALPIDE_NUM+1];
	for(int det=2; det<=ALPIDE_NUM; ++det) {
		colCal[det]->Write();
		rowCal[det]->Write();
		colCalRes[det]->Write();
		rowCalRes[det]->Write();
		
		c[det] = new TCanvas(TString::Format("ALPIDE%d vs. ref", det), TString::Format("ALPIDE%d vs. ref", det), 1024, 1024);
		c[det]->Divide(1,2);
		c[det]->cd(1); colCal[det]->Draw();
		c[det]->cd(2); rowCal[det]->Draw();
		printf("Det %d :: Col coeffs :: Slope: %.5f +- %.5f | Offset: %.5f +- %.5f\n", det,slopeColFine[det],slopeColFineSig[det],offsetColFine[det],offsetColFineSig[det]);
		printf("         Row coeffs :: Slope: %.5f +- %.5f | Offset: %.5f +- %.5f",slopeRowFine[det],slopeRowFineSig[det],offsetRowFine[det],offsetRowFineSig[det]);
		printf("\n");
	}
	out->Write();	
	out->Close();	
	app->Run();
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
		CalibrateAll(fileName.c_str(), outFile.c_str(), firstEvent, maxEvents);	
	}
	
    cout<<endl; return 0;
}


