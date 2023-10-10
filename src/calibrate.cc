#include "AlpideClustering.h"
#include "CMDLineParser.h"
#include "AuxFunctions.h"
#include "libs.hh"
#include "spec.hh"

extern const std::string calibrate_help; 

using namespace std;
using namespace AlpideClustering;
using namespace AlpideAuxFunctions;

/* Read containers */
uint cNum;
uint AlpideID[MAX_CLUSTERS];
uint cSize[MAX_CLUSTERS];
double uCol[MAX_CLUSTERS];
double uRow[MAX_CLUSTERS];

void SetAllBranchAddress(TTree* h101, uint& cNum, uint* AlpideID, uint* cSize, double* uCol, double* uRow) {
	if(!h101 || h101->IsZombie()) return;
	h101->SetBranchAddress("CL_NUM", &cNum);
	h101->SetBranchAddress("ALPIDE_ID", AlpideID);
	h101->SetBranchAddress("CL_uCOL", uCol);
	h101->SetBranchAddress("CL_uROW", uRow);
}

void Calibrate(const char* fileName, const char* outFile, ulong firstEvent=0, ulong maxEvents=0) {
	/* It's better if file contains info from all the possible alpides */ 

	TFile* in = new TFile(fileName,"READ");
	if(!in || in->IsZombie()) {cerr << "Can't open rootfile with name: " << fileName << "\n"; exit(EXIT_FAILURE);}
    TApplication* app = new TApplication("myApp", 0, 0);
	TTree* h101 = dynamic_cast<TTree*>(in->Get("h101"));

	TFile* out = new TFile(outFile, "RECREATE");
	
	auto t1 = timeNow();
	TH2D* HitMapC[ALPIDE_NUM+1];
	TH2D* HitMapR[ALPIDE_NUM+1];
	TH1D* HitMapPC[ALPIDE_NUM+1][1024];
	TH1D* HitMapPR[ALPIDE_NUM+1][512];

	TGraph* colCal[ALPIDE_NUM+1];
	TGraph* colCalRes[ALPIDE_NUM+1];
	TGraph* rowCal[ALPIDE_NUM+1];
	TGraph* rowCalRes[ALPIDE_NUM+1];

	TF1* rowLine;
	TF1* colLine;

	for(int i=1; i<=ALPIDE_NUM; i++) {
		HitMapC[i] = new TH2D(TString::Format("HitMapC%d",i),TString::Format("HitMapC%d",i),1024,0,1024,1024,0,1024);
		HitMapR[i] = new TH2D(TString::Format("HitMapR%d",i),TString::Format("HitMapR%d",i),512,0,512,512,0,512);
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
	colLine = new TF1("colLine","[0]*x+[1]",0,1024); 
	rowLine = new TF1("rowLine","[0]*x+[1]",0,512); 

	SetAllBranchAddress(h101, cNum, AlpideID, cSize, uCol, uRow);

	ulong lastEvent = SortEntries(firstEvent, maxEvents, h101);
	printf("Entries in file: %lld\n", h101->GetEntries());
	ulong evCounter{0};
	/* MARK: event loop */	
	for(ulong evNum = firstEvent; evNum < lastEvent; ++evNum) {
		++evCounter; if(evCounter%100 == 0) PrintProgress((double)evCounter/maxEvents);		
		h101->GetEntry(evNum);

		for(uint i=0; i<cNum; ++i) {
			if(AlpideID[i] != 1) break; // id[i] == 1 must-have
			for(uint j=i+1; j<cNum; ++j) {
				int id = AlpideID[j];
				if(id == 1) continue; //don't correlate 1st vs. 1st
				HitMapC[id]->Fill(uCol[i], uCol[j]);
				HitMapR[id]->Fill(uRow[i], uRow[j]);
			}
		}
	}

	auto t2 = timeNow();
	cout << "\nTime taken for looping: " << duration_cast<seconds>(t2-t1).count() << "s\n";
	cout << "Doing analysis now, hol'up.. \n";
	double slopeCol[ALPIDE_NUM+1];
	double slopeColFine[ALPIDE_NUM+1]={0}; slopeColFine[1] = 1;
	double slopeColFineSig[ALPIDE_NUM+1] = {0};
	double slopeRow[ALPIDE_NUM+1];
	double slopeRowFine[ALPIDE_NUM+1]={0}; slopeRowFine[1] = 1;
	double slopeRowFineSig[ALPIDE_NUM+1] = {0};

	double offsetCol[ALPIDE_NUM+1] = {0};
	double offsetColFine[ALPIDE_NUM+1] = {0};
	double offsetColFineSig[ALPIDE_NUM+1] = {0};
	double offsetRow[ALPIDE_NUM+1] = {0};
	double offsetRowFine[ALPIDE_NUM+1] = {0};
	double offsetRowFineSig[ALPIDE_NUM+1] = {0};

	/* MARK: slicing and fitting */
	for(int det=2; det<=ALPIDE_NUM; ++det) {
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
		colCal[det]->SetNameTitle(TString::Format("ColCL%d",det));
		colCal[det]->SetTitle(TString::Format("ColCL%d",det));
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
			colCalRes[det]->AddPoint(col, (gradCol*col + intCol - maxBin)); //fit residues
			HitMapPC[det][col]->Delete();
		}
		
		// ---------------- //
		// ----- ROWS ----- //
		// ---------------- //
		for(int row=0; row<512; ++row) {
			HitMapPR[det][row] = (TH1D*)HitMapR[det]->ProjectionX(TString::Format("det%d; row%d",det,row), row, row+1);
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
		rowCal[det]->SetNameTitle(TString::Format("RowCL%d",det));
		rowCal[det]->SetTitle(TString::Format("RowCL%d",det));
		/* Do a finer calibration ... */
		for(int row=0; row<512; ++row) {
			HitMapPR[det][row] = (TH1D*)HitMapR[det]->ProjectionX(TString::Format("det%d; row%d",det,row), row, row+1);
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
			rowCalRes[det]->AddPoint(row,(gradRow*row + intRow - maxBin)); //fit residues  
			HitMapPR[det][row]->Delete();
		}
	}
	/* MARK: Writing output & plotting */
	cout << std::fixed << std::setprecision(5);
	TCanvas* c[ALPIDE_NUM+1];
	
	printf("\n----- %s -----\n", fileName);
	for(int det=2; det<=ALPIDE_NUM; ++det) {
		colCal[det]->Write();
		rowCal[det]->Write();
		/* colCalRes[det]->Write(); */
		/* rowCalRes[det]->Write(); */
		
		c[det] = new TCanvas(TString::Format("ALPIDE%d vs. ref", det), TString::Format("ALPIDE%d vs. ref", det), 1024, 1024);
		c[det]->Divide(1,2);
		c[det]->cd(1); colCal[det]->Draw();
		c[det]->cd(2); rowCal[det]->Draw();
		printf("Det %d :: Col coeffs :: Slope: %.5f +- %.5f | Offset: %.5f +- %.5f\n", det,slopeColFine[det],slopeColFineSig[det],offsetColFine[det],offsetColFineSig[det]);
		printf("         Row coeffs :: Slope: %.5f +- %.5f | Offset: %.5f +- %.5f",slopeRowFine[det],slopeRowFineSig[det],offsetRowFine[det],offsetRowFineSig[det]);
		printf("\n");
	}
	/* TVectorD objects to write to the output ROOT file */
	/* zero-th index should be filled with 0's or NAN's */
	/* arr[1] = 1. to each array as a 'reminder' that ALPIDE1 is the reference. */
	TVectorD* aCol = new TVectorD(ALPIDE_NUM+1,	 slopeColFine);	    aCol->Write("aCol");
	TVectorD* aColSig = new TVectorD(ALPIDE_NUM+1, slopeColFineSig);  aColSig->Write("aColSig");
	TVectorD* bCol = new TVectorD(ALPIDE_NUM+1,    offsetColFine);	bCol->Write("bCol");
	TVectorD* bColSig = new TVectorD(ALPIDE_NUM+1, offsetColFineSig); bColSig->Write("bColSig");
	
	TVectorD* aRow = new TVectorD(ALPIDE_NUM+1,    slopeRowFine);		aRow->Write("aRow");
	TVectorD* aRowSig = new TVectorD(ALPIDE_NUM+1, slopeRowFineSig);  aRowSig->Write("aRowSig");
	TVectorD* bRow = new TVectorD(ALPIDE_NUM+1,    offsetRowFine);	bRow->Write("bRow");
	TVectorD* bRowSig = new TVectorD(ALPIDE_NUM+1, offsetRowFineSig); bRowSig->Write("bRowSig");	

	out->Close();
	
	auto t3 = timeNow();
	cout << "\nTotal runtime: " << duration_cast<seconds>(t3-t1).count() << "s\n";
	
	app->Run();
	in->Close();
}

auto main(int argc, char* argv[]) -> int {
    if(IsCmdArg("help", argc, argv)) {cout << calibrate_help; return 0;}
	
	string pStr;
    string fileName;
	string outFile;
    ulong firstEvent(0);
    ulong maxEvents(0);
    
    if(!ParseCmdLine("file", fileName, argc, argv)) {
		cerr << "No file specified!\n";
		cout << clusterise_help; return 0;
	}

    if(!ParseCmdLine("output", outFile, argc, argv)) {
		outFile = fileName.substr(0, fileName.find('.')) + "_calib.root";
		cout << "No output file specified. Writing into file: " << outFile << endl;
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
	
	Calibrate(fileName.c_str(), outFile.c_str(), firstEvent, maxEvents);	

	cout<<endl; return 0;
}
