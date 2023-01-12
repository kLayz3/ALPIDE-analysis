#include "AlpideClustering.h"
#include "AlpideTrack.h"
#include "CMDLineParser.h"
#include "AuxFunctions.h"
#include "libs.hh"
#include "spec.hh"

#define LEN(x) (sizeof x / sizeof *x)
#define timeNow() std::chrono::high_resolution_clock::now()
#define MAX_TRACKS 128

typedef uint32_t uint;
typedef uint64_t ulong;
extern const std::string analyse_help; 

using namespace std;
using namespace AlpideClustering;
using std::chrono::duration_cast;
using std::chrono::seconds;

void SetOneBranchAddressRaw(TTree* h101, int x, uint* Col, uint* Row, uint& colM, uint& rowM, uint& tHi, uint& tLo) {
    assert(x<=ALPIDE_NUM && x>=1);
    h101->SetBranchAddress(TString::Format("ALPIDE%dCOLv", x), Col);
    h101->SetBranchAddress(TString::Format("ALPIDE%dROWv", x), Row);
    h101->SetBranchAddress(TString::Format("ALPIDE%dCOL", x), &colM);
    h101->SetBranchAddress(TString::Format("ALPIDE%dROW", x), &rowM);
    h101->SetBranchAddress(TString::Format("ALPIDE%dT_HI", x), &tHi);
    h101->SetBranchAddress(TString::Format("ALPIDE%dT_LO", x), &tLo);
}

	void SetAllBranchAddressRaw(TTree* h101, uint (*Col)[MAX_HITS], uint (*Row)[MAX_HITS], uint* colM, uint* rowM, uint& tHi, uint& tLo) {
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

	void SetAllBranchAddressClust(TTree* h101, uint& cNum, uint* AlpideID, uint* cSize, float* uCol, float* uRow, uint& tHi, uint& tLo) {
		if(!h101 || h101->IsZombie()) return;
		h101->SetBranchAddress("T_HI", &tHi);
		h101->SetBranchAddress("T_LO", &tLo);
		h101->SetBranchAddress("CL_NUM", &cNum);
		h101->SetBranchAddress("ALPIDE_ID", AlpideID);
		h101->SetBranchAddress("CL_SIZE", cSize);
		h101->SetBranchAddress("CL_uCOL", uCol);
		h101->SetBranchAddress("CL_uROW", uRow);
	}

	void RawHitMap(const char*, int, ulong, ulong);
	void RawHitMapAll(const char*, ulong, ulong);
	void RawCorrelation(const char*,  int, int, ulong, ulong);
	void ClustHitMap(const char*, int, ulong, ulong);
	void ClustHitMapAll(const char*, int, ulong, ulong);
	void Hitmap(const char*, int, ulong, ulong);
	void HitmapAll(const char*, ulong, ulong);
	void TrackAnalysis(const char*, const char*, ulong, ulong, const char*);

	/* MARK: main */
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
		
		if(IsCmdArg("hitmap", argc, argv)) {
			HitmapAll(fileName.c_str(), firstEvent, maxEvents);
		} 
		else if(ParseCmdLine("hitmap", pStr, argc, argv)) {
			try {
				int x = stoi(pStr);
				Hitmap(fileName.c_str(), x, firstEvent, maxEvents);
			}
			catch(exception& e) {}
		}

		if(ParseCmdLine("corr", pStr, argc, argv)) {
			auto split = SplitStringToVector(pStr, ',');
			try {
				int x = stoi(split[0]);
				int y = stoi(split[1]);
				RawCorrelation(fileName.c_str(), x,y, firstEvent, maxEvents);
			}
			catch(exception& e) {}
		}
		
		if(IsCmdArg("track", argc, argv)) {
			string calFile; string saveFile;
			if(!ParseCmdLine("cal", calFile, argc, argv)) {
				cerr << __PRETTY_FUNCTION__ << " -> provide calibration file with --cal=FILENAME.root"; return -1;
			}

			if(ParseCmdLine("save", saveFile, argc, argv))
				TrackAnalysis(fileName.c_str(), calFile.c_str(), firstEvent, maxEvents, saveFile.c_str()); // saves tracks to a file	
			else
				TrackAnalysis(fileName.c_str(), calFile.c_str(), firstEvent, maxEvents, nullptr); // doesn't save tracks to a file 
		}
		cout<<endl; return 0;
	}

	void RawHitMap(const char* fileName, int x, ulong firstEvent=0, ulong maxEvents=0) {
		assert(x>=1 && x<=ALPIDE_NUM);
		TApplication* app = new TApplication("myApp", 0, 0);
		TFile* in = new TFile(fileName,"READ");
		if(!in || in->IsZombie()) {cerr << "Can't open rootfile with name: " << fileName << "\n"; return;}
		TTree* h101 = dynamic_cast<TTree*>(in->Get("h101"));
		auto t1 = timeNow();
		
		/* Read Containers */
		constexpr size_t MALLOC_SIZE = MAX_HITS * sizeof(uint);
		uint colM = 0; uint rowM = 0; 
		uint tHi, tLo;
		uint* colV = (uint*)malloc(MALLOC_SIZE);
		uint* rowV = (uint*)malloc(MALLOC_SIZE);
		SetOneBranchAddressRaw(h101, x, colV, rowV, colM, rowM, tHi, tLo);
		
		TH2D* hRawHit = new TH2D(TString::Format("Raw Hitmap ALPIDE%d", x), TString::Format("Raw Hitmap ALPIDE%d", x), 1024,0,1024,512,0,512);  

		ulong lastEvent = SortEntries(firstEvent, maxEvents, h101);
		printf("Entries in file: %lld\n", h101->GetEntries());
		ulong evCounter{0};
		
		for(ulong evNum = firstEvent; evNum < lastEvent; ++evNum) {
			++evCounter; if(evCounter%100 == 0) PrintProgress((float)evCounter/maxEvents);		
			h101->GetEntry(evNum);
			if(colM == 0 || colM != rowM) continue;
			for(int i=0; i<colM; ++i) {
				hRawHit->Fill(colV[i], rowV[i]);
			}
		}
		hRawHit->Draw("colz");
		
		auto t2 = timeNow();
		cout << "\nTime taken: " << duration_cast<seconds>(t2-t1).count() << "s\n";
		app->Run();

		in->Close();
	}

	void RawHitMapAll(const char* fileName, ulong firstEvent=0, ulong maxEvents=0) {
		TApplication* app = new TApplication("myApp", 0, 0);
		TFile* in = new TFile(fileName,"READ");
		if(!in || in->IsZombie()) {cerr << "Can't open rootfile with name: " << fileName << "\n"; return;}
		TTree* h101 = dynamic_cast<TTree*>(in->Get("h101"));
		auto t1 = timeNow();
		
		/* Read Containers */
		uint tHi, tLo;
		uint colV[ALPIDE_NUM+1][MAX_HITS];
		uint rowV[ALPIDE_NUM+1][MAX_HITS];
		uint colM[ALPIDE_NUM+1];
		uint rowM[ALPIDE_NUM+1];
		
		SetAllBranchAddressRaw(h101, colV, rowV, colM, rowM, tHi, tLo);
		
		TH2D* hRawHit[ALPIDE_NUM+1];
		TCanvas *canvas[ALPIDE_NUM+1];
		for(int x=1; x<=ALPIDE_NUM; ++x) {
			hRawHit[x]  = new TH2D(TString::Format("Raw Hitmap ALPIDE%d", x), TString::Format("Raw Hitmap ALPIDE%d", x), 1024,0,1024,512,0,512);  
			canvas[x] = new TCanvas(TString::Format("Raw Hitmap ALPIDE%d", x), TString::Format("Raw Hitmap ALPIDE%d", x), 1200,1200);
		}
		ulong lastEvent = SortEntries(firstEvent, maxEvents, h101);
		printf("Entries in file: %lld\n", h101->GetEntries());
		ulong evCounter{0};
		
		for(ulong evNum = firstEvent; evNum < lastEvent; ++evNum) {
			++evCounter; if(evCounter%100 == 0) PrintProgress((float)evCounter/maxEvents);		
			h101->GetEntry(evNum);
			for(int x=1; x<=ALPIDE_NUM; ++x) {
				if(colM[x] == 0 || colM[x] != rowM[x]) continue;
				for(int i=0; i<colM[x]; ++i) {
					hRawHit[x]->Fill(colV[x][i], rowV[x][i]);
				}
			}
		}
		for(int x=1; x<=ALPIDE_NUM; ++x) {
			canvas[x]->cd();
			hRawHit[x]->Draw("colz");
		}

		auto t2 = timeNow();
		cout << "\nTime taken: " << duration_cast<seconds>(t2-t1).count() << "s\n";

		app->Run();
		in->Close();
	}

	void RawCorrelation(const char* fileName, int x, int y, ulong firstEvent, ulong maxEvents) {
		assert(x>=1 && x<=ALPIDE_NUM && y>=1 && y<=ALPIDE_NUM && x!=y);
		TApplication* app = new TApplication("myApp", 0, 0);
		TFile* in = new TFile(fileName,"READ");
		if(!in || in->IsZombie()) {cerr << "Can't open rootfile with name: " << fileName << "\n"; return;}
		TTree* h101 = dynamic_cast<TTree*>(in->Get("h101"));
		auto t1 = timeNow();

		uint tHi[MAX_HITS], tLo[MAX_HITS];
		uint colV[ALPIDE_NUM+1][MAX_HITS];
		uint rowV[ALPIDE_NUM+1][MAX_HITS];
		uint colM[ALPIDE_NUM+1];
		uint rowM[ALPIDE_NUM+1];

		SetOneBranchAddressRaw(h101, x, colV[x], rowV[x], colM[x], rowM[x], tHi[x], tLo[x]);
		SetOneBranchAddressRaw(h101, y, colV[y], rowV[y], colM[y], rowM[y], tHi[y], tLo[y]);
		
		TH2D* hColCorr = new TH2D(TString::Format("Column corr ALPIDE%d.vs.%d", y, x), TString::Format("Column corr ALPIDE%d.vs.%d", y,x), 1024,0,1024,1024,0,1024);
		TH2D* hRowCorr = new TH2D(TString::Format("Row corr ALPIDE%d.vs.%d", y, x), TString::Format("Row corr ALPIDE%d.vs.%d", y,x), 512,0,512,512,0,512);
		
		ulong lastEvent = SortEntries(firstEvent, maxEvents, h101);
		printf("Entries in file: %lld\n", h101->GetEntries());
		ulong evCounter{0};
		
		for(ulong evNum = firstEvent; evNum < lastEvent; ++evNum) {
			++evCounter; if(evCounter%100 == 0) PrintProgress((float)evCounter/maxEvents);		
			h101->GetEntry(evNum);
			if(colM[x] == 0 || colM[y] == 0 || colM[x]!=rowM[x] || colM[y]!=rowM[y]) continue;
			for(int i=0; i<colM[x]; ++i) {
				for(int j=0; j<colM[y]; ++j) {
					hColCorr->Fill(colV[x][i], colV[y][j]);
					hRowCorr->Fill(rowV[x][i], rowV[y][j]);
				}
			}
		}

		TCanvas* canvas = new TCanvas(TString::Format("Raw corr ALPIDE%d.vs.%d", y,x), TString::Format("Raw corr ALPIDE%d.vs.%d", y,x), 1200,1600);
		canvas->Divide(1,2);
		canvas->cd(1); hColCorr->Draw("colz");
		canvas->cd(2); hRowCorr->Draw("colz");
		
		auto t2 = timeNow();
		cout << "\nTime taken: " << duration_cast<seconds>(t2-t1).count() << "s\n";

		app->Run();
		in->Close();
	}

	void ClustHitMap(const char* fileName, int x, ulong firstEvent=0, ulong maxEvents=0) {
		assert(x>=1 && x<=ALPIDE_NUM);
		TApplication* app = new TApplication("myApp", 0, 0);
		TFile* in = new TFile(fileName,"READ");
		if(!in || in->IsZombie()) {cerr << "Can't open rootfile with name: " << fileName << "\n"; return;}
		TTree* h101 = dynamic_cast<TTree*>(in->Get("h101"));
		auto t1 = timeNow();

		/* Read containers */
		uint cNum, tHi, tLo;
		uint AlpideID[MAX_CLUSTERS];
		uint cSize[MAX_CLUSTERS];
		float uCol[MAX_CLUSTERS];
		float uRow[MAX_CLUSTERS];
		SetAllBranchAddressClust(h101, cNum, AlpideID, cSize, uCol, uRow, tHi, tLo);

		TH2D* hRawHit = new TH2D(TString::Format("Clust Hitmap ALPIDE%d", x), TString::Format("Clust Hitmap ALPIDE%d", x), 1024,0,1024,512,0,512);  

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

	void ClustHitMapAll(const char* fileName, int x, ulong firstEvent=0, ulong maxEvents=0) {
		TApplication* app = new TApplication("myApp", 0, 0);
		TFile* in = new TFile(fileName,"READ");
		if(!in || in->IsZombie()) {cerr << "Can't open rootfile with name: " << fileName << "\n"; return;}
		TTree* h101 = dynamic_cast<TTree*>(in->Get("h101"));
		auto t1 = timeNow();

		/* Read containers */
		uint cNum, tHi, tLo;
		uint AlpideID[MAX_CLUSTERS];
		uint cSize[MAX_CLUSTERS];
		float uCol[MAX_CLUSTERS];
		float uRow[MAX_CLUSTERS];
		SetAllBranchAddressClust(h101, cNum, AlpideID, cSize, uCol, uRow, tHi, tLo);

		TH2D* hRawHit[ALPIDE_NUM+1];
		TCanvas *canvas[ALPIDE_NUM+1];

		for(int x=1; x<=ALPIDE_NUM; ++x) {
			hRawHit[x]  = new TH2D(TString::Format("Clust Hitmap ALPIDE%d", x), TString::Format("Clust Hitmap ALPIDE%d", x), 1024,0,1024,512,0,512);  
			canvas[x] = new TCanvas(TString::Format("Clust Hitmap ALPIDE%d", x), TString::Format("Clust Hitmap ALPIDE%d", x), 1200,1200);
		}
		ulong lastEvent = SortEntries(firstEvent, maxEvents, h101);
		printf("Entries in file: %lld\n", h101->GetEntries());
		ulong evCounter{0};
		
		for(ulong evNum = firstEvent; evNum < lastEvent; ++evNum) {
			++evCounter; if(evCounter%100 == 0) PrintProgress((float)evCounter/maxEvents);		
			h101->GetEntry(evNum);

			for(uint c=0; c<cNum; ++c) {
				int x = AlpideID[c];
				hRawHit[x]->Fill(uCol[c], uRow[c], (double)cSize[c]);
			}
		}
		
		for(int x=1; x<=ALPIDE_NUM; ++x) {
			canvas[x]->cd();
			hRawHit[x]->Draw("colz");
		}

		auto t2 = timeNow();
		cout << "\nTime taken: " << duration_cast<seconds>(t2-t1).count() << "s\n";
		app->Run();

		in->Close();
	}

	void Hitmap(const char* fileName, int x, ulong firstEvent, ulong maxEvents) {
		TFile* in = new TFile(fileName, "READ");
		if(!in || in->IsZombie()) {cerr << "Can't open rootfile with name: " << fileName << "\n"; return;}
		TTree* h101 = dynamic_cast<TTree*>(in->Get("h101"));

		auto br = h101->FindBranch("uCol");
		in->Close();
		
		if(br) ClustHitMap(fileName, x, firstEvent, maxEvents);
		else RawHitMap(fileName, x, firstEvent, maxEvents);
	}

	void HitmapAll(const char* fileName, ulong firstEvent, ulong maxEvents) {
		TFile* in = new TFile(fileName, "READ");
		if(!in || in->IsZombie()) {cerr << "Can't open rootfile with name: " << fileName << "\n"; return;}
		TTree* h101 = dynamic_cast<TTree*>(in->Get("h101"));
		auto br = (ulong)h101->FindBranch("CL_SIZE");
		in->Close();

		if(br) ClustHitMapAll(fileName, firstEvent, maxEvents);
		else RawHitMapAll(fileName, firstEvent, maxEvents);
	}
	/* MARK: tracks */
	void TrackAnalysis(const char* fileName, const char* calibFile, ulong firstEvent=0, ulong maxEvents=0, const char* saveFile=nullptr) {
		TApplication* app = new TApplication("myApp", 0, 0);
		TFile* in = new TFile(fileName,"READ");
		if(!in || in->IsZombie()) {cerr << "Can't open rootfile with name: " << fileName << "\n"; return;}
		
		AlpideTrack::LoadCal(calibFile); // loads the static fields for AlpideTrack class
		
		TTree* h101 = dynamic_cast<TTree*>(in->Get("h101"));
		auto t1 = timeNow();

		uint cNum, tHi, tLo;
		uint AlpideID[MAX_CLUSTERS];
		uint cSize[MAX_CLUSTERS];
		float uCol[MAX_CLUSTERS];
		float uRow[MAX_CLUSTERS];

		SetAllBranchAddressClust(h101, cNum, AlpideID, cSize, uCol, uRow, tHi, tLo);
		TH1D* hPhiX = new TH1D("PhiX ALPIDE", "PhiX ALPIDE", 300, -0.03, 0.03); 
		TH1D* hPhiY = new TH1D("PhiY ALPIDE", "PhiY ALPIDE",300, -0.03, 0.03);
		TH2D* hXY = new TH2D("Hitmap XY", "Hitmap XY", 500, -5, 35, 500, -5, 20);
		
		TFile* out = nullptr; TTree* tree = nullptr;
		/* Save containers */
		double x0,y0,phiX,phiY; // fit parameters
		int nDet; // number of detectors hit in a track
		double tZ[MAX_TRACKS];
		double tX[MAX_TRACKS];
		double tY[MAX_TRACKS];
		uint tSize[MAX_TRACKS];

		if(saveFile) {
			out = new TFile(saveFile, "RECREATE");
			if(out->IsZombie()) {cerr << __PRETTY_FUNCTION__ << " - can't open rootfile with name: " << saveFile << "\n"; return;}
			tree = new TTree("h101", "h101");
			tree->Branch("T_HI", &tHi);
			tree->Branch("T_LO", &tLo);
			tree->Branch("N_HITS", &nDet);
			tree->Branch("Z", tZ, "Z[N_HITS]/D");
			tree->Branch("X", tX, "X[N_HITS]/D");
			tree->Branch("Y", tY, "Y[N_HITS]/D");
			tree->Branch("CL_SIZE", tSize, "CL_SIZE[N_HITS]/i");
			tree->Branch("X0", &x0);
			tree->Branch("Y0", &y0);
			tree->Branch("PHI_X", &phiX);
			tree->Branch("PHI_Y", &phiY);
		}
		ulong lastEvent = SortEntries(firstEvent, maxEvents, h101);
		printf("Entries in file: %lld\n", h101->GetEntries());
		ulong evCounter{0};
		
		for(ulong evNum = firstEvent; evNum < lastEvent; ++evNum) {
			++evCounter; if(evCounter%100 == 0) PrintProgress((float)evCounter/maxEvents);		
			h101->GetEntry(evNum);
			if(cNum<3) continue;

			vector<int> detV(AlpideID, AlpideID + cNum);
			vector<float> colV(uCol, uCol + cNum);
			vector<float> rowV(uRow, uRow + cNum);
			vector<uint> sizeV(cSize, cSize + cNum);

		while(detV.size()) {
			AlpideTrack tr(detV, colV, rowV, sizeV);
			if(!tr.IsValidTrack()) continue;
			/* Track points (X,Y,Z,Size) are in vector fields of tr
			 * However it's convenient to work with c-arrays in root branches. */
			std::copy(tr.tZ.begin(), tr.tZ.end(), tZ);
			std::copy(tr.tX.begin(), tr.tX.end(), tX);
			std::copy(tr.tY.begin(), tr.tY.end(), tY);
			std::copy(tr.tSize.begin(), tr.tSize.end(), tSize);
			nDet = tr.GetSize();

			tr.GetFit(x0,y0,phiX,phiY); // saves the fit parameters in the arguments

			/* Filling hists and writing to output tree */
			hPhiX->Fill(phiX);
			hPhiY->Fill(phiY);
			hXY->Fill(x0,y0);
			if(saveFile) tree->Fill();
		}
	}

	TCanvas* ang = new TCanvas("Angular Distribution", "Angular Distribution", 1800,1200);
	ang->Divide(2,1);
	ang->cd(1); hPhiX->Draw();
	ang->cd(2); hPhiY->Draw();

	TCanvas* posXY = new TCanvas("Track hitmap in ALPIDE1", "Track hitmap in ALPIDE1", 1800,1200);
	hXY->Draw("colz");
	
	out->Write();
	out->Close();
	
	auto t2 = timeNow();
    cout << "\nTime taken: " << duration_cast<seconds>(t2-t1).count() << "s\n";
	app->Run();

	in->Close();
}



