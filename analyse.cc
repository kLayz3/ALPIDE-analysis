#include "includes/AlpideClustering.cpp"
#include "libs.hh"
#include "spec.hh"

#define LEN(x) (sizeof x / sizeof *x)
#define timeNow() std::chrono::high_resolution_clock::now()

typedef uint32_t uint;
typedef uint64_t ulong;
extern const std::string analyse_help;
extern const std::vector<std::string> maskFile; 

using namespace std;
using std::chrono::duration_cast;
using std::chrono::seconds;
using namespace AlpideClustering;

/* Global containers */
UInt_t Row[ALPIDE_NUM+1][MAX_HITS];
UInt_t Col[ALPIDE_NUM+1][MAX_HITS];
UInt_t rowM[ALPIDE_NUM+1];
UInt_t colM[ALPIDE_NUM+1];
/* Handle [0] doesn't represent a chip so is discarded */

const unordered_set<string> cmdLineArgs{"help", "file", "veto", "first-event", "max-events"};
/* vector<UInt_t*> colHandle{nullptr, &col1, &col2, &col3, &col4, &col5, &col6}; */
/* vector<UInt_t*> rowHandle{nullptr, &row1, &row2, &row3, &row4, &row5, &row6}; */
/* vector<UInt_t*> ColHandle{nullptr, Col1, Col2, Col3, Col4, Col5, Col6}; */
/* vector<UInt_t*> RowHandle{nullptr, Row1, Row2, Row3, Row4, Row5, Row6}; */

bool IsHelpArg(int, char**);
bool ParseCmdLine(const char*, std::string&, int, char**);
void printProgress(double);

void SetOneBranchAddress(TTree* h101, int x) {
    assert(x<=ALPIDE_NUM && x>=1); 
    h101->SetBranchAddress(Form("ALPIDE%dCOLv", x), Col[x]);
    h101->SetBranchAddress(Form("ALPIDE%dROWv", x), Row[x]);
    h101->SetBranchAddress(Form("ALPIDE%dCOL", x), colM[x]);
    h101->SetBranchAddress(Form("ALPIDE%dROW", x), rowM[x]);
}

void SetAllBranchAddress(TTree* h101) {
	if(!h101) return;
	for(int x=1; x<=ALPIDE_NUM; ++x) {
		h101->SetBranchAddress(Form("ALPIDE%dCOLv", x), Col[x]);
		h101->SetBranchAddress(Form("ALPIDE%dROWv", x), Row[x]);
		h101->SetBranchAddress(Form("ALPIDE%dCOL", x), colM[x]);
		h101->SetBranchAddress(Form("ALPIDE%dROW", x), rowM[x]);
	}
} 

void RawHitMap(const char* fileName, int x, ulong firstEvent=0, ulong maxEvents=0) {
    assert(x<=6 && x>=1);

    TApplication* app = new TApplication("myApp", 0, 0);
    TFile* in = new TFile(fileName,"READ");
    if(!in || in->IsZombie()) {cerr << "Can't open rootfile with name: " << fileName << "\n"; exit(EXIT_FAILURE);}
    TTree* h101 = (TTree*)in->Get("h101");
    auto t1 = timeNow(); 
    SetOneBranchAddress(h101, x);

    TH2d *hCorrRow, *hCorrCol;
    hCorrRow = new TH2D(Form("Corr ALPIDE%d vs. ALPIDE%d (ROWS)", x, y), Form("Corr ALPIDE%d vs. ALPIDE%d (ROWS)", x, y), 512,0,512,512,512,0,512);  
    hCorrCol = new TH2D(Form("Corr ALPIDE%d vs. ALPIDE%d (COLS)", x, y), Form("Corr ALPIDE%d vs. ALPIDE%d (COLS)", x, y), 1024,0,1024,1024,1024,0,1024);  

    firstEvent = std::min(firstEvent, (ulong)h101->GetEntries());
    ulong n = (maxEvents==0 || firstEvent+maxEvents > h101->GetEntries()) ? (h101->GetEntries()) : (firstEvent+maxEvents);
    maxEvents = n - firstEvent;
    ulong evCounter(0);

    printf("Entries in file: %d\n", h101->GetEntries()); 
    
    for(uint64_t evNum = firstEvent; evNum < n; ++evNum) {
        ++evCounter;
        if(evCounter%100 == 0) {
            cout << "\rEv#: " << evCounter << " Completed: " << Form("%.1f",(float)evCounter/maxEvents*100) <<"%"<< flush;
        }
        
        h101->GetEntry(evNum);
        
    }
    
    auto t2 = timeNow();
    cout << "\nTime taken: " << duration_cast<seconds>(t2-t1).count() << "s\n";
    app->Run();           
}

void RawCorrelations(const char* fileName, int det1, int det2, ulong firstEvent=0, ulong maxEvents=0) {
    assert(det1<=6 && det1>=1 && det2<=6 && det2>=1);

    TApplication* app = new TApplication("myApp", 0, 0);
    TFile* in = new TFile(fileName,"READ");
    if(!in || in->IsZombie()) {cerr << "Can't open rootfile with name: " << fileName << "\n"; exit(EXIT_FAILURE);}
    TTree* h101 = (TTree*)in->Get("h101");
    auto t1 = timeNow(); 
    
    SetOneBranchAddress(h101, x);
    SetOneBranchAddress(h101, y);
    
    UInt_t* mult1 = rowHandle[det1];
    UInt_t* mult2 = rowHandle[det2];
    UInt_t* RowV1   = RowHandle[det1];
    UInt_t* ColV1   = ColHandle[det1];
    UInt_t* RowV2   = RowHandle[det1];
    UInt_t* ColV2   = ColHandle[det1];

    TH2D* hCorrRow = new TH2D(Form("Corr ALPIDE%d vs. ALPIDE%d (ROWS)", det1, det2), Form("Corr ALPIDE%d vs. ALPIDE%d (ROWS)", det1, det2), 512,0,512,512,512,0,512);  
    TH2D* hCorrCol = new TH2D(Form("Corr ALPIDE%d vs. ALPIDE%d (COLS)", det1, det2), Form("Corr ALPIDE%d vs. ALPIDE%d (COLS)", det1, det2), 1024,0,1024,1024,1024,0,1024);  
    TCanvas* canvasRow = new TCanvas(Form("ALPIDE%d vs. ALPIDE%d (ROWS)", det1, det2), Form("ALPIDE%d vs. ALPIDE%d (ROWS)"), 2100,700);
    TCanvas* canvasCol = new TCanvas(Form("ALPIDE%d vs. ALPIDE%d (COLS)", det1, det2), Form("ALPIDE%d vs. ALPIDE%d (COLS)"), 2100,700);

    firstEvent = std::min(firstEvent, (ulong)h101->GetEntries());
    uint64_t n = (maxEvents==0 || firstEvent+maxEvents > h101->GetEntries()) ? (h101->GetEntries()) : (firstEvent+maxEvents);
    maxEvents = n - firstEvent;
    uint64_t evCounter(0);

    printf("Entries in file: %d\n", h101->GetEntries()); 
    for(uint64_t evNum = firstEvent; evNum < n; ++evNum) {
        ++evCounter;
        if(evCounter%100 == 0) {
            cout << "\rEv#: " << evCounter << " Completed: " << Form("%.1f",(float)evCounter/maxEvents*100) <<"%"<< flush;
            //printProgress((evNum+1.)/n);
        }
        
        h101->GetEntry(evNum);
        for(uint i1=0; i1<mult1; ++i1) {
            for(uint i2=0; i2<mult2; ++i2) {
                hCorrCol->Fill(ColV1[i1], ColV2[i2]);
                hCorrRow->Fill(RowV1[i1], RowV2[i2]);
            }
        } 
    }
    
    auto t2 = timeNow();
    cout << "\nTime taken: " << duration_cast<seconds>(t2-t1).count() << "s\n";
    app->Run();           

}
void AnalyseWithClustering(const char* fileName, ulong firstEvent=0, ulong maxEvents=0, int veto=0) { 
    TApplication* app = new TApplication("myApp", 0, 0);
    TFile* in = new TFile(fileName,"READ");
    if(!in || in->IsZombie()) {cerr << "Can't open rootfile with name: " << fileName << "\n"; exit(EXIT_FAILURE);}
    TTree* h101 = (TTree*)in->Get("h101");
    auto t1 = timeNow();

    /* WR Container 
    UInt_t hi1,lo1;
    UInt_t hi2,lo2;
    UInt_t hi3,lo3;
    UInt_t hi4,lo4;
    UInt_t hi5,lo5;
    UInt_t hi6,lo6; */
    /*
    h101->SetBranchAddress("ALPIDE1ROWv",Row1);
    h101->SetBranchAddress("ALPIDE1COLv",Col1);
    h101->SetBranchAddress("ALPIDE1ROW",&row1);
    h101->SetBranchAddress("ALPIDE1COL",&col1);
    
    h101->SetBranchAddress("ALPIDE2ROWv",Row2);
    h101->SetBranchAddress("ALPIDE2COLv",Col2);
    h101->SetBranchAddress("ALPIDE2ROW",&row2);
    h101->SetBranchAddress("ALPIDE2COL",&col2);

    h101->SetBranchAddress("ALPIDE3ROWv",Row3);
    h101->SetBranchAddress("ALPIDE3COLv",Col3);
    h101->SetBranchAddress("ALPIDE3ROW",&row3);
    h101->SetBranchAddress("ALPIDE3COL",&col3);
    
    h101->SetBranchAddress("ALPIDE4ROWv",Row4);
    h101->SetBranchAddress("ALPIDE4COLv",Col4);
    h101->SetBranchAddress("ALPIDE4ROW",&row4);
    h101->SetBranchAddress("ALPIDE4COL",&col4);
    
    h101->SetBranchAddress("ALPIDE5ROWv",Row5);
    h101->SetBranchAddress("ALPIDE5COLv",Col5);
    h101->SetBranchAddress("ALPIDE5ROW",&row5);
    h101->SetBranchAddress("ALPIDE5COL",&col5);
    
    h101->SetBranchAddress("ALPIDE6ROWv",Row6);
    h101->SetBranchAddress("ALPIDE6COLv",Col6);
    h101->SetBranchAddress("ALPIDE6ROW",&row6);
    h101->SetBranchAddress("ALPIDE6COL",&col6);
    */
    SetAllBranchAddress(h101);
    
    TH2D *hMeanCluster[7];
    TH2D *hSigmaCluster[7];
    TH1I *hMultpCluster[7];
    TCanvas *canvasCluster[7];
    
    for(int i=1; i<7; ++i) {
        hMeanCluster[i]  = new TH2D(Form("Cluster Mean #%d", i), Form("Cluster Mean #%d", i), 1024,0,1024,512,0,512);
        hSigmaCluster[i] = new TH2D(Form("Cluster Sigma #%d", i), Form("Cluster Sigma #%d", i), 20,0,5,20,0,5); 
        hMultpCluster[i] = new TH1I(Form("Cluster Multp #%d", i), Form("Cluster Multp #%d", i), 10,0,10);
        canvasCluster[i] = new TCanvas(Form("Cl-Analysis #%d", i),Form("Cl-Analysis #%d", i), 2100,700);
        canvasCluster[i]->Divide(3,1);
    }
    
    #if 0
    linexz = new TF1("linexz","[0]*x + [1]",0,15);
    lineyz = new TF1("lineyz","[0]*x + [1]",0,15);
    for(int i = 1 ; i < 7 ; i ++){
        HITMAP[i] = new TH2D(Form("HITMAP%d",i),Form("HITMAP%d",i),1024,0,1024,512,0,512);
    }
    ANG = new TH2D("ANG","Angle",1000,-0.05,0.05,1000,-0.05,0.05);
    TFile*ofile = new TFile("out.root","RECREATE");
    #endif
    
    // Container for clusters //
    vector<vector<Point>> clusters[7];

    firstEvent = std::min(firstEvent, (ulong)h101->GetEntries());
    uint64_t n = (maxEvents==0 || firstEvent+maxEvents > h101->GetEntries()) ? (h101->GetEntries()) : (firstEvent+maxEvents);
    maxEvents = n - firstEvent;
    uint64_t evCounter(0);

    printf("Entries in file: %d\n", h101->GetEntries()); 
    for(uint64_t evNum = firstEvent; evNum < n; ++evNum) {
        ++evCounter;
        if(evCounter%100 == 0) {
            cout << "\rEv#: " << evCounter << " Completed: " << Form("%.1f",(float)evCounter/maxEvents*100) <<"%"<< flush;
            //printProgress((evNum+1.)/n);
        }
        
        h101->GetEntry(evNum);
        /* ~~~ MARTIN PART ~~~ */
        // PART 1 : Construct clusters
        if(row1 > 1 && row1 == col1) {
            clusters[1] = ConstructClusters(Col1, Row1, row1);
        } else clusters[1].clear();

        if(row2 > 1 && row2 == col2) {
            clusters[2] = ConstructClusters(Col2, Row2, row2);
        } else clusters[2].clear();

        if(row3 > 1 && row3 == col3) {
            clusters[3] = ConstructClusters(Col3, Row3, row3);
        } else clusters[3].clear();

        if(row4 > 1 && row4 == col4) {
            clusters[4] = ConstructClusters(Col4, Row4, row4);
        } else clusters[4].clear();

        if(row5 > 1 && row5 == col5) {
            clusters[5] = ConstructClusters(Col5, Row5, row5);
        } else clusters[5].clear();

        if(row6 > 1 && row6 == col6) {
            clusters[6] = ConstructClusters(Col6, Row6, row6);
        } else clusters[6].clear();

        // PART 2 : Go over constructed clusters and feed into histograms

        for(int i=1; i<7; ++i) {
            if(clusters[i].size() == 0) continue;
            for(auto& cluster : clusters[i]) {
                auto fit = FitCluster(cluster);
                int multpCl = std::get<4>(fit);
                if(multpCl > veto) {
                    hMeanCluster[i]->Fill(std::get<0>(fit), std::get<1>(fit));
                    hSigmaCluster[i]->Fill(std::get<2>(fit), std::get<3>(fit));
                    hMultpCluster[i]->Fill(std::get<4>(fit));
                }
            }
        }
    }

    for(int i=1; i<7; ++i) {
        canvasCluster[i]->cd(1);
        hMeanCluster[i]->Draw("colz");
        canvasCluster[i]->cd(2);
        hSigmaCluster[i]->Draw("colz");
        canvasCluster[i]->cd(3);
        hMultpCluster[i]->Draw();
    }
    auto t2 = timeNow();
    cout << "\nTime taken: " << duration_cast<seconds>(t2-t1).count() << "s\n";
    app->Run();           
}

int main(int argc, char* argv[]) {
    if(IsHelpArg(argc, argv)) {cout << analyse_help; return 0;}
    
    string pStr;
    string fileName("");
    ulong firstEvent(0);
    ulong maxEvents(0);
    int veto(0);
    
    ParseCmdLine("file", fileName, argc, argv);
    
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
    
    AnalyseWithClustering(fileName.c_str(), firstEvent, maxEvents, veto);
    cout<<endl; return 0;
}

/* AUX FUNCTIONS */
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

