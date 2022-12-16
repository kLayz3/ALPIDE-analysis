/* This class plays around the constructor which collects all suitable entries from (detV, colV, rowV, sizeV) event data. 
 * The passed-by-ref vectors have the collected elements removed. Call iteratively to construct 'all' tracks from 
 * an event data. A track is constructed such that all collected hits make a corresponding line with respect to the 
 * first entry (detV[0], ... , sizeV[0]).
 * The (colV,rowV) vectors have to be calibrated first, by passing them first to AlpideTrack::RowVCal(..) 
 * and AlpideTrack::ColVCal() methods. Calibration is loaded by passing a file name to AlpideTrack::LoadCalibration() call.
 * -- Written by Martin Bajzek (M.Bajzek@gsi.de) */

#include "AlpideTrack.h"
#include "AuxFunctions.h"

#define PIXEL_SIZE_X 0.02924 // in mm
#define PIXEL_SIZE_Y 0.02688 
#define TOLERANCE 50

using namespace std;

double AlpideTrack::ds = DELTA_S; // from spec.hh
float AlpideTrack::tolCol = TOLERANCE;
float AlpideTrack::tolRow = TOLERANCE;

AlpideTrack::AlpideTrack(vector<int>& detV, vector<float>& colV, vector<float>& rowV, vector<uint>& sizeV) {
	bool isCollected[ALPIDE_NUM+1] = {0}; // indicator whether this detector is already included in the track
	
	auto d0 = detV[0]; auto c0 = colV[0]; auto r0 = rowV[0]; auto s0 = sizeV[0];
	isCollected[d0] = 1;
	PushData(d0, c0, r0, s0);
	QuickErase(detV,0); QuickErase(colV,0); QuickErase(rowV,0); QuickErase(sizeV,0);
		
	for(int i=0; i<detV.size(); ++i) {
		int d = detV[i];
		if(isCollected[d]) continue;
		if(IsInTrack(d0,c0,r0, detV[i],colV[i],rowV[i])) {
			PushData(detV[i], colV[i], rowV[i], sizeV[i]); // pushes to the fields
			
			isCollected[d] = 1;
			QuickErase(detV,i); QuickErase(colV,i); QuickErase(rowV,i); QuickErase(sizeV,i);
			--i;
		}
	}
}
void AlpideTrack::GetFit(double& x0, double& y0, double& theta, double& phi) {
	if(det.size()<2) {cerr<<__PRETTY_FUNCTION__<< " , vector size is < 2."; return;}
	TF1 lineX("lineX","[0]*x+[1]", det.front(), det.back());
	TF1 lineY("lineY","[0]*x+[1]", det.front(), det.back());
	TGraph gX; TGraphFromVector(gX, det, tX);
	TGraph gY; TGraphFromVector(gX, det, tY);
	gX.Fit(&lineX);
	gY.Fit(&lineY); 
	double aX = lineX.GetParameter(0), bX = lineX.GetParameter(1);
	double aY = lineY.GetParameter(0), bY = lineY.GetParameter(1);
	x0 = aX + bX;
	y0 = aY + bY;
	theta = sqrt(aX*aX + aY*aY);
	phi = atan(aY/aX);
}

int AlpideTrack::GetSize() {return tX.size();}
bool AlpideTrack::IsValidTrack() {return (tX.size()>0) ? 1 : 0;}

void AlpideTrack::LoadCalibration(const char* fileName) {
	TFile* in = new TFile(fileName,"READ");
	if(!in || in->IsZombie()) {
		cerr << "Can't open rootfile: " << fileName << __PRETTY_FUNCTION__ <<"\n";
		exit(EXIT_FAILURE);
	}
	TVectorD* slopeCol{in->Get<TVectorD>("aCol")};
	TVectorD* offsetCol{in->Get<TVectorD>("bCol")};
	TVectorD* slopeRow{in->Get<TVectorD>("aRow")};
	TVectorD* offsetRow{in->Get<TVectorD>("bRow")};
	
	aCol.clear(); bCol.clear(); aRow.clear(); bRow.clear();
	for(int i=0; i<=ALPIDE_NUM; ++i) {
		aCol.push_back(slopeCol->operator[](i)); 	
		bCol.push_back(offsetCol->operator[](i)); 	
		aRow.push_back(slopeRow->operator[](i)); 	
		bRow.push_back(offsetRow->operator[](i)); 	
	}
	in->Close();
}

void AlpideTrack::ColVCal(vector<uint>& detV, vector<float>& colV) {
	for(int i=0; i<detV.size(); ++i) {
		auto x = detV[i];
		colV[i] = aCol[x]*colV[i] + bCol[x];
	}
}
void AlpideTrack::RowVCal(vector<uint>& detV, vector<float>& rowV) {
	for(int i=0; i<detV.size(); ++i) {
		auto x = detV[i];
		rowV[i] = aRow[x]*rowV[i] + bRow[x];
	}
}
double AlpideTrack::ColToX(float col) {
	return (double)col * PIXEL_SIZE_X;
}
double AlpideTrack::RowToY(float row) {
	return  (double)row * PIXEL_SIZE_Y;
}

void AlpideTrack::PushData(int det, float col, float row, uint size) {
	(this->det).push_back(det);
	tX.push_back(ColToX(col));
	tY.push_back(RowToY(row));
	tSize.push_back(size);
}

bool AlpideTrack::IsInTrack(int d0, float c0, float r0, int d1, int c1, int r1) {
	int d = d1 - d0;
	if((abs(c1-c0) < d*tolCol) && (abs(r1-r0) < d*tolRow)) return true;
	else return false;
}
