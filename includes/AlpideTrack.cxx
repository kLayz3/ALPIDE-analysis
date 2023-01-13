/* This class plays around the constructor which collects all suitable entries from (detV, colV, rowV, sizeV) event data. 
 * The passed-by-ref vectors have the collected elements removed. Call iteratively to construct 'all' tracks from 
 * an event data. A track is constructed such that all collected hits make a corresponding line with respect to the 
 * first entry (detV[0], colV[0], rowV[0], sizeV[0]).
 * The (colV,rowV) vectors HAVE to be calibrated first, by passing them first to AlpideTrack::RowVCal()
 * and AlpideTrack::ColVCal() methods. 
 * Calibration is loaded by passing a file name to AlpideTrack::LoadCal() call.
 * -- Written by Martin Bajzek (M.Bajzek@gsi.de) */

#include "AlpideTrack.h"
#include "AuxFunctions.h"

#define TOLERANCE 30

using namespace std;

const double AlpideTrack::minZ = *std::min_element(Z_Alpide+1, Z_Alpide+ALPIDE_NUM+1);
const double AlpideTrack::maxZ = *std::max_element(Z_Alpide+1, Z_Alpide+ALPIDE_NUM+1);
const double AlpideTrack::meanZ = std::accumulate(Z_Alpide+1, Z_Alpide+ALPIDE_NUM+1, 0.)/ALPIDE_NUM;
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
			PushData(detV[i], colV[i], rowV[i], sizeV[i]); // pushes to the vector fields
			
			isCollected[d] = 1; // flags this detector as already collected. It can't host another hit.
			QuickErase(detV,i); QuickErase(colV,i); QuickErase(rowV,i); QuickErase(sizeV,i);
			--i;
		}
	}
}

void AlpideTrack::GetFit(double& x0, double& y0, double& phiX, double& phiY) {
	if((this->tDet).size() < 3) {cerr<<__PRETTY_FUNCTION__<< " , vector size is < 3.\n"; return;}
	
	TF1 lineX("lineX","[0]*x+[1]", minZ, maxZ);
	TF1 lineY("lineY","[0]*x+[1]", minZ, maxZ);
	TGraph gX; TGraphFromVector(gX, tZ, tX);
	TGraph gY; TGraphFromVector(gY, tZ, tY);
	gX.Fit(&lineX, "Q");
	gY.Fit(&lineY, "Q");
	double aX = lineX.GetParameter(0); double bX = lineX.GetParameter(1);
	double aY = lineY.GetParameter(0); double bY = lineY.GetParameter(1);
	x0 = bX + aX*Z_Alpide[1]; // x position on ALPIDE1
	y0 = bY + aY*Z_Alpide[1]; // y position on ALPIDE1
	phiX = aX;
	phiY = aY;
}
vector<double> AlpideTrack::XToCol() {
	vector<double> vec(tX);
	for_each(vec.begin(),  vec.end(), [](auto& x){
			x -= PIXEL_SIZE_X2;
			x /= PIXEL_SIZE_X;});
	return vec;
}
vector<double> AlpideTrack::YToRow() {
	vector<double> vec(tY);
	for_each(vec.begin(),  vec.end(), [](auto& y){
			y -= PIXEL_SIZE_Y2;
			y /= PIXEL_SIZE_Y;});
	return vec;
}

void AlpideTrack::LoadCal(const char* fileName) {
	TFile* in = new TFile(fileName,"READ");
	if(!in || in->IsZombie()) {
		cerr << "Can't open rootfile: " << fileName << __PRETTY_FUNCTION__ <<"\n";
		aCol.fill(1.); aRow.fill(1.);
		bCol.fill(0.); bRow.fill(0.);
		return;
	}
	TVectorD* slopeCol(in->Get<TVectorD>("aCol"));
	TVectorD* offsetCol(in->Get<TVectorD>("bCol"));
	TVectorD* slopeRow(in->Get<TVectorD>("aRow"));
	TVectorD* offsetRow(in->Get<TVectorD>("bRow"));
	for(int i=0; i<=ALPIDE_NUM; ++i) {
		aCol[i] = slopeCol->operator[](i);
		bCol[i] = offsetCol->operator[](i);
		aRow[i] = slopeRow->operator[](i);
		bRow[i] = offsetRow->operator[](i);
	}
	/* memcpy(aCol.data(), slopeCol->fArray, sizeof(double)*(ALPIDE_NUM+1)); */
	/* memcpy(bCol.data(), offsetCol->fArray, sizeof(double)*(ALPIDE_NUM+1)); */
	/* memcpy(aRow.data(), slopeRow->fArray, sizeof(double)*(ALPIDE_NUM+1)); */
	/* memcpy(bRow.data(), offsetRow->fArray, sizeof(double)*(ALPIDE_NUM+1)); */

	in->Close(); delete in;
}

/* Calibrates colV and rowV vectors */
void AlpideTrack::ColVCal(vector<int>& detV, vector<float>& colV) {
	for(int i=0; i<detV.size(); ++i) {
		auto d = detV[i]; // 1,2,3, ... ALPIDE_NUM
		colV[i] = aCol[d] * colV[i] + bCol[d];
	}
}
void AlpideTrack::RowVCal(vector<int>& detV, vector<float>& rowV) {
	for(int i=0; i<detV.size(); ++i) {
		auto d = detV[i];
		rowV[i] = aRow[d] * rowV[i] + bRow[d];
	}
}

void AlpideTrack::PushData(int det, float col, float row, uint size) {
	(this->tDet).push_back(det);
	tZ.push_back(DetToZ(det));
	tX.push_back(ColToX(col));
	tY.push_back(RowToY(row));
	tSize.push_back(size);
}

bool AlpideTrack::IsInTrack(int d0, float c0, float r0, int d1, int c1, int r1) {
	int d = abs(d1 - d0);
	if((abs(c1-c0) < d*tolCol) && (abs(r1-r0) < d*tolRow)) return true;
	else return false;
}

