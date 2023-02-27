#ifndef ALPIDE_TRACK_H
#define ALPIDE_TRACK_H

#include "libs.hh"
#include "spec.hh"

struct AlpideTrack {
	std::vector<int> tDet;
	std::vector<double> tZ;  // Z positions of each hit, in mm
	std::vector<double> tX;  // X positions of each hit, in mm
	std::vector<double> tY;  // Y positions of each hit, in mm
	std::vector<uint> tSize; // cluster size of each hit in the track

	AlpideTrack() = default;
	AlpideTrack(std::vector<int>& detV, std::vector<double>& colV, std::vector<double>& rowV, std::vector<uint>& sizeV);
	
	void GetFit(double& x0, double& y0, double& phiX, double& phiY);

	int GetSize() const {return tX.size();}
	bool IsValidTrack() const {return (tX.size()>2) ? 1 : 0;}
	void PrintTrack() const {}

	std::vector<double> XToCol();
	std::vector<double> YToRow();

	// z-coordinates given in spec.hh: Z_Alpide[] ... however, useful to have minZ, maxZ, meanZ;
	static const double minZ; 
	static const double maxZ;
	static const double meanZ;

	// Calibration fields //
	inline static std::array<double,ALPIDE_NUM+1> aCol;    // X slope  of each individual detector (relative to first)
	inline static std::array<double,ALPIDE_NUM+1> bCol;    // X offset of each individual detector (relative to first), in pixels
	inline static std::array<double,ALPIDE_NUM+1> aRow;    // Y slope  of each individual detector (relative to first)
	inline static std::array<double,ALPIDE_NUM+1> bRow;    // Y offset of each individual detector (relative to first), in pixels
	
	static double tolCol; /* How many columns (in pixels) difference (x-direction) between neighbouring */ 
	static double tolRow; /* detectors can be tolerated and still call the hits part of a track, per meanZ */
						 /* -||- but for rows (y-direction); */
	
	/* static methods to be used in main functions */
	static void LoadCal(const char* fileName = nullptr);
	static void ColVCal(std::vector<int>& detV, std::vector<double>& colV);
	static void RowVCal(std::vector<int>& detV, std::vector<double>& rowV);

	static double ColToX(double col) {return (double)col*PIXEL_SIZE_X + PIXEL_SIZE_X2;}
	static double RowToY(double row) {return (double)row*PIXEL_SIZE_Y + PIXEL_SIZE_Y2;}
	static double DetToZ(int tDet) {return Z_Alpide[tDet];};

	static void SetColTolerance(double tol) {tolCol = tol;}
	static void SetRowTolerance(double tol) {tolRow = tol;}
	static void AngleXTolerance(double angle);
	static void AngleYTolerance(double angle);

protected:	
	void PushData(int det, double col, double row, uint size);
	bool IsInTrack(int d0, double c0, double r0, int d1, double c1, double r1);

};

#endif 
