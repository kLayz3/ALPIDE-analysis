#ifndef ALPIDE_TRACK_HH
#define ALPIDE_TRACK_HH

#include "libs.hh"
#include "spec.hh"

typedef uint32_t uint;

struct AlpideTrack {
	std::vector<int> tDet;
	std::vector<double> tZ;  // Z positions of each hit, in mm
	std::vector<double> tX;  // X positions of each hit, in mm
	std::vector<double> tY;  // Y positions of each hit, in mm
	std::vector<uint> tSize; // cluster size of each hit in the track

	AlpideTrack() = default;
	AlpideTrack(std::vector<int>& detV, std::vector<float>& colV, std::vector<float>& rowV, std::vector<uint>& sizeV);
	
	void GetFit(double& x0, double& y0, double& phiX, double& phiY);

	int GetSize() const {return tX.size();}
	bool IsValidTrack() const {return (tX.size()>2) ? 1 : 0;}
	void PrintTrack() const {}

	std::vector<double> XToCol();
	std::vector<double> YToRow();

	/* static fields */
	// z-coordinates given in spec.hh: Z_Alpide[] ... however, useful to have minZ, maxZ, meanZ;
	static const double minZ; 
	static const double maxZ;
	static const double meanZ;
	inline static std::array<double,ALPIDE_NUM+1> aCol;    // X slope  of each individual detector (relative to first)
	inline static std::array<double,ALPIDE_NUM+1> bCol;    // X offset of each individual detector (relative to first), in pixels
	inline static std::array<double,ALPIDE_NUM+1> aRow;    // Y slope  of each individual detector (relative to first)
	inline static std::array<double,ALPIDE_NUM+1> bRow;    // Y offset of each individual detector (relative to first), in pixels
	
	static float tolCol;				/* How many columns (in pixels) difference (x-direction) between neighbouring 
										 * detectors can be tolerated and still call the hits part of a track */
	static float tolRow;				// -||- but for rows (y-direction);
	
	/* static methods to be used in main functions */
	static void LoadCal(const char* fileName);
	static void ColVCal(std::vector<int>& detV, std::vector<float>& colV);
	static void RowVCal(std::vector<int>& detV, std::vector<float>& rowV);

	static double ColToX(float col) {return (double)col*PIXEL_SIZE_X + PIXEL_SIZE_X2;}
	static double RowToY(float row) {return (double)row*PIXEL_SIZE_Y + PIXEL_SIZE_Y2;}
	static double DetToZ(int tDet) {return Z_Alpide[tDet];};

	static void SetColTolerance(float tol) {tolCol = tol;}
	static void SetRowTolerance(float tol) {tolRow = tol;}

protected:
	void PushData(int det, float col, float row, uint size);
	bool IsInTrack(int d0, float c0, float r0, int d1, int c1, int r1);

};
			

#endif 
