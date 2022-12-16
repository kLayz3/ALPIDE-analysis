#ifndef ALPIDE_TRACK
#define ALPIDE_TRACK

#include "libs.hh"
#include "spec.hh"

using std::abs;
typedef uint32_t uint;

struct AlpideTrack {
	std::vector<int> det;
	std::vector<double> tX;
	std::vector<double> tY;
	std::vector<uint> tSize;
	
	double fitX;
	double fitY;

	AlpideTrack() = default;
	AlpideTrack(std::vector<int>& detV, std::vector<float>& colV, std::vector<float>& rowV, std::vector<uint>& sizeV);
	
	void GetFit(double& x0, double& y0, double& theta, double& phi);
	virtual double X0()=0;
	virtual double Y0()=0;
	virtual double Theta()=0;
	virtual double Phi()=0;
	
	int GetSize();
	bool IsValidTrack();

	static double ds;                 // separation between two detectors in millimeters
	inline static std::vector<double> aCol;  // X slope  of each individual detector (relative to first)
	inline static std::vector<double> bCol;  // X offset of each individual detector (relative to first), in pixels
	inline static std::vector<double> aRow;  // Y slope  of each individual detector (relative to first)
	inline static std::vector<double> bRow;  // Y offset of each individual detector (relative to first), in pixels
	
	static float tolCol;				/* How many columns (in pixels) difference (x-direction) between neighbouring 
										 * detectors can be tolerated o still call the hits part of a track */
	static float tolRow;				// -||- for rows (y-direction);
		
	static void LoadCalibration(const char* fileName);
	static void ColVCal(std::vector<uint>& detV, std::vector<float>& colV);
	static void RowVCal(std::vector<uint>& detV, std::vector<float>& rowV);
	static double ColToX(float col);
	static double RowToY(float row);
	static void SetDetectorSeparation(double dx) {ds = dx;}
	static void SetColTolerance(float tol) {tolCol = tol;}
	static void SetRowTolerance(float tol) {tolRow = tol;}

protected:
	void PushData(int det, float col, float row, uint size);
	bool IsInTrack(int d0, float c0, float r0, int d1, int c1, int r1);
};
			

#endif 
