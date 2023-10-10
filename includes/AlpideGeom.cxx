#include "AlpideGeom.h"
#include "cmath"
#include "algorithm"

AlpideGeom::AlpideGeom() {
	_end_ptr = (double*)pos + (MAX_BOARD_ID+1)*(MAX_CHIP_ID+1)*3;
	std::fill((double*)pos, _end_ptr, NAN);
}

// This will read the default positions from the spec.hh file,
// where users can write inputs to the array. This method needs to be called 
// only once
void AlpideGeom::SetDefaultPositions() {
	alpide_default_positions();
	memcpy(pos, alpide_def_pos, (MAX_BOARD_ID+1)*(MAX_CHIP_ID+1)*3*sizeof(double));
}


template<uint32_t n>
std::vector<double> AlpideGeom::Get(uint32_t boardId) {
	static_assert(n < 3, __PRETTY_FUNCTION__ " template parameter must be one of (0,1,2) = (X,Y,Z).\n");
	assert(boardId <= MAX_BOARD_ID && 
			__PRETTY_FUNCTION__ " , argument but be <= " std::to_string(MAX_BOARD_ID).c_str());
	std::vector<double> v;
	for(int chip=0; chip <= MAX_CHIP_ID; ++chip) {
		v.push_back(pos[boardId][chip][n]);
	}
	return v;
}

template<uint32_t n>
std::vector<double> AlpideGeom::Get() {
	static_assert(n < 3, __PRETTY_FUNCTION__ " template parameter must be one of (0,1,2) = (X,Y,Z).\n");
	std::vector<double> v;
	double* slider = (double*)pos + n;
	double* end = (double*)pos + (MAX_BOARD_ID+1)*(MAX_CHIP_ID+1)*3;
	while(slider < end) {
		if(isnormal(*slider)) v.push_back(*slider);
	}

	return v;
}
