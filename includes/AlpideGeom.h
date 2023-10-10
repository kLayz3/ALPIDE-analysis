#ifndef ALPIDE_GEOM_H
#define ALPIDE_GEOM_H

#include "spec.hh"
#include <cstdint>

/* This is a class wrapping the container 
 * of positions of individual chips. */

class AlpideGeom {
public:
	double pos[MAX_BOARD_ID+1][MAX_CHIP_ID+1][3];
//  boardId --------^            ^        ^---- X=0,Y=1,Z=2
//                chipId --------^
	double* _end_ptr;

	AlpideGeom();
	void SetDefaultPositions();
	
	/* Non-trivial methods to extract slices,
	 * where n=0,1,2 representing X,Y,Z respectively */
	template<uint32_t n>
	std::vector<double> Get(uint32_t boardId);

	template<uint32_t n>
	std::vector<double> Get();
};

#endif
