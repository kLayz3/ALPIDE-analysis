#ifndef SPEC_HH
#define SPEC_HH

#include "cmath"
#include "algorithm"

/* Chip specific defines */
/* Pixel size in mm */
//#define PIXEL_SIZE_X 0.02924 // column  
#define PIXEL_SIZE_X 0.02688 // column  --found lrose
#define PIXEL_SIZE_Y 0.02924 // row 
#define PIXEL_SIZE_X2 PIXEL_SIZE_X/2
#define PIXEL_SIZE_Y2 PIXEL_SIZE_Y/2

#define MAX_HITS 2048 /* max hits in raw root file per an event */
#define MAX_CLUSTERS 4096 /* max clusters per an event */

/* deprecated */
#define ALPIDE_NUM 6 /* number of sensors in the telescope */

#define MAX_BOARD_ID 4
#define MAX_CHIP_ID 6
/* X,Y,Z-coordinates of each ALPIDE chip detector in the setup. in mm */

inline double alpide_def_pos[MAX_BOARDS][MAX_CHIP_ID+1][3];
//           boardId -----^            ^         ^---- X=0,Y=1,Z=2
//                         chipId -----^

inline void alpide_default_positions() { 
	std::fill((double*)alpide_def_pos, (double*)alpide_def_pos + (MAX_BOARD_ID)*(MAX_CHIP_ID+1)*3, NAN);
	alpide_def_pos[1][1][0] = 0;
	alpide_def_pos[1][1][1] = 0;
	alpide_def_pos[1][1][2] = 305;

	alpide_def_pos[1][2][0] = 0;
	alpide_def_pos[1][2][1] = 0;
	alpide_def_pos[1][2][2] = 350;

	alpide_def_pos[1][5][0] = 0;
	alpide_def_pos[1][5][1] = 0;
	alpide_def_pos[1][5][2] = 390;

	alpide_def_pos[2][6][0] = 0;
	alpide_def_pos[2][6][1] = 0;
	alpide_def_pos[2][6][2] = 1097;

	alpide_def_pos[3][4][0] = 0;
	alpide_def_pos[3][4][1] = 0;
	alpide_def_pos[3][4][2] = 1142;

	alpide_def_pos[4][3][0] = 0;
	alpide_def_pos[4][3][1] = 0;
	alpide_def_pos[4][3][2] = 1182;
}


#endif 

