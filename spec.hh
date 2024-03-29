#ifndef SPEC_HH
#define SPEC_HH

#include "cmath"

/* Chip specific defines */
/* Pixel size in mm */
//#define PIXEL_SIZE_X 0.02924 // column  
#define PIXEL_SIZE_X 0.02688 // column  --found lrose
#define PIXEL_SIZE_Y 0.02924 // row 

#define MAX_HITS 2048 /* max hits in raw root file per an event */
#define MAX_CLUSTERS 1024 /* max clusters per an event */
#define ALPIDE_NUM 6 /* number of sensors in the telescope */

/* Z-coordinates of each ALPIDE detector in the telescope */
const double Z_Alpide[] = {NAN,
	0,
	20,
	50,
	70,
	100,
	120
};

/* In case future experiments will start with ALPIDE0 as 'first' */
#define ALPIDE_FIRST_ID 1

#define PIXEL_SIZE_X2 PIXEL_SIZE_X/2
#define PIXEL_SIZE_Y2 PIXEL_SIZE_Y/2

#endif 

