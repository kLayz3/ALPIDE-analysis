#include <bits/stdc++.h>
using namespace std;

#define NO_CHIP {NAN,NAN,NAN}

int main() {
	/* X,Y,Z-coordinates of each ALPIDE chip detector in the setup. in mm */
	double XYZ_Alpide[10][7][3];
	std::fill((double*)XYZ_Alpide, (double*)XYZ_Alpide + 10*7*3, NAN);

	XYZ_Alpide[1][1][0] = 0;
	XYZ_Alpide[1][1][1] = 0;
	XYZ_Alpide[1][1][2] = 0;

	XYZ_Alpide[1][2][0] = 0;
	XYZ_Alpide[1][2][1] = 0;
	XYZ_Alpide[1][2][2] = 200;

	XYZ_Alpide[1][5][0] = 0;
	XYZ_Alpide[1][5][1] = 0;
	XYZ_Alpide[1][5][2] = 400;

	XYZ_Alpide[2][3][0] = 0;
	XYZ_Alpide[2][3][1] = 0;
	XYZ_Alpide[2][3][2] = 600;

	for(int i=0;i<10;++i) {
		printf("\n###### BOARD: %d ###########\n", i);
		for(int j=0;j<7;++j) {
			printf("ChipID: %d -> (", j);
			for(int k=0;k<3;++k) {
				printf("%.1f ", XYZ_Alpide[i][j][k]);
			}
			printf(")\n");
		}
	}
}
