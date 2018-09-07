#include "CalcDistanceMatrixMain.h"

int main(int argc, char *argv[]) {
	int vecDim			= 80;
	int numRowsA		= 6000;
	int numRowsB		= 4000;
	int numInIter		= 100;
	int unitTestFlag;
	
	if (argc == 1) {
		unitTestFlag = UNIT_TEST_CALC_DISTANCE_MATRIX_VANILLA;
	}
	else if (argc > 2) {
		printf("Error: Number of Argument Is Larger Than 1");
		return 0;
	}
	else { // argc ==2 -> Single Input
		unitTestFlag = atoi(argv[1]);
	}

	switch (unitTestFlag)
	{
	case UNIT_TEST_CALC_DISTANCE_MATRIX_VANILLA:
		CalcDistanceMatrixVanillaUnitTest(vecDim, numRowsA, numRowsB, numInIter);
		break;
	case UNIT_TEST_CALC_DISTANCE_MATRIX_SSE:
		CalcDistanceMatrixSseUnitTest(vecDim, numRowsA, numRowsB, numInIter);
		break;
	case UNIT_TEST_CALC_DISTANCE_MATRIX_AVX:
		CalcDistanceMatrixAvxUnitTest(vecDim, numRowsA, numRowsB, numInIter);
		break;
	case UNIT_TEST_CALC_DISTANCE_MATRIX_EIGEN:
		CalcDistanceMatrixEigenUnitTest(vecDim, numRowsA, numRowsB, numInIter);
		break;
	case UNIT_TEST_CALC_DISTANCE_MATRIX_EIGEN_M:
		CalcDistanceMatrixEigenMUnitTest(vecDim, numRowsA, numRowsB, numInIter);
		break;
	default:
		printf("Error: Invalid Test Number");
		break;
	}


	return 0;
}

