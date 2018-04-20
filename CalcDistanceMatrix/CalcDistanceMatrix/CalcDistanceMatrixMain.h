#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <stdint.h>
#include "CalcDistanceMatrix.h"

#define UNIT_TEST_CALC_DISTANCE_MATRIX_VANILLA		1
#define UNIT_TEST_CALC_DISTANCE_MATRIX_SSE			2
#define UNIT_TEST_CALC_DISTANCE_MATRIX_AVX			3
#define UNIT_TEST_CALC_DISTANCE_MATRIX_EIGEN		4
#define UNIT_TEST_CALC_DISTANCE_MATRIX_EIGEN_M		5


// FA Library Functions
void CalcDistanceMatrixVanilla(float* mD, float* mA, float* mB, int vecDim, int numRowsA, int numRowsB);
void CalcDistanceMatrixSse(float* mD, float* mA, float* mB, int vecDim, int numRowsA, int numRowsB);
void CalcDistanceMatrixAvx(float* mD, float* mA, float* mB, int vecDim, int numRowsA, int numRowsB);
void CalcDistanceMatrixEigen(float* mD, float* mA, float* mB, int vecDim, int numRowsA, int numRowsB);
void CalcDistanceMatrixEigenM(float* mD, float* mA, float* mB, int vecDim, int numRowsA, int numRowsB);
void CalcDistanceMatrixRefTime(float* mD, float* mA, float* mB, int vecDim, int numRowsA, int numRowsB);

// Unit Test Functions 
void CalcDistanceMatrixVanillaUnitTest(int vecDim, int numRowsA, int numRowsB, int numIter);
void CalcDistanceMatrixSseUnitTest(int vecDim, int numRowsA, int numRowsB, int numIter);
void CalcDistanceMatrixAvxUnitTest(int vecDim, int numRowsA, int numRowsB, int numIter);
void CalcDistanceMatrixEigenUnitTest(int vecDim, int numRowsA, int numRowsB, int numIter);
void CalcDistanceMatrixEigenMUnitTest(int vecDim, int numRowsA, int numRowsB, int numIter);

// Auxiliary Functions
void PrintRunTimeData(double *vRunTime, int numIter);
int QSortCompFunDouble(const void * a, const void * b);

