
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <stdint.h>
#include "SleefVsSvml.h"

#define UNIT_TEST_SINE						1
#define UNIT_TEST_COSINE					2
#define UNIT_TEST_EXP						3

#define NUM_FUN_FLAVOR 4

// Main Function
void SineSleefSse(float* vO, float* vI, int numElements);
void SineSleefAvx(float* vO, float* vI, int numElements);
void SineSvmlSse(float* vO, float* vI, int numElements);
void SineSvmlAvx(float* vO, float* vI, int numElements);
//void CosineSleefSse(float* vO, float* vI, int numElements);
//void CosineSleefAvx(float* vO, float* vI, int numElements);
//void CosineSvmlSse(float* vO, float* vI, int numElements);
//void CosineSvmlAvx(float* vO, float* vI, int numElements);
void ExpSleefSse(float* vO, float* vI, int numElements);
void ExpSleefAvx(float* vO, float* vI, int numElements);
void ExpSvmlSse(float* vO, float* vI, int numElements);
void ExpSvmlAvx(float* vO, float* vI, int numElements);

// Unit Test Functions 
void UnitTestSine(int numElements, int numIter, float maxVal, int negValFlag);
//void UnitTestCosine(int numElements, int numIter, float maxVal, int negativeValSupport);
//void UnitTestExp(int numElements, int numIter, float maxVal, int negativeValSupport);


// Auxiliary Functions
void PrintRunTimeData(double *vRunTime, int numIter);
int QSortCompFunDouble(const void * a, const void * b);

