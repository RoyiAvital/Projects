#include "CalcDistanceMatrixMain.h"

void CalcDistanceMatrixVanillaUnitTest(int vecDim, int numRowsA, int numRowsB, int numIter) {
	int ii;
	int inputSizeA, inputSizeB, inputSizeD;
	float *mA;
	float *mB;
	float *mD;
	clock_t clockStart, clockEnd;
	double *vRunTime;

	vRunTime = (double*)_mm_malloc(numIter * sizeof(double), SSE_ALIGNMENT);

	inputSizeA = numRowsA * vecDim;
	inputSizeB = numRowsB * vecDim;
	inputSizeD = numRowsA * numRowsB;

	mA = (float*)_mm_malloc(inputSizeA * sizeof(float), SSE_ALIGNMENT);
	mB = (float*)_mm_malloc(inputSizeB * sizeof(float), SSE_ALIGNMENT);
	mD = (float*)_mm_malloc(inputSizeD * sizeof(float), SSE_ALIGNMENT);

	for (ii = 0; ii < inputSizeA; ii++) {
		mA[ii] = (float)(rand()) / (float)(RAND_MAX);
	}

	for (ii = 0; ii < inputSizeB; ii++) {
		mB[ii] = (float)(rand()) / (float)(RAND_MAX);
	}

	clockStart = clock();

	for (ii = 0; ii < numIter; ii++)
	{
		clockStart = clock();

		CalcDistanceMatrixVanilla(mD, mA, mB, vecDim, numRowsA, numRowsB);

		clockEnd = clock();

		vRunTime[ii] = (double)(clockEnd - clockStart) / CLOCKS_PER_SEC;
	}
	
	printf("\n");
	printf("CalcDistanceMatrixVanilla Unit Test\n");
	PrintRunTimeData(vRunTime, numIter);

	_mm_free(vRunTime);
	_mm_free(mA);
	_mm_free(mB);
	_mm_free(mD);

}


void CalcDistanceMatrixSseUnitTest(int vecDim, int numRowsA, int numRowsB, int numIter) {
	int ii;
	int inputSizeA, inputSizeB, inputSizeD;
	float *mA;
	float *mB;
	float *mD;
	clock_t clockStart, clockEnd;
	double *vRunTime;

	vRunTime = (double*)_mm_malloc(numIter * sizeof(double), SSE_ALIGNMENT);

	inputSizeA = numRowsA * vecDim;
	inputSizeB = numRowsB * vecDim;
	inputSizeD = numRowsA * numRowsB;

	mA = (float*)_mm_malloc(inputSizeA * sizeof(float), SSE_ALIGNMENT);
	mB = (float*)_mm_malloc(inputSizeB * sizeof(float), SSE_ALIGNMENT);
	mD = (float*)_mm_malloc(inputSizeD * sizeof(float), SSE_ALIGNMENT);

	for (ii = 0; ii < inputSizeA; ii++) {
		mA[ii] = (float)(rand()) / (float)(RAND_MAX);
	}

	for (ii = 0; ii < inputSizeB; ii++) {
		mB[ii] = (float)(rand()) / (float)(RAND_MAX);
	}

	clockStart = clock();

	for (ii = 0; ii < numIter; ii++)
	{
		clockStart = clock();

		CalcDistanceMatrixSse(mD, mA, mB, vecDim, numRowsA, numRowsB);

		clockEnd = clock();

		vRunTime[ii] = (double)(clockEnd - clockStart) / CLOCKS_PER_SEC;
	}

	printf("\n");
	printf("CalcDistanceMatrixSse Unit Test\n");
	PrintRunTimeData(vRunTime, numIter);

	_mm_free(vRunTime);
	_mm_free(mA);
	_mm_free(mB);
	_mm_free(mD);

}


void CalcDistanceMatrixAvxUnitTest(int vecDim, int numRowsA, int numRowsB, int numIter) {
	int ii;
	int inputSizeA, inputSizeB, inputSizeD;
	float *mA;
	float *mB;
	float *mD;
	clock_t clockStart, clockEnd;
	double *vRunTime;

	vRunTime = (double*)_mm_malloc(numIter * sizeof(double), SSE_ALIGNMENT);

	inputSizeA = numRowsA * vecDim;
	inputSizeB = numRowsB * vecDim;
	inputSizeD = numRowsA * numRowsB;

	mA = (float*)_mm_malloc(inputSizeA * sizeof(float), AVX_ALIGNMENT);
	mB = (float*)_mm_malloc(inputSizeB * sizeof(float), AVX_ALIGNMENT);
	mD = (float*)_mm_malloc(inputSizeD * sizeof(float), AVX_ALIGNMENT);

	for (ii = 0; ii < inputSizeA; ii++) {
		mA[ii] = (float)(rand()) / (float)(RAND_MAX);
	}

	for (ii = 0; ii < inputSizeB; ii++) {
		mB[ii] = (float)(rand()) / (float)(RAND_MAX);
	}

	clockStart = clock();

	for (ii = 0; ii < numIter; ii++)
	{
		clockStart = clock();

		CalcDistanceMatrixAvx(mD, mA, mB, vecDim, numRowsA, numRowsB);

		clockEnd = clock();

		vRunTime[ii] = (double)(clockEnd - clockStart) / CLOCKS_PER_SEC;
	}

	printf("\n");
	printf("CalcDistanceMatrixAvx Unit Test\n");
	PrintRunTimeData(vRunTime, numIter);

	_mm_free(vRunTime);
	_mm_free(mA);
	_mm_free(mB);
	_mm_free(mD);

}


void CalcDistanceMatrixEigenUnitTest(int vecDim, int numRowsA, int numRowsB, int numIter) {
	int ii;
	int inputSizeA, inputSizeB, inputSizeD;
	float *mA;
	float *mB;
	float *mD;
	clock_t clockStart, clockEnd;
	double *vRunTime;

	vRunTime = (double*)_mm_malloc(numIter * sizeof(double), SSE_ALIGNMENT);

	inputSizeA = numRowsA * vecDim;
	inputSizeB = numRowsB * vecDim;
	inputSizeD = numRowsA * numRowsB;

	mA = (float*)_mm_malloc(inputSizeA * sizeof(float), SSE_ALIGNMENT);
	mB = (float*)_mm_malloc(inputSizeB * sizeof(float), SSE_ALIGNMENT);
	mD = (float*)_mm_malloc(inputSizeD * sizeof(float), SSE_ALIGNMENT);

	for (ii = 0; ii < inputSizeA; ii++) {
		mA[ii] = (float)(rand()) / (float)(RAND_MAX);
	}

	for (ii = 0; ii < inputSizeB; ii++) {
		mB[ii] = (float)(rand()) / (float)(RAND_MAX);
	}

	clockStart = clock();

	for (ii = 0; ii < numIter; ii++)
	{
		clockStart = clock();

		CalcDistanceMatrixEigen(mD, mA, mB, vecDim, numRowsA, numRowsB);

		clockEnd = clock();

		vRunTime[ii] = (double)(clockEnd - clockStart) / CLOCKS_PER_SEC;
	}

	printf("\n");
	printf("CalcDistanceMatrixEigen Unit Test\n");
	PrintRunTimeData(vRunTime, numIter);

	_mm_free(vRunTime);
	_mm_free(mA);
	_mm_free(mB);
	_mm_free(mD);

}


void CalcDistanceMatrixEigenMUnitTest(int vecDim, int numRowsA, int numRowsB, int numIter) {
	int ii;
	int inputSizeA, inputSizeB, inputSizeD;
	float *mA;
	float *mB;
	float *mD;
	clock_t clockStart, clockEnd;
	double *vRunTime;

	vRunTime = (double*)_mm_malloc(numIter * sizeof(double), SSE_ALIGNMENT);

	inputSizeA = numRowsA * vecDim;
	inputSizeB = numRowsB * vecDim;
	inputSizeD = numRowsA * numRowsB;

	mA = (float*)_mm_malloc(inputSizeA * sizeof(float), SSE_ALIGNMENT);
	mB = (float*)_mm_malloc(inputSizeB * sizeof(float), SSE_ALIGNMENT);
	mD = (float*)_mm_malloc(inputSizeD * sizeof(float), SSE_ALIGNMENT);

	for (ii = 0; ii < inputSizeA; ii++) {
		mA[ii] = (float)(rand()) / (float)(RAND_MAX);
	}

	for (ii = 0; ii < inputSizeB; ii++) {
		mB[ii] = (float)(rand()) / (float)(RAND_MAX);
	}

	clockStart = clock();

	for (ii = 0; ii < numIter; ii++)
	{
		clockStart = clock();

		CalcDistanceMatrixEigenM(mD, mA, mB, vecDim, numRowsA, numRowsB);

		clockEnd = clock();

		vRunTime[ii] = (double)(clockEnd - clockStart) / CLOCKS_PER_SEC;
	}

	printf("\n");
	printf("CalcDistanceMatrixEigenM Unit Test\n");
	PrintRunTimeData(vRunTime, numIter);

	_mm_free(vRunTime);
	_mm_free(mA);
	_mm_free(mB);
	_mm_free(mD);

}

// Auxiliary Functions

void PrintRunTimeData(double *vRunTime, int numIter) {
	int ii;
	double runTime, medianRunTime;

	runTime = 0.0;

	for (ii = 0; ii < numIter; ii++)
	{
		runTime = runTime + vRunTime[ii]; // Total Run Time
	}

	qsort(vRunTime, numIter, sizeof(double), QSortCompFunDouble);

	if ((numIter % 2) == 0) {
		// Even number
		medianRunTime = (vRunTime[numIter / 2] + vRunTime[(numIter / 2) - 1]) / 2.0;
	}
	else
	{
		medianRunTime = vRunTime[(numIter / 2)];
	}

	printf("Total Run Time %2.10f [Sec]\n", runTime);
	printf("Function Minimum Run Time %2.10f [Sec]\n", vRunTime[0]);
	printf("Function Mean Run Time %2.10f [Sec]\n", (runTime / (double)(numIter)));
	printf("Function Median Run Time %2.10f [Sec]\n", medianRunTime);
	printf("Function Maximum Run Time %2.10f [Sec]\n", vRunTime[numIter - 1]);
	printf("\n");
}


int QSortCompFunDouble(const void * a, const void * b) {
	if (*(double*)a < *(double*)b) {
		return -1;
	} 
	else if (*(double*)a >  *(double*)b) {
		return 1;
	}
	else {
		// Namely (*(double*)a == *(double*)b)
		return 0;
	}
}

