#include "ImageBilateralFilterMain.h"

void ImageConvolutionGaussianKernelUnitTest(int numRows, int numCols, int numIter, float gaussianStd, int stdToRadiusFactor) {
	int ii;
	int inputSize;
	float *mI;
	float *mO;
	float *mTmp;
	clock_t clockStart, clockEnd;
	double *vRunTime;

	vRunTime = (double*)_mm_malloc(numIter * sizeof(double), SSE_ALIGNMENT);

	inputSize = numRows * numCols;

	mI = (float*)_mm_malloc(inputSize * sizeof(float), SSE_ALIGNMENT);
	mO = (float*)_mm_malloc(inputSize * sizeof(float), SSE_ALIGNMENT);
	mTmp = (float*)_mm_malloc(inputSize * sizeof(float), SSE_ALIGNMENT);

	for (ii = 0; ii < inputSize; ii++) {
		mI[ii] = (float)(rand()) / (float)(RAND_MAX);
		mO[ii] = 0.0f;
	}


	clockStart = clock();

	for (ii = 0; ii < numIter; ii++)
	{
		clockStart = clock();

		ImageConvolutionGaussianKernel(mO, mI, mTmp, numRows, numCols, gaussianStd, stdToRadiusFactor);

		clockEnd = clock();

		vRunTime[ii] = (double)(clockEnd - clockStart) / CLOCKS_PER_SEC;
	}

	printf("\n");
	printf("ImageConvolutionGaussianKernel Unit Test\n");
	PrintRunTimeData(vRunTime, numIter);

	_mm_free(mI);
	_mm_free(mO);
	_mm_free(mTmp);
	_mm_free(vRunTime);

	// getchar();

}

void BilateralFilterFastCompressiveUnitTest(int numRows, int numCols, int numIter, float spatialStd, float rangeStd, int paramK) {
	int ii;
	int inputSize;
	float *mI;
	float *mO;
	clock_t clockStart, clockEnd;
	double *vRunTime;

	vRunTime = (double*)_mm_malloc(numIter * sizeof(double), SSE_ALIGNMENT);

	inputSize	= numRows * numCols;

	mI = (float*)_mm_malloc(inputSize * sizeof(float), SSE_ALIGNMENT);
	mO = (float*)_mm_malloc(inputSize * sizeof(float), SSE_ALIGNMENT);

	for (ii = 0; ii < inputSize; ii++) {
		mI[ii] = (float)(rand()) / (float)(RAND_MAX);
		mO[ii] = 0.0f;
	}


	clockStart = clock();

	for (ii = 0; ii < numIter; ii++)
	{
		clockStart = clock();
		
		BilateralFilterFastCompressive(mO, mI, numRows, numCols, spatialStd, rangeStd, paramK);

		clockEnd = clock();

		vRunTime[ii] = (double)(clockEnd - clockStart) / CLOCKS_PER_SEC;
	}
	
	printf("\n");
	printf("BilateralFilterFastCompressive Unit Test\n");
	PrintRunTimeData(vRunTime, numIter);

	_mm_free(mI);
	_mm_free(mO);
	_mm_free(vRunTime);

	// getchar();

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
	else if (*(double*)a > *(double*)b) {
		return 1;
	}
	else {
		// Namely (*(double*)a == *(double*)b)
		return 0;
	}
}



