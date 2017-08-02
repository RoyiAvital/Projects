#include <time.h>
#include <stdlib.h>
#include "ImageConvolutionMain.h"
#include "ImageConvolution.h"

int main() {
	int ii, jj, numRows, numCols, numIter, numElements, stdToRadiusFactor;
	float *mO, *mI, *mTmp;
	float gaussianStd;
	clock_t clockStart, clockEnd;
	double runTime;

	numRows = 5000;
	numCols = 5000;

	numIter = 100;

	gaussianStd			= 5;
	stdToRadiusFactor	= 5;

	numElements = numRows * numCols;

	mO		= (float*)_mm_malloc(numElements * sizeof(float), SSE_ALIGNMENT);
	mI		= (float*)_mm_malloc(numElements * sizeof(float), SSE_ALIGNMENT);
	mTmp	= (float*)_mm_malloc(numElements * sizeof(float), SSE_ALIGNMENT);

	for (ii = 0; ii < numElements; ii++) {
		mI[ii] = (float)rand();
	}

	clockStart = clock();

	clockStart = clock();
	for (jj = 0; jj < numIter; jj++) {
		ImageConvolutionGaussianKernel(mO, mI, mTmp, numRows, numCols, gaussianStd, stdToRadiusFactor);
		// ImageConvolutionGaussianKernelOpt(mO, mI, mTmp, numRows, numCols, gaussianStd, stdToRadiusFactor); // In practice it is slower!
	}
	clockEnd = clock();

	runTime = (double)(clockEnd - clockStart) / CLOCKS_PER_SEC;

	printf("Run Time %2.10f\n", runTime);

	_mm_free(mO);
	_mm_free(mI);
	_mm_free(mTmp);

	getchar();

	return 0;
}

