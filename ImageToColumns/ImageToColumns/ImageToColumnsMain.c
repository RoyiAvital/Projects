#include <time.h>
#include <stdlib.h>
#include "ImageToColumnsMain.h"
#include "ImageToColumns.h"

int main() {
	int ii, jj, numRows, numCols, numElements, numIter, blockRadius, blockSize, numPixelsToProcess;
	DATA_TYPE *mO, *mI;
	clock_t clockStart, clockEnd;
	double runTime;

	blockRadius = 5;

	numRows = 2000;
	numCols = 2000;

	numIter = 10;

	numElements			= numRows * numCols;
	blockSize			= (2 * blockRadius) + 1;
	numPixelsToProcess	= (numRows - blockSize + 1) * (numCols - blockSize + 1);

	mO		= (DATA_TYPE*)_mm_malloc(numPixelsToProcess * blockSize * blockSize * sizeof(DATA_TYPE), SSE_ALIGNMENT);
	mI		= (DATA_TYPE*)_mm_malloc(numElements * sizeof(DATA_TYPE), SSE_ALIGNMENT);

	if (mO == NULL) {
		return 10;
	}
	if (mI == NULL) {
		return 10;
	}

	for (ii = 0; ii < numElements; ii++) {
		mI[ii] = (DATA_TYPE)rand();
	}

	clockStart = clock();

	clockStart = clock();
	for (jj = 0; jj < numIter; jj++) {
		ImageToColumns(mO, mI, numRows, numCols, blockRadius);
	}
	clockEnd = clock();

	runTime = (double)(clockEnd - clockStart) / CLOCKS_PER_SEC;

	printf("Run Time %2.10f\n", runTime);

	_mm_free(mO);
	_mm_free(mI);

	getchar();

	return 0;
}

