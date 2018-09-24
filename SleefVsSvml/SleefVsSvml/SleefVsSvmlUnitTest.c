#include "SleefVsSvmlMain.h"

void UnitTestSine(int numElements, int numIter, float maxVal, int negValFlag) {
	int ii, jj;
	float *vI;
	float *vO;
	void(*funPtr[NUM_FUN_FLAVOR])(float *, float *, int) = { SineSleefSse, SineSvmlSse, SineSleefAvx, SineSvmlAvx };
	char *funNameSleefSse[] = {"Sine SLEEF SSE Unit Test\n\0"};
	char *funNameSleefAvx[] = { "Sine SVML SSE Unit Test\n\0" };
	char *funNameSvmlSse[] = { "Sine SLEEF SSE Unit Test\n\0" };
	char *funNameSvmlAvx[] = { "Sine SVML AVX Unit Test\n\0" };
	char **funName[NUM_FUN_FLAVOR] = { funNameSleefSse, funNameSleefAvx, funNameSvmlSse, funNameSvmlAvx };
	clock_t clockStart, clockEnd;
	double *vRunTime;

	vRunTime = (double*)_mm_malloc(numIter * sizeof(double), SSE_ALIGNMENT);

	vI = (float*)_mm_malloc(numElements * sizeof(float), SSE_ALIGNMENT);
	vO = (float*)_mm_malloc(numElements * sizeof(float), SSE_ALIGNMENT);

	for (ii = 0; ii < numElements; ii++) {
		if (negValFlag == ON) {
			vI[ii] = 2.0f * maxVal * ((float)(rand()) / (float)(RAND_MAX) - 0.5f);
		}
		else {
			vI[ii] = maxVal * (float)(rand()) / (float)(RAND_MAX);
		}
	}

	for (jj = 0; jj < NUM_FUN_FLAVOR; jj++)
	{
		clockStart = clock();

		for (ii = 0; ii < numIter; ii++)
		{
			clockStart = clock();

			(*funPtr[jj])(vO, vI, numElements);

			clockEnd = clock();

			vRunTime[ii] = (double)(clockEnd - clockStart) / CLOCKS_PER_SEC;
		}

		printf("\n");
		printf("%s", *funName[jj]);
		PrintRunTimeData(vRunTime, numIter);
	}

	_mm_free(vI);
	_mm_free(vO);
	_mm_free(vRunTime);

}


void UnitTestExp(int numElements, int numIter, float maxVal, int negValFlag) {
	int ii, jj;
	float *vI;
	float *vO;
	void(*funPtr[NUM_FUN_FLAVOR])(float *, float *, int) = { ExpSleefSse, ExpSvmlSse, ExpSleefAvx, ExpSvmlAvx };
	char *funNameSleefSse[] = { "Exp SLEEF SSE Unit Test\n\0" };
	char *funNameSleefAvx[] = { "Exp SVML SSE Unit Test\n\0" };
	char *funNameSvmlSse[] = { "Exp SLEEF SSE Unit Test\n\0" };
	char *funNameSvmlAvx[] = { "Exp SVML AVX Unit Test\n\0" };
	char **funName[NUM_FUN_FLAVOR] = { funNameSleefSse, funNameSleefAvx, funNameSvmlSse, funNameSvmlAvx };
	clock_t clockStart, clockEnd;
	double *vRunTime;

	vRunTime = (double*)_mm_malloc(numIter * sizeof(double), SSE_ALIGNMENT);

	vI = (float*)_mm_malloc(numElements * sizeof(float), SSE_ALIGNMENT);
	vO = (float*)_mm_malloc(numElements * sizeof(float), SSE_ALIGNMENT);

	for (ii = 0; ii < numElements; ii++) {
		if (negValFlag == ON) {
			vI[ii] = 2.0f * maxVal * ((float)(rand()) / (float)(RAND_MAX)-0.5f);
		}
		else {
			vI[ii] = maxVal * (float)(rand()) / (float)(RAND_MAX);
		}
	}

	for (jj = 0; jj < NUM_FUN_FLAVOR; jj++)
	{
		clockStart = clock();

		for (ii = 0; ii < numIter; ii++)
		{
			clockStart = clock();

			(*funPtr[jj])(vO, vI, numElements);

			clockEnd = clock();

			vRunTime[ii] = (double)(clockEnd - clockStart) / CLOCKS_PER_SEC;
		}

		printf("\n");
		printf("%s", *funName[jj]);
		PrintRunTimeData(vRunTime, numIter);
	}

	_mm_free(vI);
	_mm_free(vO);
	_mm_free(vRunTime);

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

