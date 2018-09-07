#ifdef _USRDLL
#define EXPORT_FCNS
#include "ImageBilateralFilterDll.h"
#endif // _USRDLL

#ifndef _USRDLL
#include "ImageBilateralFilterMain.h"
#endif // !_USRDLL

#include "ImageBilateralFilter.h"

#define M_PIf (float)(M_PI)

void InitOmegaArrays(float* mCOmega, float* mSOmega, float* mI, int numRows, int numCols, float paramOmega);
void UpdateArrays(float* mO, float* mZ, float* mC, float* mS, float* mCFiltered, float* mSFiltered, int numRows, int numCols, int iterationIdx, float paramD);
void InitArraysSC(float* mC, float* mS, float* mCOmega, float* mSOmega, int numRows, int numCols);
void UpdateArraysSC(float* mC, float* mS, float* mT, float* mCOmega, float* mSOmega, int numRows, int numCols);
void UpdtaeOutput(float* mO, float* mZ, float* mI, int numRows, int numCols, float rangeStd, float paramL);

// ------------------------------- BilateralFilterFastCompressive ------------------------------- //
/*
Applies Convolution on an Image (2D Array) using given 2D Kernel.
Input:
- mO			-	Output Image.
					Structure: Image Array (Single Channel).
					Type : 'Single'.
					Range : [0, 1].
- mI			-	Input Image.
					Structure: Image Matrix (Single Channel).
					Type : 'Single'.
					Range : [0, 1].
- numRows		-	Number of Rows.
					Structure: Scalar.
					Type : 'Int'.
					Range : {1, 2, ...}.
- numCols		-	Number of Columns.
					Assumbed to be a factor of 4 (SSE Stride).
					Structure: Scalar.
					Type : 'Int'.
					Range : {1, 2, ...}.
- mConvKernel	-	Convolution Kernel.
					Structure: 2D Array.
					Type : 'Single'.
					Range : (-inf, inf).
- kernelNumRows	-	Kernel Number of Rows.
					Structure: Scalar.
					Type : 'Int'.
					Range : {1, 2, ...}.
- kernelNumCols -	Kernel Number of Columns.
					Assumbed to be a factor of 4 (SSE Stride).
					Structure: Scalar.
					Type : 'Int'.
					Range : {1, 2, ...}.
Reference:
1.	Fast Compressive Bilateral Filter (https://ieeexplore.ieee.org/document/7843844/).
Remarks:
1.	The input size must have 'numCols' which is a multiplication of 16 for vectorization.
2.	Requires CPU with AVX2 support.
3.	It is assumed that 'mO' is initialized with zeros.
TODO :
1.	C
Release Notes:
-	1.0.000	07/09/2017	Royi Avital
*   First release version
*/
// ------------------------------- BilateralFilterFastCompressive ------------------------------- //
void BilateralFilterFastCompressive(float* mO, float* mI, int numRows, int numCols, float spatialStd, float rangeStd, int paramK)
{
	int ii, paramN;
	float paramL, paramTau, *vParamD, *mZ, *mT, paramOmega, *mCOmega, *mSOmega, *mC, *mS, *mCFiltered, *mSFiltered;

	mZ = (float*)_mm_malloc(numRows * numCols * sizeof(float), AVX_ALIGNMENT); // Should be initialized to Zero
	mT = (float*)_mm_malloc(numRows * numCols * sizeof(float), AVX_ALIGNMENT); // Buffer
	mC = (float*)_mm_malloc(numRows * numCols * sizeof(float), AVX_ALIGNMENT);
	mS = (float*)_mm_malloc(numRows * numCols * sizeof(float), AVX_ALIGNMENT);
	mCOmega = (float*)_mm_malloc(numRows * numCols * sizeof(float), AVX_ALIGNMENT);
	mSOmega = (float*)_mm_malloc(numRows * numCols * sizeof(float), AVX_ALIGNMENT);
	mCFiltered = (float*)_mm_malloc(numRows * numCols * sizeof(float), AVX_ALIGNMENT);
	mSFiltered = (float*)_mm_malloc(numRows * numCols * sizeof(float), AVX_ALIGNMENT);

	memset(mZ, 0.0f, numRows * numCols * sizeof(float));

	// Init Parameters

	paramL		= paramK * rangeStd;
	paramTau	= paramK / M_PIf;
	paramN		= ceilf((paramK * paramK) / M_PIf);
	paramOmega	= M_PIf / paramL;

	vParamD = (float*)_mm_malloc(paramN * sizeof(float), AVX_ALIGNMENT);
#ifdef __GCC__
#pragma vector aligned always
#pragma ivdep
#endif
	for (ii = 1; ii <= paramN; ii++)
	{
		vParamD[ii - 1] = 2 * expf(-(ii * ii) / (2 * paramTau * paramTau));
	}

	InitOmegaArrays(mCOmega, mSOmega, mI, numRows, numCols, paramOmega);

	// Iteration Number 1
	ii = 1;
	// ImageConvolutionGaussianKernel(mCFiltered, mCOmega, mT, numRows, numCols, spatialStd, paramK);
	// ImageConvolutionGaussianKernel(mSFiltered, mSOmega, mT, numRows, numCols, spatialStd, paramK);
	UpdateArrays(mO, mZ, mCOmega, mSOmega, mCFiltered, mSFiltered, numRows, numCols, ii, vParamD[ii - 1]);

	// Iteration Number 2
	ii = 2;
	InitArraysSC(mC, mS, mCOmega, mSOmega, numRows, numCols);
	/*ImageConvolutionGaussianKernel(mCFiltered, mC, mT, numRows, numCols, spatialStd, paramK);
	ImageConvolutionGaussianKernel(mSFiltered, mS, mT, numRows, numCols, spatialStd, paramK);*/
	UpdateArrays(mO, mZ, mC, mS, mCFiltered, mSFiltered, numRows, numCols, ii, vParamD[ii - 1]);

	for (ii = 3; ii <= paramN; ii++)
	{
		UpdateArraysSC(mC, mS, mT, mCOmega, mSOmega, numRows, numCols);
		/*ImageConvolutionGaussianKernel(mCFiltered, mC, mT, numRows, numCols, spatialStd, paramK);
		ImageConvolutionGaussianKernel(mSFiltered, mS, mT, numRows, numCols, spatialStd, paramK);*/
		UpdateArrays(mO, mZ, mC, mS, mCFiltered, mSFiltered, numRows, numCols, ii, vParamD[ii - 1]);
	}

	UpdtaeOutput(mO, mZ, mI, numRows, numCols, rangeStd, paramL);

	_mm_free(mZ);
	_mm_free(mT);
	_mm_free(mC);
	_mm_free(mS);
	_mm_free(mCOmega);
	_mm_free(mSOmega);
	_mm_free(mCFiltered);
	_mm_free(mSFiltered);
	_mm_free(vParamD);

}

// Auxiliary Functions
void InitOmegaArrays(float* mCOmega, float* mSOmega, float* mI, int numRows, int numCols, float paramOmega) {

	int ii;

#pragma omp parallel for simd
#pragma vector aligned always
#pragma ivdep
	for (ii = 0; ii < numRows * numCols; ii++)
	{
		mCOmega[ii] = cosf(paramOmega * mI[ii]);
		mSOmega[ii] = sinf(paramOmega * mI[ii]);
	}

}


void UpdateArrays(float* mO, float* mZ, float* mC, float* mS, float* mCFiltered, float* mSFiltered, int numRows, int numCols, int iterationIdx, float paramD) {
	
	int ii;

#pragma omp parallel for simd
#pragma vector aligned always
#pragma ivdep
	for (ii = 0; ii < numRows * numCols; ii++)
	{
		mO[ii] += (iterationIdx * paramD) * (mC[ii] * mSFiltered[ii] - mS[ii] * mCFiltered[ii]);
		mZ[ii] += paramD * (mC[ii] * mCFiltered[ii] + mS[ii] * mSFiltered[ii]);
	}

}


void InitArraysSC(float* mC, float* mS, float* mCOmega, float* mSOmega, int numRows, int numCols) {
	
	int ii;

#pragma omp parallel for simd
#pragma vector aligned always
#pragma ivdep
	for (ii = 0; ii < numRows * numCols; ii++)
	{
		mS[ii] = 2.0f * mCOmega[ii] * mSOmega[ii];
		mC[ii] = 2.0f * mCOmega[ii] * mCOmega[ii] - 1.0f;
	}

}


void UpdateArraysSC(float* mC, float* mS, float* mT, float* mCOmega, float* mSOmega, int numRows, int numCols) {

	int ii;

#pragma omp parallel for simd
#pragma vector aligned always
#pragma ivdep
	for (ii = 0; ii < numRows * numCols; ii++)
	{
		mT[ii] = mC[ii] * mSOmega[ii] + mS[ii] * mCOmega[ii];
		mC[ii] = mC[ii] * mCOmega[ii] - mS[ii] * mSOmega[ii];
		mS[ii] = mT[ii];
	}

}

void UpdtaeOutput(float* mO, float* mZ, float* mI, int numRows, int numCols, float rangeStd, float paramL) {

	int ii;
	float outFactor;

	outFactor = (M_PIf * rangeStd * rangeStd) / paramL;

#pragma omp parallel for simd
#pragma vector aligned always
#pragma ivdep
	for (ii = 0; ii < numRows * numCols; ii++)
	{
		mO[ii] = mI[ii] + (outFactor * (mO[ii] / (1.0f + mZ[ii])));
	}

}