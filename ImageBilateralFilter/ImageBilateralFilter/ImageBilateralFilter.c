#ifdef _USRDLL
#define EXPORT_FCNS
#include "ImageBilateralFilterDll.h"
#endif // _USRDLL

#ifndef _USRDLL
#include "ImageBilateralFilterMain.h"
#endif // !_USRDLL

#include "ImageBilateralFilter.h"

#if defined(__ICC) || defined(__ICL) || defined(__INTEL_COMPILER)
#define ENABLE_INTEL_COMPILER_OPTIMIZATIONS ON
#endif
#define ENABLE_GAUSSIAN_BLUR ON

#define M_PIf (float)(M_PI)

void InitOmegaArrays(float* mCOmega, float* mSOmega, float* mI, int numElements, float paramOmega);
void UpdateArrays(float* mO, float* mZ, float* mC, float* mS, float* mCFiltered, float* mSFiltered, int numElements, int iterationIdx, float paramD);
void InitArraysSC(float* mC, float* mS, float* mCOmega, float* mSOmega, int numElements);
void UpdateArraysSC(float* mC, float* mS, float* mCOmega, float* mSOmega, int numElements);
void UpdateOutput(float* mO, float* mZ, float* mI, int numElements, float rangeStd, float paramL);

// ------------------------------- BilateralFilterFastCompressive ------------------------------- //
/*
Applies Gaussian Bilateral Filter using approximation made by Fourier Series.
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
					Assumbed to be a factor of 16 (AVX 512 Stride).
					Structure: Scalar.
					Type : 'Int'.
					Range : {1, 2, ...}.
- spatialStd	-	Spatial STD.
					The STD of the Gaussian Kernel applied on Spatial data.
					Structure: 2D Array.
					Type : 'Float32'.
					Range : (0, inf).
- rangeStd		-	Range STD.
					The STD of the Gaussian Kernel applied on Range data.
					Structure: Scalar.
					Type : 'Float32'.
					Range : (0, inf).
- paramK		-	Parameter K.
					Sets the quality of the aprpxomation.
					Higher number yields better approximation.
					Structure: Scalar.
					Type : 'Int'.
					Range : {2, 3, ...}.
Reference:
1.	Fast Compressive Bilateral Filter (https://ieeexplore.ieee.org/document/7843844/).
	The code is actually a simple improvement of the idea to make it more efficient and faster.
Remarks:
1.	The input size must have 'numCols' which is a multiplication of 16 for vectorization.
2.	Requires CPU with AVX2 support.
3.	It is assumed that 'mO' is initialized with zeros.
TODO :
1.	C
Release Notes:
-	1.0.001	09/09/2018	Royi Avital
	*   Using 'numElements' instead of 'numRows * numCols'.
	*	Fixed typo in the name of 'UpdateOutput'.
-	1.0.000	08/09/2018	Royi Avital
	*   First release version.
*/
// ------------------------------- BilateralFilterFastCompressive ------------------------------- //
void BilateralFilterFastCompressive(float* mO, float* mI, int numRows, int numCols, float spatialStd, float rangeStd, int paramK)
{
	int ii, numElements, paramN;
	float paramL, paramTau, *vParamD, *mZ, *mT, paramOmega, *mCOmega, *mSOmega, *mC, *mS, *mCFiltered, *mSFiltered;

	numElements = numRows * numCols;

	mZ = (float*)_mm_malloc(numElements * sizeof(float), AVX_ALIGNMENT); // Should be initialized to Zero
	mT = (float*)_mm_malloc(numElements * sizeof(float), AVX_ALIGNMENT); // Buffer
	mC = (float*)_mm_malloc(numElements * sizeof(float), AVX_ALIGNMENT);
	mS = (float*)_mm_malloc(numElements * sizeof(float), AVX_ALIGNMENT);
	mCOmega = (float*)_mm_malloc(numElements * sizeof(float), AVX_ALIGNMENT);
	mSOmega = (float*)_mm_malloc(numElements * sizeof(float), AVX_ALIGNMENT);
	mCFiltered = (float*)_mm_malloc(numElements * sizeof(float), AVX_ALIGNMENT);
	mSFiltered = (float*)_mm_malloc(numElements * sizeof(float), AVX_ALIGNMENT);

	memset(mZ, 0.0f, numElements * sizeof(float));

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

	InitOmegaArrays(mCOmega, mSOmega, mI, numElements, paramOmega);

	// Iteration Number 1
	ii = 1;
#if (ENABLE_GAUSSIAN_BLUR == ON)
	 ImageConvolutionGaussianKernel(mCFiltered, mCOmega, mT, numRows, numCols, spatialStd, paramK);
	 ImageConvolutionGaussianKernel(mSFiltered, mSOmega, mT, numRows, numCols, spatialStd, paramK);
#endif
	UpdateArrays(mO, mZ, mCOmega, mSOmega, mCFiltered, mSFiltered, numElements, ii, vParamD[ii - 1]);

	// Iteration Number 2
	ii = 2;
	InitArraysSC(mC, mS, mCOmega, mSOmega, numElements);
#if (ENABLE_GAUSSIAN_BLUR == ON)
	ImageConvolutionGaussianKernel(mCFiltered, mC, mT, numRows, numCols, spatialStd, paramK);
	ImageConvolutionGaussianKernel(mSFiltered, mS, mT, numRows, numCols, spatialStd, paramK);
#endif
	UpdateArrays(mO, mZ, mC, mS, mCFiltered, mSFiltered, numElements, ii, vParamD[ii - 1]);

	for (ii = 3; ii <= paramN; ii++)
	{
		UpdateArraysSC(mC, mS, mCOmega, mSOmega, numElements);
#if (ENABLE_GAUSSIAN_BLUR == ON)
		ImageConvolutionGaussianKernel(mCFiltered, mC, mT, numRows, numCols, spatialStd, paramK);
		ImageConvolutionGaussianKernel(mSFiltered, mS, mT, numRows, numCols, spatialStd, paramK);
#endif
		UpdateArrays(mO, mZ, mC, mS, mCFiltered, mSFiltered, numElements, ii, vParamD[ii - 1]);
	}

	UpdateOutput(mO, mZ, mI, numElements, rangeStd, paramL);

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
void InitOmegaArrays(float* mCOmega, float* mSOmega, float* mI, int numElements, float paramOmega) {

	int ii;

#if (ENABLE_INTEL_COMPILER_OPTIMIZATIONS == ON)

	__m256 m256ParamOmega;
	__m256 m256Tmp;
	__m256 m256Cos;
	__m256 m256Sin;

	m256ParamOmega = _mm256_set1_ps(paramOmega);

#pragma omp parallel for private(m256Tmp, m256Sin, m256Cos)
	for (ii = 0; ii < numElements; ii += AVX_STRIDE)
	{
		m256Tmp = _mm256_mul_ps(m256ParamOmega, _mm256_load_ps(&mI[ii]));
		m256Sin = _mm256_sincos_ps(&m256Cos, m256Tmp);

		_mm256_store_ps(&mCOmega[ii], m256Cos);
		_mm256_store_ps(&mSOmega[ii], m256Sin);
	}

#else

#ifdef __GCC__
#pragma omp parallel for simd
#pragma vector aligned always
#pragma ivdep
#else
#pragma omp parallel for
#endif
	for (ii = 0; ii < numElements; ii++)
	{
		mCOmega[ii] = cosf(paramOmega * mI[ii]);
		mSOmega[ii] = sinf(paramOmega * mI[ii]);
	}

#endif

}


void UpdateArrays(float* mO, float* mZ, float* mC, float* mS, float* mCFiltered, float* mSFiltered, int numElements, int iterationIdx, float paramD) {
	
	int ii;

#if (ENABLE_INTEL_COMPILER_OPTIMIZATIONS == ON)

	__m256 m256IterIdx;
	__m256 m256paramD;

	__m256 m256mO;
	__m256 m256mZ;
	__m256 m256mC;
	__m256 m256mS;
	__m256 m256mCF;
	__m256 m256mSF;

	m256IterIdx = _mm256_set1_ps(iterationIdx); 
	m256paramD = _mm256_set1_ps(paramD);

#pragma omp parallel for private(m256mO, m256mZ, m256mC, m256mS, m256mCF, m256mSF)
	for (ii = 0; ii < numElements; ii += AVX_STRIDE)
	{
		m256mO = _mm256_load_ps(&mO[ii]);
		m256mZ = _mm256_load_ps(&mZ[ii]);
		m256mC = _mm256_load_ps(&mC[ii]);
		m256mS = _mm256_load_ps(&mS[ii]);
		m256mCF = _mm256_load_ps(&mCFiltered[ii]);
		m256mSF = _mm256_load_ps(&mSFiltered[ii]);
		_mm256_store_ps(&mO[ii], _mm256_add_ps(m256mO, _mm256_mul_ps(_mm256_mul_ps(m256IterIdx, m256paramD), _mm256_sub_ps(_mm256_mul_ps(m256mC, m256mSF), _mm256_mul_ps(m256mS, m256mCF)))));
		_mm256_store_ps(&mZ[ii], _mm256_add_ps(m256mZ, _mm256_mul_ps(m256paramD, _mm256_add_ps(_mm256_mul_ps(m256mC, m256mCF), _mm256_mul_ps(m256mS, m256mSF)))));
	}

#else

#ifdef __GCC__
#pragma omp parallel for simd
#pragma vector aligned always
#pragma ivdep
#else
#pragma omp parallel for
#endif
	for (ii = 0; ii < numElements; ii++)
	{
		mO[ii] += (iterationIdx * paramD) * (mC[ii] * mSFiltered[ii] - mS[ii] * mCFiltered[ii]);
		mZ[ii] += paramD * (mC[ii] * mCFiltered[ii] + mS[ii] * mSFiltered[ii]);
	}

#endif

}


void InitArraysSC(float* mC, float* mS, float* mCOmega, float* mSOmega, int numElements) {
	
	int ii;

#if (ENABLE_INTEL_COMPILER_OPTIMIZATIONS == ON)

	__m256 m256Val2_0;
	__m256 m256Val1_0;

	__m256 m256mCO;
	__m256 m256mSO;

	m256Val2_0 = _mm256_set1_ps(2.0f);
	m256Val1_0 = _mm256_set1_ps(1.0f);

#pragma omp parallel for private(m256mCO, m256mSO)
	for (ii = 0; ii < numElements; ii += AVX_STRIDE)
	{
		m256mCO = _mm256_load_ps(&mCOmega[ii]);
		m256mSO = _mm256_load_ps(&mSOmega[ii]);
		_mm256_store_ps(&mS[ii], _mm256_mul_ps(m256Val2_0, _mm256_mul_ps(m256mCO, m256mSO)));
		_mm256_store_ps(&mC[ii], _mm256_sub_ps(_mm256_mul_ps(m256Val2_0, _mm256_mul_ps(m256mCO, m256mCO)), m256Val1_0));
	}

#else

#ifdef __GCC__
#pragma omp parallel for simd
#pragma vector aligned always
#pragma ivdep
#else
#pragma omp parallel for
#endif
	for (ii = 0; ii < numElements; ii++)
	{
		mS[ii] = 2.0f * mCOmega[ii] * mSOmega[ii];
		mC[ii] = 2.0f * mCOmega[ii] * mCOmega[ii] - 1.0f;
	}

#endif

}


void UpdateArraysSC(float* mC, float* mS, float* mCOmega, float* mSOmega, int numElements) {

	int ii;
	float varTmp;

#if (ENABLE_INTEL_COMPILER_OPTIMIZATIONS == ON)

	__m256 m256mC;
	__m256 m256mS;
	__m256 m256mCO;
	__m256 m256mSO;
	__m256 m256mT;

#pragma omp parallel for private(m256mC, m256mS, m256mCO, m256mSO, m256mT)
	for (ii = 0; ii < numElements; ii += AVX_STRIDE)
	{
		m256mC = _mm256_load_ps(&mC[ii]);
		m256mS = _mm256_load_ps(&mS[ii]);
		m256mCO = _mm256_load_ps(&mCOmega[ii]);
		m256mSO = _mm256_load_ps(&mSOmega[ii]);

		m256mT = _mm256_add_ps(_mm256_mul_ps(m256mC, m256mSO), _mm256_mul_ps(m256mS, m256mCO));

		_mm256_store_ps(&mC[ii], _mm256_sub_ps(_mm256_mul_ps(m256mC, m256mCO), _mm256_mul_ps(m256mS, m256mSO)));
		_mm256_store_ps(&mS[ii], m256mT);
	}

#else

#ifdef __GCC__
#pragma omp parallel for simd
#pragma vector aligned always
#pragma ivdep
#else
#pragma omp parallel for
#endif
	for (ii = 0; ii < numElements; ii++)
	{
		varTmp = mC[ii] * mSOmega[ii] + mS[ii] * mCOmega[ii];
		mC[ii] = mC[ii] * mCOmega[ii] - mS[ii] * mSOmega[ii];
		mS[ii] = varTmp;
	}

#endif

}

void UpdateOutput(float* mO, float* mZ, float* mI, int numElements, float rangeStd, float paramL) {

	int ii;

#if (ENABLE_INTEL_COMPILER_OPTIMIZATIONS == ON)

	__m256 m256Val1_0;

	__m256 m256mO;
	__m256 m256mZ;
	__m256 m256mI;
	__m256 m256OutFctr;

	m256Val1_0 = _mm256_set1_ps(1.0f);
	m256OutFctr = _mm256_set1_ps((M_PIf * rangeStd * rangeStd) / paramL);

#pragma omp parallel for private(m256mO, m256mZ, m256mI)
	for (ii = 0; ii < numElements; ii += AVX_STRIDE)
	{
		m256mO = _mm256_load_ps(&mO[ii]);
		m256mZ = _mm256_load_ps(&mZ[ii]);
		m256mI = _mm256_load_ps(&mI[ii]);

		_mm256_store_ps(&mO[ii], _mm256_add_ps(m256mI, _mm256_mul_ps(m256OutFctr, _mm256_div_ps(m256mO, _mm256_add_ps(m256Val1_0, m256mZ)))));
	}

#else

	float outFactor;

	outFactor = (M_PIf * rangeStd * rangeStd) / paramL;

#ifdef __GCC__
#pragma omp parallel for simd
#pragma vector aligned always
#pragma ivdep
#else
#pragma omp parallel for
#endif
	for (ii = 0; ii < numElements; ii++)
	{
		mO[ii] = mI[ii] + (outFactor * (mO[ii] / (1.0f + mZ[ii])));
	}

#endif

}