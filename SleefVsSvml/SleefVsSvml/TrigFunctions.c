#define SLEEF_LIB
#include "SleefVsSvml.h"

#ifdef _USRDLL
#define EXPORT_FCNS
	#include "SleefVsSvmlDll.h"
#endif // _USRDLL
#ifndef _USRDLL

#endif // !_USRDLL

void SineSleefSse(float* vO, float* vI, int numElements)
{

	int ii;
	__m128 elmI;

	for (ii = 0; ii < numElements; ii += SSE_STRIDE_32B) {
		elmI = _mm_loadu_ps(&vI[ii]);
		elmI = Sleef_sinf4_u10sse4(elmI);
		_mm_store_ps(&vO[ii], elmI);
	}


}


void SineSleefAvx(float* vO, float* vI, int numElements)
{

	int ii;
	__m256 elmI;

	for (ii = 0; ii < numElements; ii += AVX_STRIDE_32B) {
		elmI = _mm256_loadu_ps(&vI[ii]);
		elmI = Sleef_sinf8_u10avx2(elmI);
		_mm256_store_ps(&vO[ii], elmI);
	}


}


void SineSvmlSse(float* vO, float* vI, int numElements)
{

	int ii;
	__m128 elmI;

	for (ii = 0; ii < numElements; ii += SSE_STRIDE_32B) {
		elmI = _mm_loadu_ps(&vI[ii]);
		elmI = _mm_sin_ps(elmI);
		_mm_store_ps(&vO[ii], elmI);
	}


}


void SineSvmlAvx(float* vO, float* vI, int numElements)
{

	int ii;
	__m256 elmI;

	for (ii = 0; ii < numElements; ii += AVX_STRIDE_32B) {
		elmI = _mm256_loadu_ps(&vI[ii]);
		elmI = _mm256_sin_ps(elmI);
		_mm256_store_ps(&vO[ii], elmI);
	}


}
