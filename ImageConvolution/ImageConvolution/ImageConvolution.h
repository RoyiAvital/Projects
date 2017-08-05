#ifndef IMAGE_CONVOLUTION
#define IMAGE_CONVOLUTION

#include <stdio.h>
#include <math.h>
#include <memory.h>
#include <omp.h>

#define OFF 0
#define ON	1

// Intel Intrinsics - https://software.intel.com/sites/landingpage/IntrinsicsGuide
// Header Files for SIMD Intrinsics - http://stackoverflow.com/questions/11228855/header-files-for-simd-intrinsics
// #include <mmintrin.h>  // MMX
// #include <xmmintrin.h> // SSE
// #include <emmintrin.h> // SSE2
// #include <pmmintrin.h> // SSE3
#include <tmmintrin.h> // SSSE3
#ifdef FA_SSE_4
#include <smmintrin.h> // SSE4.1
#include <nmmintrin.h> // SSE4.2
#endif // FA_SSE_4
// #include <ammintrin.h> // SSE4A
// #include <wmmintrin.h> // AES
// #include <immintrin.h> // AVX
// #include <zmmintrin.h> // AVX512

#define SSE_STRIDE 4
#define SSE_ALIGNMENT 128

#define AVX_STRIDE 8
#define AVX_ALIGNMENT 256

// SSE - 16 Byte Alignment, AVX - 32 Byte Alignment
// Limited to 128 (Though documentation states up to 8192)
#define CACHE_ALIGN 16

#ifdef __GNUC__
    #define DECLARE_ALIGN(varType, varName, varSize) varType varName[varSize] __attribute__((aligned (CACHE_ALIGN)));
    // #define DECLARE_ALIGN __declspec(align(CACHE_ALIGN))
#endif

#ifdef _MSC_VER
    #define DECLARE_ALIGN(varType, varName, varSize) __declspec(align(CACHE_ALIGN)) varType varName[varSize];
    // #define DECLARE_ALIGN __declspec(align(CACHE_ALIGN))
#endif

// Macros
// Reference:
//	-	http://www.cprogramming.com/tutorial/cpreprocessor.html
//	-	Agner Fog's Vector Class Library
//	-	libsimdpp (https://github.com/p12tic/libsimdpp)

#define SIGN_SSE(x) _mm_or_ps(_mm_and_ps(_mm_cmpgt_ps((x), _mm_setzero_ps()), _mm_set1_ps(1.0f)), _mm_and_ps(_mm_cmplt_ps((x), _mm_setzero_ps()), _mm_set1_ps(-1.0f)))
#define CMP_GT_SSE(x, y) _mm_and_ps(_mm_cmpgt_ps((x), (y)), _mm_set1_ps(1.0f)) // This is a Floating Point Mask (1.0f for True, 0.0f for false)
#define CMP_GE_SSE(x, y) _mm_and_ps(_mm_cmpge_ps((x), (y)), _mm_set1_ps(1.0f)) // This is a Floating Point Mask (1.0f for True, 0.0f for false)
#define CMP_LT_SSE(x, y) _mm_and_ps(_mm_cmplt_ps((x), (y)), _mm_set1_ps(1.0f)) // This is a Floating Point Mask (1.0f for True, 0.0f for false)
#define CMP_LE_SSE(x, y) _mm_and_ps(_mm_cmple_ps((x), (y)), _mm_set1_ps(1.0f)) // This is a Floating Point Mask (1.0f for True, 0.0f for false)
#define CMP_POS_SSE(x) _mm_and_ps(_mm_cmpgt_ps((x), _mm_setzero_ps()), _mm_set1_ps(1.0f)) // This is a Floating Point Mask (1.0f for True, 0.0f for false)
#define CMP_NEG_SSE(x) _mm_and_ps(_mm_cmplt_ps((x), _mm_setzero_ps()), _mm_set1_ps(1.0f)) // This is a Floating Point Mask (1.0f for True, 0.0f for false)
#define CONV_INTERP_SSE(x, y, mask) _mm_add_ps(_mm_mul_ps((mask), (x)), _mm_mul_ps(_mm_sub_ps(_mm_set1_ps(1), (mask)), (y))) // (mask * x) + ((1 - mask) * y) (Requires Floating Point Mask)
#define COND_SEL_SSE(x, y, mask) _mm_or_ps(_mm_and_ps(x, mask), _mm_andnot_ps(mask, y)) // Conditional Selection x = cond ? a : b; See https://fgiesen.wordpress.com/2016/04/03/sse-mind-the-gap/ (Requires Boolean Mask)
#define HOR_SUM_SSE(x) {x = _mm_hadd_ps(x, x); x = _mm_hadd_ps(x, x);} // See http://stackoverflow.com/questions/6996764/fastest-way-to-do-horizontal-float-vector-sum-on-x86

// Inline Functions

static inline float HorizontalSumSse(__m128 x) {
	// Calculates the sum of SSE Register - http://stackoverflow.com/a/35270026/195787
	__m128 shufReg, sumsReg;
	shufReg = _mm_movehdup_ps(x);        // Broadcast elements 3,1 to 2,0
	sumsReg = _mm_add_ps(x, shufReg);
	shufReg = _mm_movehl_ps(shufReg, sumsReg); // High Half -> Low Half
	sumsReg = _mm_add_ss(sumsReg, shufReg);
	return	_mm_cvtss_f32(sumsReg); // Result in the lower part of the SSE Register
}

static inline float CalcEuclideanDistSquaredSse(__m128 x, __m128 y) {
	__m128 sqrDist, shufReg, sumsReg;
	sqrDist = _mm_sub_ps(x, y);
	sqrDist = _mm_mul_ps(sqrDist, sqrDist);

	// Calculates the sum of SSE Register - http://stackoverflow.com/a/35270026/195787
	shufReg = _mm_movehdup_ps(sqrDist);        // Broadcast elements 3,1 to 2,0
	sumsReg = _mm_add_ps(sqrDist, shufReg);
	shufReg = _mm_movehl_ps(shufReg, sumsReg); // High Half -> Low Half
	sumsReg = _mm_add_ss(sumsReg, shufReg);
	return	_mm_cvtss_f32(sumsReg); // Result in the lower part of the SSE Register
}

static inline float CalcDotProductSse(__m128 x, __m128 y) {
	// SSE Variable Dot Product - http://stackoverflow.com/questions/4120681
	__m128 mulRes, shufReg, sumsReg;
	mulRes = _mm_mul_ps(x, y);

	// Calculates the sum of SSE Register - http://stackoverflow.com/a/35270026/195787
	shufReg = _mm_movehdup_ps(mulRes);        // Broadcast elements 3,1 to 2,0
	sumsReg = _mm_add_ps(mulRes, shufReg);
	shufReg = _mm_movehl_ps(shufReg, sumsReg); // High Half -> Low Half
	sumsReg = _mm_add_ss(sumsReg, shufReg);
	return	_mm_cvtss_f32(sumsReg); // Result in the lower part of the SSE Register
}

#endif