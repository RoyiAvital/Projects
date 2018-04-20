#ifndef CALC_DISTANCE_MATRIX
#define CALC_DISTANCE_MATRIX

// Loading Standard Libraries
#include <math.h>
#include <omp.h>

#include "../../../../Libraries/EigenLib/Eigen/Dense"
// #include "Eigen/Dense"

using namespace Eigen;

typedef Map<MatrixXf> EigenMatExt;
typedef Map<VectorXf> EigenVecExt;
typedef MatrixXf EigenMat;
typedef VectorXf EigenVec;

#define OFF 0
#define ON	1

// Intel Intrinsics - https://software.intel.com/sites/landingpage/IntrinsicsGuide
// Header Files for SIMD Intrinsics - http://stackoverflow.com/questions/11228855/header-files-for-simd-intrinsics
// x64 (amd64) Intrinsics List - https://docs.microsoft.com/en-us/cpp/intrinsics/x64-amd64-intrinsics-list
// #include <mmintrin.h>  // MMX
// #include <xmmintrin.h> // SSE
// #include <emmintrin.h> // SSE2
// #include <pmmintrin.h> // SSE3
// #include <tmmintrin.h> // SSSE3
// #include <smmintrin.h> // SSE4.1
// #include <nmmintrin.h> // SSE4.2
// #include <ammintrin.h> // SSE4A
// #include <wmmintrin.h> // AES
#include <immintrin.h> // AVX


#define SSE_STRIDE			4
#define SSE_STRIDE_DOUBLE	8
#define SSE_STRIDE_TRIPLE	12
#define SSE_STRIDE_QUAD		16
#define SSE_ALIGNMENT		128

#define AVX_STRIDE			8
#define AVX_STRIDE_DOUBLE	16
#define AVX_STRIDE_TRIPLE	24
#define AVX_STRIDE_QUAD		32
#define AVX_ALIGNMENT		256

// SSE - 16 Byte Alignment, AVX - 32 Byte Alignment
// Limited to 128 (Though documentation states up to 8192)

#define CACHE_ALIGN 32 // AVX Ready
#define DECLARE_ALIGN __declspec(align(CACHE_ALIGN)) // AVX Ready

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
#define CONV_INTERP_AVX(x, y, mask) _mm256_add_ps(_mm256_mul_ps((mask), (x)), _mm256_mul_ps(_mm256_sub_ps(_mm256_set1_ps(1), (mask)), (y))) // (mask * x) + ((1 - mask) * y) (Requires Floating Point Mask)
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

static inline float HorizontalSumAvx(__m256 x) {
	// Calculates the sum of AVX Register - http://stackoverflow.com/a/35270026/195787
	__m128 lowPack, highPack, shufReg, sumsReg;
	lowPack = _mm256_castps256_ps128(x);
	highPack = _mm256_extractf128_ps(x, 1); // high 128
	lowPack = _mm_add_ps(lowPack, highPack);     // add the low 128

	shufReg = _mm_movehdup_ps(lowPack);        // Broadcast elements 3,1 to 2,0
	sumsReg = _mm_add_ps(lowPack, shufReg);
	shufReg = _mm_movehl_ps(shufReg, sumsReg); // High Half -> Low Half
	sumsReg = _mm_add_ss(sumsReg, shufReg);
	return	_mm_cvtss_f32(sumsReg); // Result in the lower part of the SSE Register
}

static inline float HorizontalMaxSse(__m128 x) {
	// Calculates the max value of SSE Register - http://stackoverflow.com/a/35270026/195787
	__m128 shufReg, sumsReg;
	shufReg = _mm_movehdup_ps(x);        // Broadcast elements 3,1 to 2,0
	sumsReg = _mm_max_ps(x, shufReg);
	shufReg = _mm_movehl_ps(shufReg, sumsReg); // High Half -> Low Half
	sumsReg = _mm_max_ss(sumsReg, shufReg);
	return	_mm_cvtss_f32(sumsReg); // Result in the lower part of the SSE Register
}

static inline float HorizontalMinSse(__m128 x) {
	// Calculates the min value of SSE Register - http://stackoverflow.com/a/35270026/195787
	__m128 shufReg, sumsReg;
	shufReg = _mm_movehdup_ps(x);        // Broadcast elements 3,1 to 2,0
	sumsReg = _mm_min_ps(x, shufReg);
	shufReg = _mm_movehl_ps(shufReg, sumsReg); // High Half -> Low Half
	sumsReg = _mm_min_ss(sumsReg, shufReg);
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

static inline __m128 exp128_ps(__m128 x) {
	// See https://stackoverflow.com/questions/48863719
	__m128   exp_hi = _mm_set1_ps(88.3762626647949f);
	__m128   exp_lo = _mm_set1_ps(-88.3762626647949f);

	__m128   cephes_LOG2EF = _mm_set1_ps(1.44269504088896341f);
	__m128   cephes_exp_C1 = _mm_set1_ps(0.693359375f);
	__m128   cephes_exp_C2 = _mm_set1_ps(-2.12194440e-4f);

	__m128   cephes_exp_p0 = _mm_set1_ps(1.9875691500E-4f);
	__m128   cephes_exp_p1 = _mm_set1_ps(1.3981999507E-3f);
	__m128   cephes_exp_p2 = _mm_set1_ps(8.3334519073E-3f);
	__m128   cephes_exp_p3 = _mm_set1_ps(4.1665795894E-2f);
	__m128   cephes_exp_p4 = _mm_set1_ps(1.6666665459E-1f);
	__m128   cephes_exp_p5 = _mm_set1_ps(5.0000001201E-1f);
	__m128   tmp = _mm_setzero_ps(), fx;
	__m128i  imm0;
	__m128   one = _mm_set1_ps(1.0f);

	x = _mm_min_ps(x, exp_hi);
	x = _mm_max_ps(x, exp_lo);

	/* express exp(x) as exp(g + n*log(2)) */
	fx = _mm_mul_ps(x, cephes_LOG2EF);
	fx = _mm_add_ps(fx, _mm_set1_ps(0.5f));
	tmp = _mm_floor_ps(fx);
	__m128  mask = _mm_cmpgt_ps(tmp, fx);
	mask = _mm_and_ps(mask, one);
	fx = _mm_sub_ps(tmp, mask);
	tmp = _mm_mul_ps(fx, cephes_exp_C1);
	__m128  z = _mm_mul_ps(fx, cephes_exp_C2);
	x = _mm_sub_ps(x, tmp);
	x = _mm_sub_ps(x, z);
	z = _mm_mul_ps(x, x);

	__m128  y = cephes_exp_p0;
	y = _mm_mul_ps(y, x);
	y = _mm_add_ps(y, cephes_exp_p1);
	y = _mm_mul_ps(y, x);
	y = _mm_add_ps(y, cephes_exp_p2);
	y = _mm_mul_ps(y, x);
	y = _mm_add_ps(y, cephes_exp_p3);
	y = _mm_mul_ps(y, x);
	y = _mm_add_ps(y, cephes_exp_p4);
	y = _mm_mul_ps(y, x);
	y = _mm_add_ps(y, cephes_exp_p5);
	y = _mm_mul_ps(y, z);
	y = _mm_add_ps(y, x);
	y = _mm_add_ps(y, one);

	/* build 2^n */
	imm0 = _mm_cvttps_epi32(fx);
	imm0 = _mm_add_epi32(imm0, _mm_set1_epi32(0x7f));
	imm0 = _mm_slli_epi32(imm0, 23);
	__m128  pow2n = _mm_castsi128_ps(imm0);
	y = _mm_mul_ps(y, pow2n);
	return y;
}

static inline __m128 exp128_ps_fast(__m128 x) {
	// See https://stackoverflow.com/questions/48863719
	__m128   exp_hi = _mm_set1_ps(88.3762626647949f);
	__m128   exp_lo = _mm_set1_ps(-88.3762626647949f);

	__m128   cephes_LOG2EF = _mm_set1_ps(1.44269504088896341f);
	__m128   inv_LOG2EF = _mm_set1_ps(0.693147180559945f);

	__m128   cephes_exp_p0 = _mm_set1_ps(1.9875691500E-4f);
	__m128   cephes_exp_p1 = _mm_set1_ps(1.3981999507E-3f);
	__m128   cephes_exp_p2 = _mm_set1_ps(8.3334519073E-3f);
	__m128   cephes_exp_p3 = _mm_set1_ps(4.1665795894E-2f);
	__m128   cephes_exp_p4 = _mm_set1_ps(1.6666665459E-1f);
	__m128   cephes_exp_p5 = _mm_set1_ps(5.0000001201E-1f);
	__m128   fx;
	__m128i  imm0;
	__m128   one = _mm_set1_ps(1.0f);

	x = _mm_min_ps(x, exp_hi);
	x = _mm_max_ps(x, exp_lo);

	/* express exp(x) as exp(g + n*log(2)) */
	fx = _mm_mul_ps(x, cephes_LOG2EF);
	fx = _mm_round_ps(fx, _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC);
	__m128  z = _mm_mul_ps(fx, inv_LOG2EF);
	x = _mm_sub_ps(x, z);
	z = _mm_mul_ps(x, x);

	__m128  y = cephes_exp_p0;
	y = _mm_mul_ps(y, x);
	y = _mm_add_ps(y, cephes_exp_p1);
	y = _mm_mul_ps(y, x);
	y = _mm_add_ps(y, cephes_exp_p2);
	y = _mm_mul_ps(y, x);
	y = _mm_add_ps(y, cephes_exp_p3);
	y = _mm_mul_ps(y, x);
	y = _mm_add_ps(y, cephes_exp_p4);
	y = _mm_mul_ps(y, x);
	y = _mm_add_ps(y, cephes_exp_p5);
	y = _mm_mul_ps(y, z);
	y = _mm_add_ps(y, x);
	y = _mm_add_ps(y, one);

	/* build 2^n */
	imm0 = _mm_cvttps_epi32(fx);
	imm0 = _mm_add_epi32(imm0, _mm_set1_epi32(0x7f));
	imm0 = _mm_slli_epi32(imm0, 23);
	__m128  pow2n = _mm_castsi128_ps(imm0);
	y = _mm_mul_ps(y, pow2n);
	return y;
}

static inline __m256 exp256_ps(__m256 x) {
	// See https://stackoverflow.com/questions/48863719
	__m256   exp_hi = _mm256_set1_ps(88.3762626647949f);
	__m256   exp_lo = _mm256_set1_ps(-88.3762626647949f);

	__m256   cephes_LOG2EF = _mm256_set1_ps(1.44269504088896341f);
	__m256   cephes_exp_C1 = _mm256_set1_ps(0.693359375f);
	__m256   cephes_exp_C2 = _mm256_set1_ps(-2.12194440e-4f);

	__m256   cephes_exp_p0 = _mm256_set1_ps(1.9875691500E-4f);
	__m256   cephes_exp_p1 = _mm256_set1_ps(1.3981999507E-3f);
	__m256   cephes_exp_p2 = _mm256_set1_ps(8.3334519073E-3f);
	__m256   cephes_exp_p3 = _mm256_set1_ps(4.1665795894E-2f);
	__m256   cephes_exp_p4 = _mm256_set1_ps(1.6666665459E-1f);
	__m256   cephes_exp_p5 = _mm256_set1_ps(5.0000001201E-1f);
	__m256   tmp = _mm256_setzero_ps(), fx;
	__m256i  imm0;
	__m256   one = _mm256_set1_ps(1.0f);

	x = _mm256_min_ps(x, exp_hi);
	x = _mm256_max_ps(x, exp_lo);

	/* express exp(x) as exp(g + n*log(2)) */
	fx = _mm256_mul_ps(x, cephes_LOG2EF);
	fx = _mm256_add_ps(fx, _mm256_set1_ps(0.5f));
	tmp = _mm256_floor_ps(fx);
	__m256  mask = _mm256_cmp_ps(tmp, fx, _CMP_GT_OS);
	mask = _mm256_and_ps(mask, one);
	fx = _mm256_sub_ps(tmp, mask);
	tmp = _mm256_mul_ps(fx, cephes_exp_C1);
	__m256  z = _mm256_mul_ps(fx, cephes_exp_C2);
	x = _mm256_sub_ps(x, tmp);
	x = _mm256_sub_ps(x, z);
	z = _mm256_mul_ps(x, x);

	__m256  y = cephes_exp_p0;
	y = _mm256_mul_ps(y, x);
	y = _mm256_add_ps(y, cephes_exp_p1);
	y = _mm256_mul_ps(y, x);
	y = _mm256_add_ps(y, cephes_exp_p2);
	y = _mm256_mul_ps(y, x);
	y = _mm256_add_ps(y, cephes_exp_p3);
	y = _mm256_mul_ps(y, x);
	y = _mm256_add_ps(y, cephes_exp_p4);
	y = _mm256_mul_ps(y, x);
	y = _mm256_add_ps(y, cephes_exp_p5);
	y = _mm256_mul_ps(y, z);
	y = _mm256_add_ps(y, x);
	y = _mm256_add_ps(y, one);

	/* build 2^n */
	imm0 = _mm256_cvttps_epi32(fx);
	imm0 = _mm256_add_epi32(imm0, _mm256_set1_epi32(0x7f));
	imm0 = _mm256_slli_epi32(imm0, 23);
	__m256  pow2n = _mm256_castsi256_ps(imm0);
	y = _mm256_mul_ps(y, pow2n);
	return y;
}

static inline __m256 exp256_ps_fast(__m256 x) {
	// See https://stackoverflow.com/questions/48863719
	__m256   exp_hi = _mm256_set1_ps(88.3762626647949f);
	__m256   exp_lo = _mm256_set1_ps(-88.3762626647949f);

	__m256   cephes_LOG2EF = _mm256_set1_ps(1.44269504088896341f);
	__m256   inv_LOG2EF = _mm256_set1_ps(0.693147180559945f);

	__m256   cephes_exp_p0 = _mm256_set1_ps(1.9875691500E-4f);
	__m256   cephes_exp_p1 = _mm256_set1_ps(1.3981999507E-3f);
	__m256   cephes_exp_p2 = _mm256_set1_ps(8.3334519073E-3f);
	__m256   cephes_exp_p3 = _mm256_set1_ps(4.1665795894E-2f);
	__m256   cephes_exp_p4 = _mm256_set1_ps(1.6666665459E-1f);
	__m256   cephes_exp_p5 = _mm256_set1_ps(5.0000001201E-1f);
	__m256   fx;
	__m256i  imm0;
	__m256   one = _mm256_set1_ps(1.0f);

	x = _mm256_min_ps(x, exp_hi);
	x = _mm256_max_ps(x, exp_lo);

	/* express exp(x) as exp(g + n*log(2)) */
	fx = _mm256_mul_ps(x, cephes_LOG2EF);
	fx = _mm256_round_ps(fx, _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC);
	__m256  z = _mm256_mul_ps(fx, inv_LOG2EF);
	x = _mm256_sub_ps(x, z);
	z = _mm256_mul_ps(x, x);

	__m256  y = cephes_exp_p0;
	y = _mm256_mul_ps(y, x);
	y = _mm256_add_ps(y, cephes_exp_p1);
	y = _mm256_mul_ps(y, x);
	y = _mm256_add_ps(y, cephes_exp_p2);
	y = _mm256_mul_ps(y, x);
	y = _mm256_add_ps(y, cephes_exp_p3);
	y = _mm256_mul_ps(y, x);
	y = _mm256_add_ps(y, cephes_exp_p4);
	y = _mm256_mul_ps(y, x);
	y = _mm256_add_ps(y, cephes_exp_p5);
	y = _mm256_mul_ps(y, z);
	y = _mm256_add_ps(y, x);
	y = _mm256_add_ps(y, one);

	/* build 2^n */
	imm0 = _mm256_cvttps_epi32(fx);
	imm0 = _mm256_add_epi32(imm0, _mm256_set1_epi32(0x7f));
	imm0 = _mm256_slli_epi32(imm0, 23);
	__m256  pow2n = _mm256_castsi256_ps(imm0);
	y = _mm256_mul_ps(y, pow2n);
	return y;
}


#endif