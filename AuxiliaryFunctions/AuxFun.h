#ifndef AUX_FUN
#define AUX_FUN

// Loading Standard Libraries
#include <stdint.h>
#include <math.h>
#include <string.h>
#ifdef _OPENMP
#include <omp.h>
#endif // _OPENMP

// Detecting Compiler
// MSVC Macros https://docs.microsoft.com/en-us/cpp/preprocessor/predefined-macros
// Intel Compiler Macros https://software.intel.com/en-us/cpp-compiler-developer-guide-and-reference-additional-predefined-macros
#if defined(__INTEL_COMPILER) || defined(__ICL) || defined(__ICC) || defined(__INTEL_LLVM_COMPILER)
#define AUX_FUN_INTEL_COMPILER
#elif defined(__clang__)
#define AUX_FUN_CLANG_COMPILER
#define __assume_aligned(inPtr, alignBytes) inPtr = __builtin_assume_aligned (inPtr, alignBytes) // https://gcc.gnu.org/legacy-ml/gcc-patches/2011-06/msg01876.html
// GNU C / C++ Compiler https://stackoverflow.com/questions/259248 / https://gcc.gnu.org/onlinedocs/cpp/Common-Predefined-Macros.html
#elif defined(__GNUC__) || defined(__GNUG__)
#define AUX_FUN__GNU__COMPILER
#define __assume_aligned(inPtr, alignBytes) inPtr = __builtin_assume_aligned (inPtr, alignBytes) // https://gcc.gnu.org/legacy-ml/gcc-patches/2011-06/msg01876.html
#elif defined(_MSC_VER)
#define AUX_FUN_MSC_COMPILER
#endif


#define FALSE	0
#define TRUE	1

#define OFF 0
#define ON	1

// Intel Intrinsics - https://software.intel.com/sites/landingpage/IntrinsicsGuide
// Header Files for SIMD Intrinsics - http://stackoverflow.com/questions/11228855/header-files-for-simd-intrinsics
// x64 (amd64) Intrinsics List - https://docs.microsoft.com/en-us/cpp/intrinsics/x64-amd64-intrinsics-list
// Compiler MACROS https://stackoverflow.com/questions/28939652
// #include <mmintrin.h>  // MMX
// #include <xmmintrin.h> // SSE
// #include <emmintrin.h> // SSE2
#include <pmmintrin.h> // SSE3 // AMD CPU's do not support SSSE3 unless they supprt AVX (FX line) - https://stackoverflow.com/questions/52858556
// #include <tmmintrin.h> // SSSE3
#ifdef AUX_FUN_SSE_4
#include <smmintrin.h> // SSE4.1
#include <nmmintrin.h> // SSE4.2
#endif // AUX_FUN_SSE_4
// #include <ammintrin.h> // SSE4A
// #include <wmmintrin.h> // AES
#ifdef AUX_FUN_AVX
#include <immintrin.h> // AVX & AVX2
#endif // AUX_FUN_AVX
// #include <zmmintrin.h> // AVX512

#define SSE_STRIDE			4
#define SSE_STRIDE_DOUBLE	8
#define SSE_STRIDE_TRIPLE	12
#define SSE_STRIDE_QUAD		16
#define SSE_ALIGNMENT		16

#define AVX_STRIDE			8
#define AVX_STRIDE_DOUBLE	16
#define AVX_STRIDE_TRIPLE	24
#define AVX_STRIDE_QUAD		32
#define AVX_ALIGNMENT		32

#define SIMD_ALIGNMENT      128 // Future proof up to 1024 Bit SIMD

// SSE - 16 Byte Alignment, AVX - 32 Byte Alignment
// Example: DECLARE_ALIGN_PRE int array[5] DECLARE_ALIGN_POST = { 0, 1, 2, 3, 4 };
// Example: DECLARE_ALIGNED( int array[5] ) = { 0, 1, 2, 3, 4 };
#if defined(AUX_FUN_INTEL_COMPILER) || defined(__GNUC__) || defined(__GNUG__) || defined(__clang__)
#define DECLARE_ALIGN_PRE
#define DECLARE_ALIGN_POST __attribute__((aligned(SIMD_ALIGNMENT)))
#define DECLARE_ALIGNED(ARRAY_DEC) ARRAY_DEC __attribute__((aligned ( SIMD_ALIGNMENT )))
#elif defined(_MSC_VER)
#define DECLARE_ALIGN_PRE __declspec(align(SIMD_ALIGNMENT))
#define DECLARE_ALIGN_POST
#define DECLARE_ALIGNED(ARRAY_DEC) __declspec(align( SIMD_ALIGNMENT )) ARRAY_DEC
#else
#define DECLARE_ALIGN_PRE
#define DECLARE_ALIGN_POST
#define DECLARE_ALIGNED(ARRAY_DEC) ARRAY_DEC
#endif // defined(AUX_FUN_INTEL_COMPILER) || defined(__GNUC__) || defined(__GNUC__) || defined(__clang__)

#if defined(AUX_FUN_INTEL_COMPILER) || defined(__GNUC__) || defined(__GNUC__) || defined(__clang__)
#define RESTRICT restrict
#elif defined(_MSC_VER)
#define RESTRICT __restrict
#else
#define RESTRICT
#endif // defined(AUX_FUN_INTEL_COMPILER) || defined(__GNUC__) || defined(__GNUC__) || defined(__clang__)


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
#define CLAMP_SSE(x, minVal, maxVal) _mm_max_ps(_mm_min_ps(x, maxVal), minVal)

#if defined(AUX_FUN_INTEL_COMPILER)
#define CBRT_SSE(x) _mm_cbrt_ps(x)
#define EXP_SSE(x) _mm_exp_ps(x)
#define POWER_SSE(x, y) _mm_pow_ps(x, y)
#else
#define CBRT_SSE(x) Sleef_cbrtf4_u10sse4(x)
//#define CBRT_SSE(x) Sleef_cbrtf4_u35sse4(x) // 3.5 ULP variants
#define EXP_SSE(x) Sleef_expf4_u10sse4(x)
#define POWER_SSE(x, y) Sleef_powf4_u10sse4(x, y)
#endif // AUX_FUN_INTEL_COMPILER


#ifdef AUX_FUN_AVX
// See https://stackoverflow.com/questions/16988199
#define SIGN_AVX(x) _mm256_or_ps(_mm256_and_ps(_mm256_cmp_ps((x), _mm256_setzero_ps(), _CMP_GT_OQ), _mm256_set1_ps(1.0f)), _mm256_and_ps(_mm256_cmp_ps((x), _mm256_setzero_ps(), _CMP_LT_OQ), _mm256_set1_ps(-1.0f)))
#define CMP_GT_AVX(x, y) _mm256_and_ps(_mm256_cmp_ps((x), (y), _CMP_GT_OQ), _mm256_set1_ps(1.0f)) // This is a Floating Point Mask (1.0f for True, 0.0f for false)
#define CMP_GE_AVX(x, y) _mm256_and_ps(_mm256_cmp_ps((x), (y), _CMP_GE_OQ), _mm256_set1_ps(1.0f)) // This is a Floating Point Mask (1.0f for True, 0.0f for false)
#define CMP_LT_AVX(x, y) _mm256_and_ps(_mm256_cmp_ps((x), (y), _CMP_LT_OQ), _mm256_set1_ps(1.0f)) // This is a Floating Point Mask (1.0f for True, 0.0f for false)
#define CMP_LE_AVX(x, y) _mm256_and_ps(_mm256_cmp_ps((x), (y), _CMP_LE_OQ), _mm256_set1_ps(1.0f)) // This is a Floating Point Mask (1.0f for True, 0.0f for false)
#define CMP_POS_AVX(x) _mm256_and_ps(_mm256_cmp_ps((x), _mm256_setzero_ps(), _CMP_GT_OQ), _mm256_set1_ps(1.0f)) // This is a Floating Point Mask (1.0f for True, 0.0f for false)
#define CMP_NEG_AVX(x) _mm256_and_ps(_mm256_cmp_ps((x), _mm256_setzero_ps(), _CMP_LT_OQ), _mm256_set1_ps(1.0f)) // This is a Floating Point Mask (1.0f for True, 0.0f for false)
#define CONV_INTERP_AVX(x, y, mask) _mm256_add_ps(_mm256_mul_ps((mask), (x)), _mm256_mul_ps(_mm256_sub_ps(_mm256_set1_ps(1), (mask)), (y))) // (mask * x) + ((1 - mask) * y) (Requires Floating Point Mask)
#define CLAMP_AVX(x, minVal, maxVal) _mm256_max_ps(_mm256_min_ps(x, maxVal), minVal)

#if defined(AUX_FUN_INTEL_COMPILER)
#define CBRT_AVX(x) _mm256_cbrt_ps(x)
#define EXP_AVX(x) _mm256_exp_ps(x)
#define POWER_AVX(x, y) _mm256_pow_ps(x, y)
#else
#define CBRT_AVX(x) Sleef_cbrtf8_u10avx2(x)
//#define CBRT_AVX(x) Sleef_cbrtf8_u35avx2(x) // 3.5 ULP variants
#define EXP_AVX(x) Sleef_expf8_u10avx2(x)
#define POWER_AVX(x, y) Sleef_powf8_u10avx2(x, y)
#endif // AUX_FUN_INTEL_COMPILER
#endif //AUX_FUN_AVX

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

static inline __m128 FloorSse(const __m128 x) {
	// http://dss.stephanierct.com/DevBlog/?p=8 (Also in Reference folder)
	__m128i v0 = _mm_setzero_si128();
	__m128i v1 = _mm_cmpeq_epi32(v0, v0);
	__m128i ji = _mm_srli_epi32(v1, 25);
	__m128i tmp = _mm_slli_epi32(ji, 23);
	__m128 j = *(__m128*)(&tmp); //create vector 1.0f
	__m128i i = _mm_cvttps_epi32(x);
	__m128 fi = _mm_cvtepi32_ps(i);
	__m128 igx = _mm_cmpgt_ps(fi, x);
	j = _mm_and_ps(igx, j);
	return _mm_sub_ps(fi, j);
}

static inline __m128 CeilSse(const __m128 x) {
	// http://dss.stephanierct.com/DevBlog/?p=8 (Also in Reference folder)
	__m128i v0 = _mm_setzero_si128();
	__m128i v1 = _mm_cmpeq_epi32(v0, v0);
	__m128i ji = _mm_srli_epi32(v1, 25);
	__m128i tmp = _mm_slli_epi32(ji, 23);
	__m128 j = *(__m128*)(&tmp); //create vector 1.0f
	__m128i i = _mm_cvttps_epi32(x);
	__m128 fi = _mm_cvtepi32_ps(i);
	__m128 igx = _mm_cmplt_ps(fi, x);
	j = _mm_and_ps(igx, j);
	return _mm_add_ps(fi, j);
}

static inline __m128 RoundSse(const __m128 a) {
	// http://dss.stephanierct.com/DevBlog/?p=8 (Also in Reference folder)
	__m128 v0 = _mm_setzero_ps();             //generate the highest value &lt; 2
	__m128 v1 = _mm_cmpeq_ps(v0, v0);
	__m128i tmp = *(__m128i*)(&v1);
	tmp = _mm_srli_epi32(tmp, 2);
	__m128 vNearest2 = *(__m128*)(&tmp);
	__m128i i = _mm_cvttps_epi32(a);
	__m128 aTrunc = _mm_cvtepi32_ps(i);        // truncate a
	__m128 rmd = _mm_sub_ps(a, aTrunc);        // get remainder
	__m128 rmd2 = _mm_mul_ps(rmd, vNearest2); // mul remainder by near 2 will yield the needed offset
	__m128i rmd2i = _mm_cvttps_epi32(rmd2);    // after being truncated of course
	__m128 rmd2Trunc = _mm_cvtepi32_ps(rmd2i);
	__m128 r = _mm_add_ps(aTrunc, rmd2Trunc);
	return r;
}


inline __m128 ModSee(const __m128 a, const __m128 aDiv) {
	// http://dss.stephanierct.com/DevBlog/?p=8 (Also in Reference folder)
	__m128 c = _mm_div_ps(a, aDiv);
	__m128i i = _mm_cvttps_epi32(c);
	__m128 cTrunc = _mm_cvtepi32_ps(i);
	__m128 base = _mm_mul_ps(cTrunc, aDiv);
	__m128 r = _mm_sub_ps(a, base);
	return r;
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
#ifdef AUX_FUN_SSE_4
	tmp = _mm_floor_ps(fx); // SSE 4.1 Instruction
#else
	tmp = FloorSse(fx); // SSE 4.1 Instruction
#endif // AUX_FUN_SSE_4
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
#ifdef AUX_FUN_SSE_4
	fx = _mm_round_ps(fx, _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC); // SSE 4.1 Instruction
#else
	fx = RoundSse(fx); // SSE 4.1 Instruction
#endif // AUX_FUN_SSE_4
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

#ifdef AUX_FUN_AVX

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

#endif //AUX_FUN_AVX

// Auxiliary Function
#ifdef CPU_INSTRUCTION_SET_LIB
static int cpuSupportsSSE2() {
	int32_t reg[4];
	Sleef_x86CpuID(reg, 1, 0);
	return (reg[3] & (1 << 26)) != 0;
}

static int cpuSupportsSSE3() {
	int32_t reg[4];
	Sleef_x86CpuID(reg, 1, 0);
	return (reg[2] & (1 << 0)) != 0;
}

static int cpuSupportsSSE4_1() {
	int32_t reg[4];
	Sleef_x86CpuID(reg, 1, 0);
	return (reg[2] & (1 << 19)) != 0;
}

static int cpuSupportsAVX() {
	int32_t reg[4];
	Sleef_x86CpuID(reg, 1, 0);
	return (reg[2] & (1 << 28)) != 0;
}

static int cpuSupportsFMA4() {
	int32_t reg[4];
	Sleef_x86CpuID(reg, 0x80000001, 0);
	return (reg[3] & (1 << 16)) != 0;
}

static int cpuSupportsAVX2() {
	int32_t reg[4];
	Sleef_x86CpuID(reg, 7, 0);
	return (reg[1] & (1 << 5)) != 0;
}

static int cpuSupportsFMA() {
	int32_t reg[4];
	Sleef_x86CpuID(reg, 1, 0);
	return (reg[2] & (1 << 12)) != 0;
}

static int CpuInstSet() {
	// See helpersse2.h helperavx.h helperavx2.h in Sleef
	if (cpuSupportsAVX2() && cpuSupportsFMA()) {
		return AVX_ISA_FLAG;
	}
	else {
		return SSE_4_ISA_FLAG;
	}
}

int CpuInstSetStatic() {
	static int cpuInsSet = -1;

	if (cpuInsSet != -1) return cpuInsSet;

	// See helpersse2.h helperavx.h helperavx2.h in Sleef
	if (cpuSupportsAVX2() && cpuSupportsFMA()) {
		cpuInsSet = AVX_ISA_FLAG;
	}
	else {
		cpuInsSet = SSE_4_ISA_FLAG;
	}

	return cpuInsSet;
}
#endif // CPU_INSTRUCTION_SET_LIB


#endif