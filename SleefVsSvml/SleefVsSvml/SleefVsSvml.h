#ifndef SLEEF_VS_SVML
#define SLEEF_VS_SVML

#define PROJECT_VER 2018_09_24_001

// Loading Standard Libraries
#include <stdint.h>
#include <math.h>
#include <string.h>
#ifdef _OPENMP
#include <omp.h>
#endif // _OPENMP

#ifdef SLEEF_LIB
#define SLEEF_STATIC_LIBS
#define __SSE2__
#define __AVX__
#include "../include/sleef.h"
#endif // SLEEF_LIB

#define OFF 0
#define ON	1

// It doesn't expand to multi line so it doesn't work
//# define LOOP_OPT_PRAGMA \
//#if defined(__ICC) || defined(__ICL) || defined(__INTEL_COMPILER) || defined(_GCC_) || defined(__clang__) \
//__assume_aligned(vO, AVX_ALIGNMENT); \
//__assume_aligned(vI, AVX_ALIGNMENT); \
//__assume(numElements % AVX_ALIGNMENT == 0); \
//#pragma ivdep \
//#endif

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
// #include <zmmintrin.h> // AVX512

#define SSE_4_ISA_FLAG	5
#define AVX_ISA_FLAG	8 // Includes AVX2 & FMA. AVX is fast only on Haswell and Above which is AVX2.

#define SSE_STRIDE_08B		16
#define SSE_STRIDE_16B		8
#define SSE_STRIDE_32B		4
#define SSE_STRIDE_64B		2
#define SSE_STRIDE_DOUBLE	8
#define SSE_STRIDE_TRIPLE	12
#define SSE_STRIDE_QUAD		16
#define SSE_ALIGNMENT		16

#define AVX_STRIDE_08B		32
#define AVX_STRIDE_16B		16
#define AVX_STRIDE_32B		8
#define AVX_STRIDE_64B		4
#define AVX_STRIDE_DOUBLE	16
#define AVX_STRIDE_TRIPLE	24
#define AVX_STRIDE_QUAD		32
#define AVX_ALIGNMENT		32

// SSE - 16 Byte Alignment, AVX - 32 Byte Alignment
#if defined(__ICC) || defined(__ICL) || defined(__INTEL_COMPILER) || defined(_GCC_) || defined(__clang__)
#define DECLARE_ALIGN_PRE 
#define DECLARE_ALIGN_POST __attribute__((aligned(AVX_ALIGNMENT)))
#elif defined(_MSC_VER)
#define DECLARE_ALIGN_PRE __declspec(align(AVX_ALIGNMENT))
#define DECLARE_ALIGN_POST
#else
#define DECLARE_ALIGN_PRE 
#define #define DECLARE_ALIGN_POST
#endif // defined(__ICC) || defined(__ICL) || defined(__INTEL_COMPILER) || defined(_GCC_) || defined(__clang__)

#if defined(__ICC) || defined(__ICL) || defined(__INTEL_COMPILER) || defined(_GCC_) || defined(__clang__)
#define RESTRICT restrict
#elif defined(_MSC_VER)
#define RESTRICT __restrict
#else
#define RESTRICT
#endif // defined(__ICC) || defined(__ICL) || defined(__INTEL_COMPILER) || defined(_GCC_) || defined(__clang__)

#endif // SLEEF_VS_SVML