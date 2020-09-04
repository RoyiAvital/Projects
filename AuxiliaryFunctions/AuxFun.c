
#include "AuxFun.h"

// Resouecse:
	// * Median Sort - https://github.com/CarlosLunaMota/MedianSort.

unsigned int MaxIdx(ELM_TYPE *RESTRICT vA, unsigned int numElements)
{
#if defined(FA_INTEL_COMPILER)
	// See https://software.intel.com/content/www/us/en/develop/articles/data-alignment-to-assist-vectorization.html
	__assume_aligned(vA, SIMD_ALIGNMENT);
#endif

#if defined(FA_CLANG_COMPILER) || defined(FA__GNU__COMPILER)
	// See https://stackoverflow.com/questions/9608171, https://clang.llvm.org/docs/LanguageExtensions.html
	vA = __builtin_assume_aligned(vA, SIMD_ALIGNMENT);
#endif

	unsigned int ii, maxIdx;
	ELM_TYPE maxVal;

	maxIdx = 0;
	maxVal = vA[0];

	for (ii = 1; ii < numElements; ii++)
	{
		if (vA[ii] > maxVal) {
			maxVal = vA[ii];
			maxIdx = ii;
		}
	}

	return maxIdx;

}

unsigned int MinIdx(ELM_TYPE *RESTRICT vA, unsigned int numElements)
{
#if defined(FA_INTEL_COMPILER)
	// See https://software.intel.com/content/www/us/en/develop/articles/data-alignment-to-assist-vectorization.html
	__assume_aligned(vA, SIMD_ALIGNMENT);
#endif

#if defined(FA_CLANG_COMPILER) || defined(FA__GNU__COMPILER)
	// See https://stackoverflow.com/questions/9608171, https://clang.llvm.org/docs/LanguageExtensions.html
	vA = __builtin_assume_aligned(vA, SIMD_ALIGNMENT);
#endif

	unsigned int ii, minIdx;
	ELM_TYPE minVal;

	minIdx = 0;
	minVal = vA[0];

	for (ii = 1; ii < numElements; ii++)
	{
		if (vA[ii] < minVal) {
			minVal = vA[ii];
			minIdx = ii;
		}
	}

	return minIdx;

}

/*
 *  This Quickselect routine is based on the algorithm described in
 *  "Numerical recipes in C", Second Edition,
 *  Cambridge University Press, 1992, Section 8.5, ISBN 0-521-43108-5
 *  This code by Nicolas Devillard - 1998. Public domain.
 *	See http://ndevilla.free.fr/median/median/index.html
 */

 /*
  *  Returns the Median in the array. It modifies the array!
  */


#define ELEM_SWAP(a,b) { register ELM_TYPE t=(a);(a)=(b);(b)=t; }

ELM_TYPE Median(ELM_TYPE arr[], unsigned int n)
{
	unsigned int low, high;
	unsigned int median;
	unsigned int middle, ll, hh;

	low = 0; high = n - 1; median = (low + high) / 2;
	for (;;) {
		if (high <= low) /* One element only */
			return arr[median];

		if (high == low + 1) {  /* Two elements only */
			if (arr[low] > arr[high])
				ELEM_SWAP(arr[low], arr[high]);
			return arr[median];
		}

		/* Find median of low, middle and high items; swap into position low */
		middle = (low + high) / 2;
		if (arr[middle] > arr[high])    ELEM_SWAP(arr[middle], arr[high]);
		if (arr[low] > arr[high])       ELEM_SWAP(arr[low], arr[high]);
		if (arr[middle] > arr[low])     ELEM_SWAP(arr[middle], arr[low]);

		/* Swap low item (now in position middle) into position (low+1) */
		ELEM_SWAP(arr[middle], arr[low + 1]);

		/* Nibble from each end towards middle, swapping items when stuck */
		ll = low + 1;
		hh = high;
		for (;;) {
			do ll++; while (arr[low] > arr[ll]);
			do hh--; while (arr[hh] > arr[low]);

			if (hh < ll)
				break;

			ELEM_SWAP(arr[ll], arr[hh]);
		}

		/* Swap middle item (in position low) back into correct position */
		ELEM_SWAP(arr[low], arr[hh]);

		/* Re-set active partition */
		if (hh <= median)
			low = ll;
		if (hh >= median)
			high = hh - 1;
	}
}

// https://www.stat.cmu.edu/~ryantibs/median/quickselect.c
ELM_TYPE QuickSelect(ELM_TYPE arr[], unsigned int n, unsigned int k) {
	unsigned int i, ir, j, l, mid;
	ELM_TYPE a;

	l = 0;
	ir = n - 1;
	for (;;) {
		if (ir <= l + 1) {
			if (ir == l + 1 && arr[ir] < arr[l]) {
				ELEM_SWAP(arr[l], arr[ir]);
			}
			return arr[k];
		}
		else {
			mid = (l + ir) >> 1;
			ELEM_SWAP(arr[mid], arr[l + 1]);
			if (arr[l] > arr[ir]) {
				ELEM_SWAP(arr[l], arr[ir]);
			}
			if (arr[l + 1] > arr[ir]) {
				ELEM_SWAP(arr[l + 1], arr[ir]);
			}
			if (arr[l] > arr[l + 1]) {
				ELEM_SWAP(arr[l], arr[l + 1]);
			}
			i = l + 1;
			j = ir;
			a = arr[l + 1];
			for (;;) {
				do i++; while (arr[i] < a);
				do j--; while (arr[j] > a);
				if (j < i) break;
				ELEM_SWAP(arr[i], arr[j]);
			}
			arr[l + 1] = arr[j];
			arr[j] = a;
			if (j >= k) ir = j - 1;
			if (j <= k) l = i;
		}
	}
}

#undef ELEM_SWAP

// -------------------------------------- ApplyGammaCurveFa ------------------------------------- //
/*
Applies Gamma Curve on the Input Image.
Input:
- mI			-	Input Image.
					Structure: Image Matrix (Single Channel).
					Type : 'Single'.
					Range : [0, 1].
- mO			-	Output Image.
					Structure: Image Matrix (Single Channel).
					Type : 'Single'.
					Range : [0, 1].
- numRows		-	Number of Rows.
					Structure: Scalar.
					Type : 'Int'.
					Range : {1, 2, ...}.
- numCols		-	Number of Columns.
					Structure: Scalar.
					Type : 'Int'.
					Range : {1, 2, ...}.
- numColsPad	-	Number of Columns Padded.
					Number of columns including padded elements.
					Must be a factor of 4 (SSE) / 8 (AVX).
					Structure: Scalar.
					Type : 'Int'.
					Range : {1, 2, ...}.
- gammaFactor	-	Gammaa Factor.
					Sets the Power.
					Structure: Scalar.
					Type : 'Single'.
					Range : [0, inf).
Reference:
 1.	See MATLAB's reference implementation.
Remarks:
 1.	Supports "In Place" operation.
 2.	Shouldn't be called if 'gammaFactor == 1' or 'gammaFactor == 0'.
 3. The function is simple and element wise hence should be fused with another such operation.
TODO :
 1.	C
Release Notes:
-	1.2.001	18/09/2018	Royi Avital
	*   Using SVML or Sleef instead of Vector Class (Removing C++ dependency).
-	1.2.000	22/02/2018	Royi Avital
	*   Using Intel Power Function.
	*   Removed loop unrolling.
-	1.1.000	28/10/2017	Royi Avital
	*   Using pointers artihmetic to access data in order to mitigate the multiplication of the index.
	*   Added manual loop unrolling (Factor x2).
 -	1.0.000	12/03/2017	Royi Avital
	*   First release version.
*/
// -------------------------------------- ApplyGammaCurveFa ------------------------------------- //

unsigned int MedianIdx(ELM_TYPE *RESTRICT vA, ELM_TYPE *RESTRICT vB, unsigned int numElements)
{
#if defined(FA_INTEL_COMPILER)
	// See https://software.intel.com/content/www/us/en/develop/articles/data-alignment-to-assist-vectorization.html
	__assume_aligned(vA, SIMD_ALIGNMENT);
#endif

#if defined(FA_CLANG_COMPILER) || defined(FA__GNU__COMPILER)
	// See https://stackoverflow.com/questions/9608171, https://clang.llvm.org/docs/LanguageExtensions.html
	vA = __builtin_assume_aligned(vA, SIMD_ALIGNMENT);
#endif

	unsigned int ii;
	ELM_TYPE medianVal;

	medianVal = Median(vB, numElements);

	for (ii = 0; ii < numElements; ii++)
	{
		if (vA[ii] == medianVal) {
			return ii;
		}
	}

	return 0; // In case of a failure

}


ELM_TYPE CalcMedianVal(ELM_TYPE *RESTRICT vA, unsigned int numElements)
{
#if defined(FA_INTEL_COMPILER)
	// See https://software.intel.com/content/www/us/en/develop/articles/data-alignment-to-assist-vectorization.html
	__assume_aligned(vA, SIMD_ALIGNMENT);
#endif

#if defined(FA_CLANG_COMPILER) || defined(FA__GNU__COMPILER)
	// See https://stackoverflow.com/questions/9608171, https://clang.llvm.org/docs/LanguageExtensions.html
	vA = __builtin_assume_aligned(vA, SIMD_ALIGNMENT);
#endif

	ELM_TYPE *vB;

	vB = (ELM_TYPE*)_mm_malloc(numElements * sizeof(ELM_TYPE), SIMD_ALIGNMENT);

	memcpy(vB, vA, numElements * sizeof(ELM_TYPE));

	if ((numElements % 2) == 0)
	{
		return (QuickSelect(vB, numElements, numElements / 2) + QuickSelect(vB, numElements, (numElements / 2) - 1)) / (ELM_TYPE)2.0;
	}
	else
	{
		return Median(vB, numElements);
	}

}


unsigned int CalcMedianIdx(ELM_TYPE *RESTRICT vA, unsigned int numElements)
{
#if defined(FA_INTEL_COMPILER)
	// See https://software.intel.com/content/www/us/en/develop/articles/data-alignment-to-assist-vectorization.html
	__assume_aligned(vA, SIMD_ALIGNMENT);
#endif

#if defined(FA_CLANG_COMPILER) || defined(FA__GNU__COMPILER)
	// See https://stackoverflow.com/questions/9608171, https://clang.llvm.org/docs/LanguageExtensions.html
	vA = __builtin_assume_aligned(vA, SIMD_ALIGNMENT);
#endif

	unsigned int ii;
	ELM_TYPE medianVal;

	medianVal = CalcMedianVal(vA, numElements);

	for (ii = 0; ii < numElements; ii++)
	{
		if (vA[ii] == medianVal) {
			return ii;
		}
	}

	return 0; // In case of a failure

}

