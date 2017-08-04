#ifdef _USRDLL
#define EXPORT_FCNS
#include "../ImageConvolutionDll/ImageConvolutionDll.h"
#endif // _USRDLL

#ifndef _USRDLL
#include "ImageConvolutionMain.h"
#endif // !_USRDLL

#include "ImageConvolution.h"
// -------------------------------------- ImageConvolution -------------------------------------- //
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
1.	[].
Remarks:
1.	Actually applies Correlation and not Convolution.
2.	Using 'Replicate' / 'Nearest Neighbor' Boundary Conditions.
2.	It seems VS 2015 generates faster code than GCC/
TODO :
1.	Why the OpenMP DLL Generated in VS 2015 crashes MATLAB (GCC's DLL Works in MATLAB).
Release Notes:
-	1.0.000	04/08/2017	Royi Avital
*   First release version
*/
// -------------------------------------- ImageConvolution -------------------------------------- //
void ImageConvolution(float* mO, float* mI, int numRows, int numCols, float* mConvKernel, int kernelNumRows, int kernelNumCols)
{
	int ii, jj, kk, ll, rowShift, colShift;
	int kernelRowRadius, kernelColRadius, sseKernelRowRadius, sseKernelColRadius;
	int rowIdx, colIdx1, colIdx2, colIdx3, colIdx4;

	__m128 currSum;
	__m128 currPx;
	__m128 kernelWeight;

	// Init Parameters

	kernelRowRadius = kernelNumRows / 2;
	kernelColRadius = kernelNumCols / 2;

	if ((kernelRowRadius % SSE_STRIDE)) {
		sseKernelRowRadius = kernelRowRadius + (SSE_STRIDE - (kernelRowRadius % SSE_STRIDE));
	}
	else {
		sseKernelRowRadius = kernelRowRadius;
	}

	if ((kernelColRadius % SSE_STRIDE)) {
		sseKernelColRadius = kernelColRadius + (SSE_STRIDE - (kernelColRadius % SSE_STRIDE));
	}
	else {
		sseKernelColRadius = kernelColRadius;
	}

	/*--- Top Rows --- */

#pragma omp parallel for private(jj, currSum, kk, rowShift, rowIdx, ll, colShift, kernelWeight, colIdx1, colIdx2, colIdx3, colIdx4, currPx)
	for (ii = 0; ii < sseKernelRowRadius; ii++) {
		/*--- Left Columns --- */
		for (jj = 0; jj < sseKernelColRadius; jj += SSE_STRIDE) {
			currSum = _mm_setzero_ps();
			for (kk = 0; kk < kernelNumRows; kk++) {
				rowShift	= kk - kernelRowRadius;
				rowIdx		= ((ii + rowShift) < 0 ? 0 : ii + rowShift);

				for (ll = 0; ll < kernelNumCols; ll++) {
					colShift = ll - kernelColRadius;
					kernelWeight = _mm_set1_ps(mConvKernel[(kk * kernelNumCols) + ll]);

					if ((jj + colShift) < -2) {
						colIdx1 = 0;
						colIdx2 = 0;
						colIdx3 = 0;
						colIdx4 = 0;
					}
					else if ((jj + colShift) < -1) {
						colIdx1 = 0;
						colIdx2 = 0;
						colIdx3 = 0;
						colIdx4 = 1;
					}
					else if ((jj + colShift) < 0) {
						colIdx1 = 0;
						colIdx2 = 0;
						colIdx3 = 1;
						colIdx4 = 2;
					}
					else {
						colIdx1 = jj + colShift;
						colIdx2 = jj + colShift + 1;
						colIdx3 = jj + colShift + 2;
						colIdx4 = jj + colShift + 3;
					}

					currPx	= _mm_set_ps(mI[(rowIdx * numCols) + colIdx4], mI[(rowIdx * numCols) + colIdx3], mI[(rowIdx * numCols) + colIdx2], mI[(rowIdx * numCols) + colIdx1]); // Using `_mm_set_ps` data is packed in reverse compared to `_mm_loadu_ps`!
					currSum = _mm_add_ps(currSum, _mm_mul_ps(kernelWeight, currPx));
				}

			}
			_mm_store_ps(&mO[(ii * numCols) + jj], currSum);
		}
		/* --- Middle Columns --- */
		for (jj = sseKernelColRadius; jj < (numCols - sseKernelColRadius); jj += SSE_STRIDE) {
			currSum = _mm_setzero_ps();
			for (kk = 0; kk < kernelNumRows; kk++) {
				rowShift	= kk - kernelRowRadius;
				rowIdx		= ((ii + rowShift) < 0 ? 0 : ii + rowShift);

				for (ll = 0; ll < kernelNumCols; ll++) {
					colShift		= ll - kernelColRadius;
					kernelWeight	= _mm_set1_ps(mConvKernel[(kk * kernelNumCols) + ll]);

					currPx	= _mm_loadu_ps(&mI[(rowIdx * numCols) + jj + colShift]);
					currSum = _mm_add_ps(currSum, _mm_mul_ps(kernelWeight, currPx));
				}

			}
			_mm_store_ps(&mO[(ii * numCols) + jj], currSum);
		}
		/* --- Right Columns --- */
		for (jj = (numCols - sseKernelColRadius); jj < numCols; jj += SSE_STRIDE) {
			currSum = _mm_setzero_ps();
			for (kk = 0; kk < kernelNumRows; kk++) {
				rowShift	= kk - kernelRowRadius;
				rowIdx		= ((ii + rowShift) < 0 ? 0 : ii + rowShift);

				for (ll = 0; ll < kernelNumCols; ll++) {
					colShift		= ll - kernelColRadius;
					kernelWeight	= _mm_set1_ps(mConvKernel[(kk * kernelNumCols) + ll]);

					if ((jj + colShift) > (numCols - 2)) {
						colIdx1 = numCols - 1;
						colIdx2 = numCols - 1;
						colIdx3 = numCols - 1;
						colIdx4 = numCols - 1;
					}
					else if ((jj + colShift) > (numCols - 3)) {
						colIdx1 = numCols - 2;
						colIdx2 = numCols - 1;
						colIdx3 = numCols - 1;
						colIdx4 = numCols - 1;
					}
					else if ((jj + colShift) > (numCols - 4)) {
						colIdx1 = numCols - 3;
						colIdx2 = numCols - 2;
						colIdx3 = numCols - 1;
						colIdx4 = numCols - 1;
					}
					else {
						colIdx1 = jj + colShift;
						colIdx2 = jj + colShift + 1;
						colIdx3 = jj + colShift + 2;
						colIdx4 = jj + colShift + 3;
					}

					currPx = _mm_set_ps(mI[(rowIdx * numCols) + colIdx4], mI[(rowIdx * numCols) + colIdx3], mI[(rowIdx * numCols) + colIdx2], mI[(rowIdx * numCols) + colIdx1]); // Using `_mm_set_ps` data is packed in reverse compared to `_mm_loadu_ps`!
					currSum = _mm_add_ps(currSum, _mm_mul_ps(kernelWeight, currPx));
				}

			}
			_mm_store_ps(&mO[(ii * numCols) + jj], currSum);
		}
	}

	/*--- Middle Rows --- */

#pragma omp parallel for private(jj, currSum, kk, rowShift, rowIdx, ll, colShift, kernelWeight, colIdx1, colIdx2, colIdx3, colIdx4, currPx)
	for (ii = sseKernelRowRadius; ii < (numRows - sseKernelRowRadius); ii++) {
		/* --- Left Columns --- */
		for (jj = 0; jj < sseKernelColRadius; jj += SSE_STRIDE) {
			currSum = _mm_setzero_ps();
			for (kk = 0; kk < kernelNumRows; kk++) {
				rowShift = kk - kernelRowRadius;

				for (ll = 0; ll < kernelNumCols; ll++) {
					colShift = ll - kernelColRadius;
					kernelWeight = _mm_set1_ps(mConvKernel[(kk * kernelNumCols) + ll]);

					if ((jj + colShift) < -2) {
						colIdx1 = 0;
						colIdx2 = 0;
						colIdx3 = 0;
						colIdx4 = 0;
					}
					else if ((jj + colShift) < -1) {
						colIdx1 = 0;
						colIdx2 = 0;
						colIdx3 = 0;
						colIdx4 = 1;
					}
					else if ((jj + colShift) < 0) {
						colIdx1 = 0;
						colIdx2 = 0;
						colIdx3 = 1;
						colIdx4 = 2;
					}
					else {
						colIdx1 = jj + colShift;
						colIdx2 = jj + colShift + 1;
						colIdx3 = jj + colShift + 2;
						colIdx4 = jj + colShift + 3;
					}

					currPx	= _mm_set_ps(mI[((ii + rowShift) * numCols) + colIdx4], mI[((ii + rowShift) * numCols) + colIdx3], mI[((ii + rowShift) * numCols) + colIdx2], mI[((ii + rowShift) * numCols) + colIdx1]); // Using `_mm_set_ps` data is packed in reverse compared to `_mm_loadu_ps`!
					currSum = _mm_add_ps(currSum, _mm_mul_ps(kernelWeight, currPx));
				}

			}
			_mm_store_ps(&mO[(ii * numCols) + jj], currSum);
		}
		/* --- Middle Columns --- */
		for (jj = sseKernelColRadius; jj < (numCols - sseKernelColRadius); jj += SSE_STRIDE) {
			currSum = _mm_setzero_ps();
			for (kk = 0; kk < kernelNumRows; kk++) {
				rowShift = kk - kernelRowRadius;
				for (ll = 0; ll < kernelNumCols; ll++) {
					colShift		= ll - kernelColRadius;
					kernelWeight	= _mm_set1_ps(mConvKernel[(kk * kernelNumCols) + ll]);
					currSum			= _mm_add_ps(currSum, _mm_mul_ps(kernelWeight, _mm_loadu_ps(&mI[((ii + rowShift) * numCols) + jj + colShift])));
				}

				//printf("Address %d\n",((int)(&mI[(ii * numCols) + jj + pxShift]) % 16));
				//printf("Address %p\n", &mI[(ii * numCols) + jj + pxShift]);

			}
			_mm_store_ps(&mO[(ii * numCols) + jj], currSum);
		}
		/* --- Right Columns --- */
		for (jj = (numCols - sseKernelColRadius); jj < numCols; jj += SSE_STRIDE) {
			currSum = _mm_setzero_ps();
			for (kk = 0; kk < kernelNumRows; kk++) {
				rowShift = kk - kernelRowRadius;

				for (ll = 0; ll < kernelNumCols; ll++) {
					colShift = ll - kernelColRadius;
					kernelWeight = _mm_set1_ps(mConvKernel[(kk * kernelNumCols) + ll]);

					if ((jj + colShift) > (numCols - 2)) {
						colIdx1 = numCols - 1;
						colIdx2 = numCols - 1;
						colIdx3 = numCols - 1;
						colIdx4 = numCols - 1;
					}
					else if ((jj + colShift) > (numCols - 3)) {
						colIdx1 = numCols - 2;
						colIdx2 = numCols - 1;
						colIdx3 = numCols - 1;
						colIdx4 = numCols - 1;
					}
					else if ((jj + colShift) > (numCols - 4)) {
						colIdx1 = numCols - 3;
						colIdx2 = numCols - 2;
						colIdx3 = numCols - 1;
						colIdx4 = numCols - 1;
					}
					else {
						colIdx1 = jj + colShift;
						colIdx2 = jj + colShift + 1;
						colIdx3 = jj + colShift + 2;
						colIdx4 = jj + colShift + 3;
					}

					currPx	= _mm_set_ps(mI[((ii + rowShift) * numCols) + colIdx4], mI[((ii + rowShift) * numCols) + colIdx3], mI[((ii + rowShift) * numCols) + colIdx2], mI[((ii + rowShift) * numCols) + colIdx1]); // Using `_mm_set_ps` data is packed in reverse compared to `_mm_loadu_ps`!
					currSum = _mm_add_ps(currSum, _mm_mul_ps(kernelWeight, currPx));
				}

			}
			_mm_store_ps(&mO[(ii * numCols) + jj], currSum);
		}
	}

	/*--- Bottom Rows --- */

#pragma omp parallel for private(jj, currSum, kk, rowShift, rowIdx, ll, colShift, kernelWeight, colIdx1, colIdx2, colIdx3, colIdx4, currPx)
	for (ii = (numRows - sseKernelRowRadius); ii < numRows; ii++) {
		/*--- Left Columns --- */
		for (jj = 0; jj < sseKernelColRadius; jj += SSE_STRIDE) {
			currSum = _mm_setzero_ps();
			for (kk = 0; kk < kernelNumRows; kk++) {
				rowShift	= kk - kernelRowRadius;
				rowIdx		= ((ii + rowShift) > (numRows - 1) ? (numRows - 1) : ii + rowShift);

				for (ll = 0; ll < kernelNumCols; ll++) {
					colShift		= ll - kernelColRadius;
					kernelWeight	= _mm_set1_ps(mConvKernel[(kk * kernelNumCols) + ll]);

					if ((jj + colShift) < -2) {
						colIdx1 = 0;
						colIdx2 = 0;
						colIdx3 = 0;
						colIdx4 = 0;
					}
					else if ((jj + colShift) < -1) {
						colIdx1 = 0;
						colIdx2 = 0;
						colIdx3 = 0;
						colIdx4 = 1;
					}
					else if ((jj + colShift) < 0) {
						colIdx1 = 0;
						colIdx2 = 0;
						colIdx3 = 1;
						colIdx4 = 2;
					}
					else {
						colIdx1 = jj + colShift;
						colIdx2 = jj + colShift + 1;
						colIdx3 = jj + colShift + 2;
						colIdx4 = jj + colShift + 3;
					}

					currPx	= _mm_set_ps(mI[(rowIdx * numCols) + colIdx4], mI[(rowIdx * numCols) + colIdx3], mI[(rowIdx * numCols) + colIdx2], mI[(rowIdx * numCols) + colIdx1]); // Using `_mm_set_ps` data is packed in reverse compared to `_mm_loadu_ps`!
					currSum = _mm_add_ps(currSum, _mm_mul_ps(kernelWeight, currPx));
				}

			}
			_mm_store_ps(&mO[(ii * numCols) + jj], currSum);
		}
		/* --- Middle Columns --- */
		for (jj = sseKernelColRadius; jj < (numCols - sseKernelColRadius); jj += SSE_STRIDE) {
			currSum = _mm_setzero_ps();
			for (kk = 0; kk < kernelNumRows; kk++) {
				rowShift	= kk - kernelRowRadius;
				rowIdx		= ((ii + rowShift) >(numRows - 1) ? (numRows - 1) : ii + rowShift);

				for (ll = 0; ll < kernelNumCols; ll++) {
					colShift		= ll - kernelColRadius;
					kernelWeight	= _mm_set1_ps(mConvKernel[(kk * kernelNumCols) + ll]);

					currPx = _mm_loadu_ps(&mI[(rowIdx * numCols) + jj + colShift]);
					currSum = _mm_add_ps(currSum, _mm_mul_ps(kernelWeight, currPx));
				}

			}
			_mm_store_ps(&mO[(ii * numCols) + jj], currSum);
		}
		/* --- Right Columns --- */
		for (jj = (numCols - sseKernelColRadius); jj < numCols; jj += SSE_STRIDE) {
			currSum = _mm_setzero_ps();
			for (kk = 0; kk < kernelNumRows; kk++) {
				rowShift	= kk - kernelRowRadius;
				rowIdx		= (ii + rowShift >(numRows - 1) ? (numRows - 1) : ii + rowShift);

				for (ll = 0; ll < kernelNumCols; ll++) {
					colShift		= ll - kernelColRadius;
					kernelWeight	= _mm_set1_ps(mConvKernel[(kk * kernelNumCols) + ll]);

					if ((jj + colShift) > (numCols - 2)) {
						colIdx1 = numCols - 1;
						colIdx2 = numCols - 1;
						colIdx3 = numCols - 1;
						colIdx4 = numCols - 1;
					}
					else if ((jj + colShift) > (numCols - 3)) {
						colIdx1 = numCols - 2;
						colIdx2 = numCols - 1;
						colIdx3 = numCols - 1;
						colIdx4 = numCols - 1;
					}
					else if ((jj + colShift) > (numCols - 4)) {
						colIdx1 = numCols - 3;
						colIdx2 = numCols - 2;
						colIdx3 = numCols - 1;
						colIdx4 = numCols - 1;
					}
					else {
						colIdx1 = jj + colShift;
						colIdx2 = jj + colShift + 1;
						colIdx3 = jj + colShift + 2;
						colIdx4 = jj + colShift + 3;
					}

					currPx	= _mm_set_ps(mI[(rowIdx * numCols) + colIdx4], mI[(rowIdx * numCols) + colIdx3], mI[(rowIdx * numCols) + colIdx2], mI[(rowIdx * numCols) + colIdx1]); // Using `_mm_set_ps` data is packed in reverse compared to `_mm_loadu_ps`!
					currSum = _mm_add_ps(currSum, _mm_mul_ps(kernelWeight, currPx));
				}

			}
			_mm_store_ps(&mO[(ii * numCols) + jj], currSum);
		}
	}


}
