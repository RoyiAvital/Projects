#ifdef _USRDLL
#define EXPORT_FCNS
#include "../ImageConvolutionDll/ImageConvolutionDll.h"
#endif // _USRDLL

#ifndef _USRDLL
#include "ImageConvolutionMain.h"
#endif // !_USRDLL

#include "ImageConvolution.h"

// -------------------------------------- ImageConvolution -------------------------------------- //
void ImageConvolution(float* mO, float* mI, int numRows, int numCols, float* mConvKernel, int kernelNumRows, int kernelNumCols)
{
	int ii, jj, kk, ll, rowShift, colShift;
	int kernelRadius, sseKernelRadius;
	int rowIdx1, rowIdx2, rowIdx3, rowIdx4, colIdx1, colIdx2, colIdx3, colIdx4;

	__m128 currSum;
	__m128 currPx;
	__m128 kernelWeight;

	// Init Parameters

	kernelRadius = kernelNumRows / 2;

	if ((kernelRadius % SSE_STRIDE)) {
		sseKernelRadius = kernelRadius + (SSE_STRIDE - (kernelRadius % SSE_STRIDE));
	}
	else {
		sseKernelRadius = kernelRadius;
	}

	/*--- Top Rows --- */

#pragma omp parallel for private(jj, currSum, currPx, kk, pxShift, kernelWeight, tmpVal)
	for (ii = 0; ii < sseKernelRadius; ii++) {
		/*--- Left Columns --- */
		for (jj = 0; jj < sseKernelRadius; jj += SSE_STRIDE) {
			currSum = _mm_setzero_ps();
			for (kk = 0; kk < kernelNumRows; kk++) {
				rowShift = kk - kernelRadius;

				if ((ii + rowShift) < -2) {
					rowIdx1 = 0;
					rowIdx2 = 0;
					rowIdx3 = 0;
					rowIdx4 = 0;
				}
				else if ((ii + rowShift) < -1) {
					rowIdx1 = 0;
					rowIdx2 = 0;
					rowIdx3 = 0;
					rowIdx4 = 1;
				}
				else if ((ii + rowShift) < 0) {
					rowIdx1 = 0;
					rowIdx2 = 0;
					rowIdx3 = 1;
					rowIdx4 = 2;
				}
				else {
					rowIdx1 = 0;
					rowIdx2 = 1;
					rowIdx3 = 2;
					rowIdx4 = 3;
				}

				for (ll = 0; ll < kernelNumCols; ll++) {
					colShift = ll - kernelRadius;
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
						colIdx1 = 0;
						colIdx2 = 1;
						colIdx3 = 2;
						colIdx4 = 3;
					}

					currPx	= _mm_set_ps(mI[(rowIdx4 * numCols) + colIdx4], mI[(rowIdx3 * numCols) + colIdx3], mI[(rowIdx2 * numCols) + colIdx2], mI[(rowIdx1 * numCols) + colIdx1]); // Using `_mm_set_ps` data is packed in reverse compared to `_mm_loadu_ps`!
					currSum = _mm_add_ps(currSum, _mm_mul_ps(kernelWeight, currPx));
				}

			}
			_mm_store_ps(&mO[(ii * numCols) + jj], currSum);
		}
		/* --- Middle Columns --- */
		for (jj = sseKernelRadius; jj < (numCols - sseKernelRadius); jj += SSE_STRIDE) {
			currSum = _mm_setzero_ps();
			for (kk = 0; kk < kernelNumRows; kk++) {
				rowShift = kk - kernelRadius;

				if ((ii + rowShift) < -2) {
					rowIdx1 = 0;
					rowIdx2 = 0;
					rowIdx3 = 0;
					rowIdx4 = 0;
				}
				else if ((ii + rowShift) < -1) {
					rowIdx1 = 0;
					rowIdx2 = 0;
					rowIdx3 = 0;
					rowIdx4 = 1;
				}
				else if ((ii + rowShift) < 0) {
					rowIdx1 = 0;
					rowIdx2 = 0;
					rowIdx3 = 1;
					rowIdx4 = 2;
				}
				else {
					rowIdx1 = 0;
					rowIdx2 = 1;
					rowIdx3 = 2;
					rowIdx4 = 3;
				}

				for (ll = 0; ll < kernelNumCols; ll++) {
					colShift		= ll - kernelRadius;
					kernelWeight	= _mm_set1_ps(mConvKernel[(kk * kernelNumCols) + ll]);

					currPx	= _mm_set_ps(mI[(rowIdx4 * numCols) + jj + colShift], mI[(rowIdx3 * numCols) + jj + colShift], mI[(rowIdx2 * numCols) + jj + colShift], mI[(rowIdx1 * numCols) + jj + colShift]); // Using `_mm_set_ps` data is packed in reverse compared to `_mm_loadu_ps`!
					currSum = _mm_add_ps(currSum, _mm_mul_ps(kernelWeight, currPx));
				}

			}
			_mm_store_ps(&mO[(ii * numCols) + jj], currSum);
		}
		/* --- Right Columns --- */
		for (jj = (numCols - sseKernelRadius); jj < numCols; jj += SSE_STRIDE) {
			currSum = _mm_setzero_ps();
			for (kk = 0; kk < kernelNumRows; kk++) {
				rowShift = kk - kernelRadius;

				if ((ii + rowShift) < -2) {
					rowIdx1 = 0;
					rowIdx2 = 0;
					rowIdx3 = 0;
					rowIdx4 = 0;
				}
				else if ((ii + rowShift) < -1) {
					rowIdx1 = 0;
					rowIdx2 = 0;
					rowIdx3 = 0;
					rowIdx4 = 1;
				}
				else if ((ii + rowShift) < 0) {
					rowIdx1 = 0;
					rowIdx2 = 0;
					rowIdx3 = 1;
					rowIdx4 = 2;
				}
				else {
					rowIdx1 = 0;
					rowIdx2 = 1;
					rowIdx3 = 2;
					rowIdx4 = 3;
				}

				for (ll = 0; ll < kernelNumCols; ll++) {
					colShift = ll - kernelRadius;
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
						colIdx1 = numCols - 4;
						colIdx2 = numCols - 3;
						colIdx3 = numCols - 2;
						colIdx4 = numCols - 1;
					}

					currPx = _mm_set_ps(mI[(rowIdx4 * numCols) + colIdx4], mI[(rowIdx3 * numCols) + colIdx3], mI[(rowIdx2 * numCols) + colIdx2], mI[(rowIdx1 * numCols) + colIdx1]); // Using `_mm_set_ps` data is packed in reverse compared to `_mm_loadu_ps`!
					currSum = _mm_add_ps(currSum, _mm_mul_ps(kernelWeight, currPx));
				}

			}
			_mm_store_ps(&mO[(ii * numCols) + jj], currSum);
		}
	}

	/*--- Middle Rows --- */
#pragma omp parallel for private(jj, currSum, kk, pxShift, kernelWeight, tmpVal)
	for (ii = sseKernelRadius; ii < (numRows - sseKernelRadius); ii++) {
		/* --- Left Columns --- */
		for (jj = 0; jj < sseKernelRadius; jj += SSE_STRIDE) {
			currSum = _mm_setzero_ps();
			for (kk = 0; kk < kernelNumRows; kk++) {
				rowShift = kk - kernelRadius;

				for (ll = 0; ll < kernelNumCols; ll++) {
					colShift = ll - kernelRadius;
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
						colIdx1 = 0;
						colIdx2 = 1;
						colIdx3 = 2;
						colIdx4 = 3;
					}

					currPx	= _mm_set_ps(mI[((ii + rowShift) * numCols) + colIdx4], mI[((ii + rowShift) * numCols) + colIdx3], mI[((ii + rowShift) * numCols) + colIdx2], mI[((ii + rowShift) * numCols) + colIdx1]); // Using `_mm_set_ps` data is packed in reverse compared to `_mm_loadu_ps`!
					currSum = _mm_add_ps(currSum, _mm_mul_ps(kernelWeight, currPx));
				}

			}
			_mm_store_ps(&mO[(ii * numCols) + jj], currSum);
		}
		/* --- Middle Columns --- */
		for (jj = sseKernelRadius; jj < (numCols - sseKernelRadius); jj += SSE_STRIDE) {
			currSum = _mm_setzero_ps();
			for (kk = 0; kk < kernelNumRows; kk++) {
				rowShift = kk - kernelRadius;
				for (ll = 0; ll < kernelNumCols; ll++) {
					colShift		= ll - kernelRadius;
					kernelWeight	= _mm_set1_ps(mConvKernel[(kk * kernelNumCols) + ll]);
					currSum			= _mm_add_ps(currSum, _mm_mul_ps(kernelWeight, _mm_loadu_ps(&mI[((ii + rowShift) * numCols) + jj + colShift])));
				}

				//printf("Address %d\n",((int)(&mI[(ii * numCols) + jj + pxShift]) % 16));
				//printf("Address %p\n", &mI[(ii * numCols) + jj + pxShift]);

			}
			_mm_store_ps(&mO[(ii * numCols) + jj], currSum);
		}
		/* --- Right Columns --- */
		for (jj = (numCols - sseKernelRadius); jj < numCols; jj += SSE_STRIDE) {
			currSum = _mm_setzero_ps();
			for (kk = 0; kk < kernelNumRows; kk++) {
				rowShift = kk - kernelRadius;

				for (ll = 0; ll < kernelNumCols; ll++) {
					colShift = ll - kernelRadius;
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
						colIdx1 = numCols - 4;
						colIdx2 = numCols - 3;
						colIdx3 = numCols - 2;
						colIdx4 = numCols - 1;
					}

					currPx	= _mm_set_ps(mI[((ii + rowShift) * numCols) + colIdx4], mI[((ii + rowShift) * numCols) + colIdx3], mI[((ii + rowShift) * numCols) + colIdx2], mI[((ii + rowShift) * numCols) + colIdx1]); // Using `_mm_set_ps` data is packed in reverse compared to `_mm_loadu_ps`!
					currSum = _mm_add_ps(currSum, _mm_mul_ps(kernelWeight, currPx));
				}

			}
			_mm_store_ps(&mO[(ii * numCols) + jj], currSum);
		}
	}

	/*--- Bottom Rows --- */
#pragma omp parallel for private(jj, currSum, currPx, kk, pxShift, kernelWeight, tmpVal)
	for (ii = (numRows - sseKernelRadius); ii < numRows; ii++) {
		/*--- Left Columns --- */
		for (jj = 0; jj < sseKernelRadius; jj += SSE_STRIDE) {
			currSum = _mm_setzero_ps();
			for (kk = 0; kk < kernelNumRows; kk++) {
				rowShift = kk - kernelRadius;

				if ((ii + rowShift) > (numRows - 2)) {
					rowIdx1 = numRows - 1;
					rowIdx2 = numRows - 1;
					rowIdx3 = numRows - 1;
					rowIdx4 = numRows - 1;
				}
				else if ((ii + rowShift) > (numRows - 3)) {
					rowIdx1 = numRows - 2;
					rowIdx2 = numRows - 1;
					rowIdx3 = numRows - 1;
					rowIdx4 = numRows - 1;
				}
				else if ((ii + rowShift) > (numRows - 4)) {
					rowIdx1 = numRows - 3;
					rowIdx2 = numRows - 2;
					rowIdx3 = numRows - 1;
					rowIdx4 = numRows - 1;
				}
				else {
					rowIdx1 = numRows - 4;
					rowIdx2 = numRows - 3;
					rowIdx3 = numRows - 2;
					rowIdx4 = numRows - 1;
				}

				for (ll = 0; ll < kernelNumCols; ll++) {
					colShift = ll - kernelRadius;
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
						colIdx1 = 0;
						colIdx2 = 1;
						colIdx3 = 2;
						colIdx4 = 3;
					}

					currPx = _mm_set_ps(mI[(rowIdx4 * numCols) + colIdx4], mI[(rowIdx3 * numCols) + colIdx3], mI[(rowIdx2 * numCols) + colIdx2], mI[(rowIdx1 * numCols) + colIdx1]); // Using `_mm_set_ps` data is packed in reverse compared to `_mm_loadu_ps`!
					currSum = _mm_add_ps(currSum, _mm_mul_ps(kernelWeight, currPx));
				}

			}
			_mm_store_ps(&mO[(ii * numCols) + jj], currSum);
		}
		/* --- Middle Columns --- */
		for (jj = sseKernelRadius; jj < (numCols - sseKernelRadius); jj += SSE_STRIDE) {
			currSum = _mm_setzero_ps();
			for (kk = 0; kk < kernelNumRows; kk++) {
				rowShift = kk - kernelRadius;

				if ((ii + rowShift) >(numRows - 2)) {
					rowIdx1 = numRows - 1;
					rowIdx2 = numRows - 1;
					rowIdx3 = numRows - 1;
					rowIdx4 = numRows - 1;
				}
				else if ((ii + rowShift) > (numRows - 3)) {
					rowIdx1 = numRows - 2;
					rowIdx2 = numRows - 1;
					rowIdx3 = numRows - 1;
					rowIdx4 = numRows - 1;
				}
				else if ((ii + rowShift) > (numRows - 4)) {
					rowIdx1 = numRows - 3;
					rowIdx2 = numRows - 2;
					rowIdx3 = numRows - 1;
					rowIdx4 = numRows - 1;
				}
				else {
					rowIdx1 = numRows - 4;
					rowIdx2 = numRows - 3;
					rowIdx3 = numRows - 2;
					rowIdx4 = numRows - 1;
				}

				for (ll = 0; ll < kernelNumCols; ll++) {
					colShift = ll - kernelRadius;
					kernelWeight = _mm_set1_ps(mConvKernel[(kk * kernelNumCols) + ll]);

					currPx = _mm_set_ps(mI[(rowIdx4 * numCols) + jj + colShift], mI[(rowIdx3 * numCols) + jj + colShift], mI[(rowIdx2 * numCols) + jj + colShift], mI[(rowIdx1 * numCols) + jj + colShift]); // Using `_mm_set_ps` data is packed in reverse compared to `_mm_loadu_ps`!
					currSum = _mm_add_ps(currSum, _mm_mul_ps(kernelWeight, currPx));
				}

			}
			_mm_store_ps(&mO[(ii * numCols) + jj], currSum);
		}
		/* --- Right Columns --- */
		for (jj = (numCols - sseKernelRadius); jj < numCols; jj += SSE_STRIDE) {
			currSum = _mm_setzero_ps();
			for (kk = 0; kk < kernelNumRows; kk++) {
				rowShift = kk - kernelRadius;

				if ((ii + rowShift) >(numRows - 2)) {
					rowIdx1 = numRows - 1;
					rowIdx2 = numRows - 1;
					rowIdx3 = numRows - 1;
					rowIdx4 = numRows - 1;
				}
				else if ((ii + rowShift) > (numRows - 3)) {
					rowIdx1 = numRows - 2;
					rowIdx2 = numRows - 1;
					rowIdx3 = numRows - 1;
					rowIdx4 = numRows - 1;
				}
				else if ((ii + rowShift) > (numRows - 4)) {
					rowIdx1 = numRows - 3;
					rowIdx2 = numRows - 2;
					rowIdx3 = numRows - 1;
					rowIdx4 = numRows - 1;
				}
				else {
					rowIdx1 = numRows - 4;
					rowIdx2 = numRows - 3;
					rowIdx3 = numRows - 2;
					rowIdx4 = numRows - 1;
				}

				for (ll = 0; ll < kernelNumCols; ll++) {
					colShift = ll - kernelRadius;
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
						colIdx1 = numCols - 4;
						colIdx2 = numCols - 3;
						colIdx3 = numCols - 2;
						colIdx4 = numCols - 1;
					}

					currPx = _mm_set_ps(mI[(rowIdx4 * numCols) + colIdx4], mI[(rowIdx3 * numCols) + colIdx3], mI[(rowIdx2 * numCols) + colIdx2], mI[(rowIdx1 * numCols) + colIdx1]); // Using `_mm_set_ps` data is packed in reverse compared to `_mm_loadu_ps`!
					currSum = _mm_add_ps(currSum, _mm_mul_ps(kernelWeight, currPx));
				}

			}
			_mm_store_ps(&mO[(ii * numCols) + jj], currSum);
		}
	}


}
