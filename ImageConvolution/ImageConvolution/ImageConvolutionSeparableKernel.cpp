#ifdef _USRDLL
	#define EXPORT_FCNS
	#include "../ImageConvolutionDll/ImageConvolutionDll.h"
#endif // _USRDLL

#ifndef _USRDLL
	#include "ImageConvolutionMain.h"
#endif // !_USRDLL

#include "ImageConvolution.h"

// ------------------------------- ImageConvolutionSeparableKernel ------------------------------ //
void ImageConvolutionSeparableKernel(float* mO, float* mI, float* mTmp, int numRows, int numCols, float* vRowKernel, int rowKernelLength, float* vColKernel, int colKernelLength)
{
	int ii, jj, kk, pxShift;
	DECLARE_ALIGN float tmpVal[SSE_STRIDE];
	int rowKernelRadius, colKernelRadius, rowSseKernelRadius, colSseKernelRadius;

	__m128 currSum;
	__m128 currPx;
	__m128 kernelWeight;

	rowKernelRadius = rowKernelLength / 2;
	colKernelRadius = colKernelLength / 2;

	if ((rowKernelRadius % SSE_STRIDE)) {
		rowSseKernelRadius = rowKernelRadius + (SSE_STRIDE - (rowKernelRadius % SSE_STRIDE));
	}
	else {
		rowSseKernelRadius = rowKernelRadius;
	}

	if ((colKernelRadius % SSE_STRIDE)) {
		colSseKernelRadius = colKernelRadius + (SSE_STRIDE - (colKernelRadius % SSE_STRIDE));
	}
	else {
		colSseKernelRadius = colKernelRadius;
	}

	/*--- Start - Filtering Rows --- */
	// Unpacking data in Transpose as Pre Processing for filtering along Columns

	/*--- Left Edge Pixels --- */
#pragma omp parallel for private(jj, currSum, kk, pxShift, kernelWeight, currPx, tmpVal)
	for (ii = 0; ii < numRows; ii++) {
		for (jj = 0; jj < rowSseKernelRadius; jj += SSE_STRIDE) {
			currSum = _mm_setzero_ps();
			for (kk = 0; kk < rowKernelLength; kk++) {
				pxShift = kk - rowKernelRadius;
				kernelWeight = _mm_set1_ps(vRowKernel[kk]);

				if ((jj + pxShift) < -2) {
					currPx = _mm_set1_ps(mI[(ii * numCols)]);
				}
				else if ((jj + pxShift) < -1) {
					// currPx = _mm_set_ps(mI[(ii * numCols)], mI[(ii * numCols)], mI[(ii * numCols)], mI[(ii * numCols) + 1]);
					currPx = _mm_set_ps(mI[(ii * numCols) + 1], mI[(ii * numCols)], mI[(ii * numCols)], mI[(ii * numCols)]); // Using set data is packed in reverse compared to load!
				}
				else if ((jj + pxShift) < 0) {
					// currPx = _mm_set_ps(mI[(ii * numCols)], mI[(ii * numCols)], mI[(ii * numCols) + 1], mI[(ii * numCols) + 2]);
					currPx = _mm_set_ps(mI[(ii * numCols) + 2], mI[(ii * numCols) + 1], mI[(ii * numCols)], mI[(ii * numCols)]);
				}
				else {
					currPx = _mm_loadu_ps(&mI[(ii * numCols) + jj + pxShift]);
				}

				currSum = _mm_add_ps(currSum, _mm_mul_ps(kernelWeight, currPx));
			}
			_mm_store_ps(tmpVal, currSum);

			// Unpack Data in Transpose
			for (kk = 0; kk < SSE_STRIDE; kk++) {
				mTmp[((jj + kk) * numRows) + ii] = tmpVal[kk];
			}
		}
	}

	/*--- Main Pixels --- */
#pragma omp parallel for private(jj, currSum, kk, pxShift, kernelWeight, tmpVal)
	for (ii = 0; ii < numRows; ii++) {
		for (jj = rowSseKernelRadius; jj < (numCols - rowSseKernelRadius); jj += SSE_STRIDE) {
			currSum = _mm_setzero_ps();
			for (kk = 0; kk < rowKernelLength; kk++) {
				pxShift = kk - rowKernelRadius;
				kernelWeight = _mm_set1_ps(vRowKernel[kk]);
				//printf("Address %d\n",((int)(&mI[(ii * numCols) + jj + pxShift]) % 16));
				//printf("Address %p\n", &mI[(ii * numCols) + jj + pxShift]);
				currSum = _mm_add_ps(currSum, _mm_mul_ps(kernelWeight, _mm_loadu_ps(&mI[(ii * numCols) + jj + pxShift])));
			}
			_mm_store_ps(tmpVal, currSum);

			// Unpack Data in Transpose
			for (kk = 0; kk < SSE_STRIDE; kk++) {
				mTmp[((jj + kk) * numRows) + ii] = tmpVal[kk];
			}
		}
	}

	/*--- Right Edge Pixels --- */
#pragma omp parallel for private(jj, currSum, kk, pxShift, kernelWeight, currPx, tmpVal)
	for (ii = 0; ii < numRows; ii++) {
		for (jj = (numCols - rowSseKernelRadius); jj < numCols; jj += SSE_STRIDE) {
			currSum = _mm_setzero_ps();
			for (kk = 0; kk < rowKernelLength; kk++) {
				pxShift = kk - rowKernelRadius;
				kernelWeight = _mm_set1_ps(vRowKernel[kk]);

				if ((jj + pxShift) > (numCols - 2)) {
					currPx = _mm_set1_ps(mI[(ii * numCols) + numCols - 1]);
				}
				else if ((jj + pxShift) > (numCols - 3)) {
					//currPx = _mm_set_ps(mI[(ii * numCols) + numCols - 2], mI[(ii * numCols) + numCols - 1], mI[(ii * numCols) + numCols - 1], mI[(ii * numCols) + numCols - 1]);
					currPx = _mm_set_ps(mI[(ii * numCols) + numCols - 1], mI[(ii * numCols) + numCols - 1], mI[(ii * numCols) + numCols - 1], mI[(ii * numCols) + numCols - 2]);
				}
				else if ((jj + pxShift) > (numCols - 4)) {
					// currPx = _mm_set_ps(mI[(ii * numCols) + numCols - 3], mI[(ii * numCols) + numCols - 2], mI[(ii * numCols) + numCols - 1], mI[(ii * numCols) + numCols - 1]);
					currPx = _mm_set_ps(mI[(ii * numCols) + numCols - 1], mI[(ii * numCols) + numCols - 1], mI[(ii * numCols) + numCols - 2], mI[(ii * numCols) + numCols - 3]);
				}
				else {
					currPx = _mm_loadu_ps(&mI[(ii * numCols) + jj + pxShift]);
				}

				currSum = _mm_add_ps(currSum, _mm_mul_ps(kernelWeight, currPx));
			}
			_mm_store_ps(tmpVal, currSum);

			// Unpack Data in Transpose
			for (kk = 0; kk < SSE_STRIDE; kk++) {
				mTmp[((jj + kk) * numRows) + ii] = tmpVal[kk];
			}
		}
	}

	/*--- Finish - Filtering Rows --- */


	/*--- Start - Filtering Columns --- */
	// Loading data from Transposed array for contiguous data

	/*--- Left Edge Pixels --- */
#pragma omp parallel for private(jj, currSum, kk, pxShift, kernelWeight, currPx, tmpVal)
	for (ii = 0; ii < numCols; ii++) {
		for (jj = 0; jj < colSseKernelRadius; jj += SSE_STRIDE) {
			currSum = _mm_setzero_ps();
			for (kk = 0; kk < colKernelLength; kk++) {
				pxShift = kk - colKernelRadius;
				kernelWeight = _mm_set1_ps(vColKernel[kk]);

				if ((jj + pxShift) < -2) {
					currPx = _mm_set1_ps(mTmp[(ii * numRows)]);
				}
				else if ((jj + pxShift) < -1) {
					// currPx = _mm_set_ps(mI[(ii * numRows)], mI[(ii * numRows)], mI[(ii * numRows)], mI[(ii * numRows) + 1]);
					currPx = _mm_set_ps(mTmp[(ii * numRows) + 1], mTmp[(ii * numRows)], mTmp[(ii * numRows)], mTmp[(ii * numRows)]);
				}
				else if ((jj + pxShift) < 0) {
					// currPx = _mm_set_ps(mI[(ii * numRows)], mI[(ii * numRows)], mI[(ii * numRows) + 1], mI[(ii * numRows) + 2]);
					currPx = _mm_set_ps(mTmp[(ii * numRows) + 2], mTmp[(ii * numRows) + 1], mTmp[(ii * numRows)], mTmp[(ii * numRows)]);
				}
				else {
					currPx = _mm_loadu_ps(&mTmp[(ii * numRows) + jj + pxShift]);
				}

				currSum = _mm_add_ps(currSum, _mm_mul_ps(kernelWeight, currPx));
			}
			_mm_store_ps(tmpVal, currSum);

			// Unpack Data in Transpose
			for (kk = 0; kk < SSE_STRIDE; kk++) {
				mO[((jj + kk) * numCols) + ii] = tmpVal[kk];
			}
		}
	}

	/*--- Main Pixels --- */
#pragma omp parallel for private(jj, currSum, kk, pxShift, kernelWeight, tmpVal)
	for (ii = 0; ii < numCols; ii++) {
		for (jj = colSseKernelRadius; jj < (numRows - colSseKernelRadius); jj += SSE_STRIDE) {
			currSum = _mm_setzero_ps();
			for (kk = 0; kk < colKernelLength; kk++) {
				pxShift = kk - colKernelRadius;
				kernelWeight = _mm_set1_ps(vColKernel[kk]);
				currSum = _mm_add_ps(currSum, _mm_mul_ps(kernelWeight, _mm_loadu_ps(&mTmp[(ii * numRows) + jj + pxShift])));
			}
			_mm_store_ps(tmpVal, currSum);

			// Unpack Data in Transpose
			for (kk = 0; kk < SSE_STRIDE; kk++) {
				mO[((jj + kk) * numCols) + ii] = tmpVal[kk];
			}
		}
	}

	/*--- Right Edge Pixels --- */
#pragma omp parallel for private(jj, currSum, kk, pxShift, kernelWeight, currPx, tmpVal)
	for (ii = 0; ii < numCols; ii++) {
		for (jj = (numRows - colSseKernelRadius); jj < numRows; jj += SSE_STRIDE) {
			currSum = _mm_setzero_ps();
			for (kk = 0; kk < colKernelLength; kk++) {
				pxShift = kk - colKernelRadius;
				kernelWeight = _mm_set1_ps(vColKernel[kk]);

				if ((jj + pxShift) > (numRows - 2)) {
					currPx = _mm_set1_ps(mTmp[(ii * numRows) + numRows - 1]);
				}
				else if ((jj + pxShift) >(numRows - 3)) {
					// currPx = _mm_set_ps(mI[(ii * numRows) + numRows - 2], mI[(ii * numRows) + numRows - 1], mI[(ii * numRows) + numRows - 1], mI[(ii * numRows) + numRows - 1]);
					currPx = _mm_set_ps(mTmp[(ii * numRows) + numRows - 1], mTmp[(ii * numRows) + numRows - 1], mTmp[(ii * numRows) + numRows - 1], mTmp[(ii * numRows) + numRows - 2]);
				}
				else if ((jj + pxShift) > (numRows - 4)) {
					//currPx = _mm_set_ps(mI[(ii * numRows) + numRows - 3], mI[(ii * numRows) + numRows - 2], mI[(ii * numRows) + numRows - 1], mI[(ii * numRows) + numRows - 1]);
					currPx = _mm_set_ps(mTmp[(ii * numRows) + numRows - 1], mTmp[(ii * numRows) + numRows - 1], mTmp[(ii * numRows) + numRows - 2], mTmp[(ii * numRows) + numRows - 3]);
				}
				else {
					currPx = _mm_loadu_ps(&mTmp[(ii * numRows) + jj + pxShift]);
				}

				currSum = _mm_add_ps(currSum, _mm_mul_ps(kernelWeight, currPx));
			}
			_mm_store_ps(tmpVal, currSum);

			// Unpack Data in Transpose
			for (kk = 0; kk < SSE_STRIDE; kk++) {
				mO[((jj + kk) * numRows) + ii] = tmpVal[kk];
			}
		}
	}

	
}
