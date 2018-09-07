#ifdef _USRDLL
	#define EXPORT_FCNS
	#include "ImageBilateralFilterDll.h"
#endif // _USRDLL

#ifndef _USRDLL
	#include "ImageBilateralFilterMain.h"
#endif // !_USRDLL

#include "ImageBilateralFilter.h"

// ------------------------------- ImageConvolutionSeparableKernel ------------------------------ //
/*
Applies Convolution on an Image (2D Array) using given 2D Separable Kernel (Given using 2 1D Kernels).
Input:
- mO			-	Output Image.
					Structure: Image Array (Single Channel).
					Type : 'Single'.
					Range : [0, 1].
- mI			-	Input Image.
					Structure: Image Matrix (Single Channel).
					Type : 'Single'.
					Range : [0, 1].
- mTmp			-	Temporary Image.
					Temporary array used for the intermediate image.
					Structure: Image Matrix (Single Channel).
					Type : 'Single'.
					Range : [0, 1].
- numRows		-	Image Number of Rows.
					Structure: Scalar.
					Type : 'Int'.
					Range : {1, 2, ...}.
- numCols		-	Image Number of Columns.
					Assumbed to be a factor of 4 (SSE Stride).
					Structure: Scalar.
					Type : 'Int'.
					Range : {1, 2, ...}.
- vRowKernel	-	1D Row Convolution Kernel.
					The kernel to applied on the rows of the image.
					Structure: 1D Array.
					Type : 'Single'.
					Range : (-inf, inf).
- rowKernelLen	-	Row Kernel Number of Elementes.
					Structure: Scalar.
					Type : 'Int'.
					Range : {1, 2, ...}.
- vColKernel	-	1D Column Convolution Kernel.
					The kernel to applied on the columns of the image.
					Structure: 1D Array.
					Type : 'Single'.
					Range : (-inf, inf).
- colKernelLen	-	Column Kernel Number of Elementes.
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
// ------------------------------- ImageConvolutionSeparableKernel ------------------------------ //
void ImageConvolutionSeparableKernel(float* mO, float* mI, float* mTmp, int numRows, int numCols, float* vRowKernel, int rowKernelLength, float* vColKernel, int colKernelLength)
{
	int ii, jj, kk, pxShift;
	// DECLARE_ALIGN float tmpVal[SSE_STRIDE];
    DECLARE_ALIGN(float, tmpVal, SSE_STRIDE);
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
					currPx = _mm_set_ps(mI[(ii * numCols) + 1], mI[(ii * numCols)], mI[(ii * numCols)], mI[(ii * numCols)]); // Using `_mm_set_ps` data is packed in reverse compared to `_mm_loadu_ps`!
				}
				else if ((jj + pxShift) < 0) {
					currPx = _mm_set_ps(mI[(ii * numCols) + 2], mI[(ii * numCols) + 1], mI[(ii * numCols)], mI[(ii * numCols)]); // Using `_mm_set_ps` data is packed in reverse compared to `_mm_loadu_ps`!
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
					currPx = _mm_set_ps(mI[(ii * numCols) + numCols - 1], mI[(ii * numCols) + numCols - 1], mI[(ii * numCols) + numCols - 1], mI[(ii * numCols) + numCols - 2]); // Using `_mm_set_ps` data is packed in reverse compared to `_mm_loadu_ps`!
				}
				else if ((jj + pxShift) > (numCols - 4)) {
					currPx = _mm_set_ps(mI[(ii * numCols) + numCols - 1], mI[(ii * numCols) + numCols - 1], mI[(ii * numCols) + numCols - 2], mI[(ii * numCols) + numCols - 3]); // Using `_mm_set_ps` data is packed in reverse compared to `_mm_loadu_ps`!
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
					currPx = _mm_set_ps(mTmp[(ii * numRows) + 1], mTmp[(ii * numRows)], mTmp[(ii * numRows)], mTmp[(ii * numRows)]); // Using `_mm_set_ps` data is packed in reverse compared to `_mm_loadu_ps`!
				}
				else if ((jj + pxShift) < 0) {
					currPx = _mm_set_ps(mTmp[(ii * numRows) + 2], mTmp[(ii * numRows) + 1], mTmp[(ii * numRows)], mTmp[(ii * numRows)]); // Using `_mm_set_ps` data is packed in reverse compared to `_mm_loadu_ps`!
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
					currPx = _mm_set_ps(mTmp[(ii * numRows) + numRows - 1], mTmp[(ii * numRows) + numRows - 1], mTmp[(ii * numRows) + numRows - 1], mTmp[(ii * numRows) + numRows - 2]); // Using `_mm_set_ps` data is packed in reverse compared to `_mm_loadu_ps`!
				}
				else if ((jj + pxShift) > (numRows - 4)) {
					currPx = _mm_set_ps(mTmp[(ii * numRows) + numRows - 1], mTmp[(ii * numRows) + numRows - 1], mTmp[(ii * numRows) + numRows - 2], mTmp[(ii * numRows) + numRows - 3]); // Using `_mm_set_ps` data is packed in reverse compared to `_mm_loadu_ps`!
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
