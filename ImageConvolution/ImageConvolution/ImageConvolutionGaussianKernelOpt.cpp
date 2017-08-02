#ifdef _USRDLL
	#define EXPORT_FCNS
	#include "../ImageConvolutionDll/ImageConvolutionDll.h"
#endif // _USRDLL

#ifndef _USRDLL
	#include "ImageConvolutionMain.h"
#endif // !_USRDLL

#include "ImageConvolution.h"

// --------------------------------------- GaussianBlurSse -------------------------------------- //
void ImageConvolutionGaussianKernelOpt(float* mO, float* mI, float* mTmp, int numRows, int numCols, float gaussianStd, int stdToRadiusFactor)
{
	int ii, jj, kk, pxShift;
	DECLARE_ALIGN float tmpVal[SSE_STRIDE];
	float* vKernelArray;
	int kernelRadius, kernelLength;
	float gaussianVar, kernelSum;

	__m128 currSum;
	__m128 currPx;
	__m128 kernelWeight;

	// Init Parameters
	kernelRadius	= (unsigned int)ceilf(stdToRadiusFactor * gaussianStd);
	kernelLength	= (2 * kernelRadius) + 1;
	gaussianVar		= gaussianStd * gaussianStd;

	// Init Kernel Array
	vKernelArray	= (float*)_mm_malloc(kernelLength * sizeof(float), SSE_ALIGNMENT);
	kernelSum		= 0;

	for (ii = 0; ii < kernelLength; ii++) {
		vKernelArray[ii]	= expf(-(ii * ii) / (2 * gaussianVar));
		kernelSum			+= vKernelArray[ii];
	}

	for (ii = 0; ii < kernelLength; ii++) {
		vKernelArray[ii] /= kernelSum;
	}

	/*--- Start - Filtering Columns --- */
	// Unpacking data in Transpose as Pre Processing for filtering along Columns

	/*--- Upper Edge Pixels --- */
#pragma omp parallel for private(jj, currSum, currPx, kk, pxShift, kernelWeight, tmpVal)
	for (ii = 0; ii < kernelRadius; ii++) {
		for (jj = 0; jj < numCols; jj += SSE_STRIDE) {
			currSum = _mm_setzero_ps();
			for (kk = 0; kk < kernelLength; kk++) {
				pxShift			= kk - kernelRadius;
				kernelWeight	= _mm_set1_ps(vKernelArray[kk]);

				if ((ii + pxShift) < 0) {
					currPx = _mm_load_ps(&mI[jj]);
				}
				else {
					currPx = _mm_load_ps(&mI[((ii + pxShift) * numCols) + jj]);
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
	for (ii = kernelRadius; ii < (numRows - kernelRadius); ii++) {
		for (jj = 0; jj < numCols; jj += SSE_STRIDE) {
			currSum = _mm_setzero_ps();
			for (kk = 0; kk < kernelLength; kk++) {
				pxShift = kk - kernelRadius;
				kernelWeight = _mm_set1_ps(vKernelArray[kk]);
				//printf("Address %d\n",((int)(&mI[(ii * numCols) + jj + pxShift]) % 16));
				//printf("Address %p\n", &mI[(ii * numCols) + jj + pxShift]);
				currSum = _mm_add_ps(currSum, _mm_mul_ps(kernelWeight, _mm_load_ps(&mI[((ii + pxShift) * numCols) + jj])));
			}
			_mm_store_ps(tmpVal, currSum);

			// Unpack Data in Transpose
			for (kk = 0; kk < SSE_STRIDE; kk++) {
				mTmp[((jj + kk) * numRows) + ii] = tmpVal[kk];
			}
		}
	}

	/*--- Lower Edge Pixels --- */
#pragma omp parallel for private(jj, currSum, currPx, kk, pxShift, kernelWeight, tmpVal)
	for (ii = (numRows - kernelRadius); ii < numRows; ii++) {
		for (jj = 0; jj < numCols; jj += SSE_STRIDE) {
			currSum = _mm_setzero_ps();
			for (kk = 0; kk < kernelLength; kk++) {
				pxShift = kk - kernelRadius;
				kernelWeight = _mm_set1_ps(vKernelArray[kk]);

				if ((ii + pxShift) > (numRows - 1)) {
					currPx = _mm_load_ps(&mI[((numRows - 1) * numCols) + jj]);
				}
				else {
					currPx = _mm_load_ps(&mI[((ii + pxShift) * numCols) + jj]);
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

	/*--- Upper Edge Pixels --- */
#pragma omp parallel for private(jj, currSum, currPx, kk, pxShift, kernelWeight, tmpVal)
	for (ii = 0; ii < kernelRadius; ii++) {
		for (jj = 0; jj < numCols; jj += SSE_STRIDE) {
			currSum = _mm_setzero_ps();
			for (kk = 0; kk < kernelLength; kk++) {
				pxShift = kk - kernelRadius;
				kernelWeight = _mm_set1_ps(vKernelArray[kk]);

				if ((ii + pxShift) < 0) {
					currPx = _mm_load_ps(&mTmp[jj]);
				}
				else {
					currPx = _mm_load_ps(&mTmp[((ii + pxShift) * numCols) + jj]);
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

	/*--- Main Pixels --- */
#pragma omp parallel for private(jj, currSum, kk, pxShift, kernelWeight, tmpVal)
	for (ii = kernelRadius; ii < (numRows - kernelRadius); ii++) {
		for (jj = 0; jj < numCols; jj += SSE_STRIDE) {
			currSum = _mm_setzero_ps();
			for (kk = 0; kk < kernelLength; kk++) {
				pxShift = kk - kernelRadius;
				kernelWeight = _mm_set1_ps(vKernelArray[kk]);
				//printf("Address %d\n",((int)(&mI[(ii * numCols) + jj + pxShift]) % 16));
				//printf("Address %p\n", &mI[(ii * numCols) + jj + pxShift]);
				currSum = _mm_add_ps(currSum, _mm_mul_ps(kernelWeight, _mm_load_ps(&mTmp[((ii + pxShift) * numCols) + jj])));
			}
			_mm_store_ps(tmpVal, currSum);

			// Unpack Data in Transpose
			for (kk = 0; kk < SSE_STRIDE; kk++) {
				mO[((jj + kk) * numRows) + ii] = tmpVal[kk];
			}
		}
	}

	/*--- Lower Edge Pixels --- */
#pragma omp parallel for private(jj, currSum, currPx, kk, pxShift, kernelWeight, tmpVal)
	for (ii = (numRows - kernelRadius); ii < numRows; ii++) {
		for (jj = 0; jj < numCols; jj += SSE_STRIDE) {
			currSum = _mm_setzero_ps();
			for (kk = 0; kk < kernelLength; kk++) {
				pxShift = kk - kernelRadius;
				kernelWeight = _mm_set1_ps(vKernelArray[kk]);

				if ((ii + pxShift) >(numRows - 1)) {
					currPx = _mm_load_ps(&mTmp[((numRows - 1) * numCols) + jj]);
				}
				else {
					currPx = _mm_load_ps(&mTmp[((ii + pxShift) * numCols) + jj]);
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

	_mm_free(vKernelArray);
}
