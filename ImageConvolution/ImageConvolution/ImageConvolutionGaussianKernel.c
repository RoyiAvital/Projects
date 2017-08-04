#ifdef _USRDLL
	#define EXPORT_FCNS
	#include "../ImageConvolutionDll/ImageConvolutionDll.h"
#endif // _USRDLL

#ifndef _USRDLL
	#include "ImageConvolutionMain.h"
#endif // !_USRDLL

#include "ImageConvolution.h"

// ------------------------------- ImageConvolutionGaussianKernel ------------------------------- //
void ImageConvolutionGaussianKernel(float* mO, float* mI, float* mTmp, int numRows, int numCols, float gaussianStd, int stdToRadiusFactor)
{
	int ii, jj, kk, pxShift;
	float* vKernelArray;
	int kernelRadius, kernelLength, sseKernelRadius;
	float gaussianVar, kernelSum;

	// Init Parameters
	kernelRadius	= (unsigned int)ceilf(stdToRadiusFactor * gaussianStd);
	kernelLength	= (2 * kernelRadius) + 1;
	gaussianVar		= gaussianStd * gaussianStd;

	// Init Kernel Array
	vKernelArray	= (float*)_mm_malloc(kernelLength * sizeof(float), SSE_ALIGNMENT);

	kernelSum		= 0;

	for (ii = -kernelRadius; ii <= kernelRadius; ii++) {
		vKernelArray[ii + kernelRadius] = expf(-(ii * ii) / (2 * gaussianVar));
		kernelSum += vKernelArray[ii + kernelRadius];
	}

	for (ii = 0; ii < kernelLength; ii++) {
		vKernelArray[ii] /= kernelSum;
	}

	ImageConvolutionSeparableKernel(mO, mI, mTmp, numRows, numCols, vKernelArray, kernelLength, vKernelArray, kernelLength);

}
