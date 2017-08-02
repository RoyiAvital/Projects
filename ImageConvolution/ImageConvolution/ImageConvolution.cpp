#ifdef _USRDLL
	#define EXPORT_FCNS
	#include "../ImageConvolutionDll/ImageConvolutionDll.h"
#endif // _USRDLL

#ifndef _USRDLL
	#include "ImageConvolutionMain.h"
#endif // !_USRDLL

#include "ImageConvolution.h"

// --------------------------------------- GaussianBlurSse -------------------------------------- //
void ImageConvolution(float* mO, float* mI, float* mTmp, int numRows, int numCols, float* mConvKernel, int kernelNumRows, int kernelNumCols)
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
	
}
