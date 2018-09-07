#include <time.h>
#include "ImageBilateralFilter.h"

#define UNIT_TEST_IMAGE_GAUSSIAN_BLUR 1
#define UNIT_TEST_IMAGE_BILATERAL_FILTER 2

// Main Functions
void ImageConvolution(float* mO, float* mI, int numRows, int numCols, float* mConvKernel, int kernelNumRows, int kernelNumCols);
void ImageConvolutionSeparableKernel(float* mO, float* mI, float* mTmp, int numRows, int numCols, float* vRowKernel, int rowKernelLength, float* vColKernel, int colKernelLength);
void ImageConvolutionGaussianKernel(float* mO, float* mI, float* mTmp, int numRows, int numCols, float gaussianStd, int stdToRadiusFactor);
void BilateralFilterFastCompressive(float* mO, float* mI, int numRows, int numCols, float spatialStd, float rangeStd, int paramK);

// Unit Test Functions 
void ImageConvolutionGaussianKernelUnitTest(int numRows, int numCols, int numIter, float gaussianStd, int stdToRadiusFactor);
void BilateralFilterFastCompressiveUnitTest(int numRows, int numCols, int numIter, float spatialStd, float rangeStd, int paramK);


// Auxiliary Functions
void PrintRunTimeData(double *vRunTime, int numIter);
int QSortCompFunDouble(const void * a, const void * b);