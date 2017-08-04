
void ImageConvolution(float* mO, float* mI, int numRows, int numCols, float* mConvKernel, int kernelNumRows, int kernelNumCols);
void ImageConvolutionSeparableKernel(float* mO, float* mI, float* mTmp, int numRows, int numCols, float* vRowKernel, int rowKernelLength, float* vColKernel, int colKernelLength);
void ImageConvolutionGaussianKernel(float* mO, float* mI, float* mTmp, int numRows, int numCols, float gaussianStd, int stdToRadiusFactor);
