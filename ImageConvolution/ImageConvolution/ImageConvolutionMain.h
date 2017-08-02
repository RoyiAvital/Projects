
void ImageConvolution(float* mO, float* mI, float* mTmp, int numRows, int numCols, float* mConvKernel, int kernelNumRows, int kernelNumCols);
void ImageConvolutionSeparableKernel(float* mO, float* mI, float* mTmp, int numRows, int numCols, float* vRowKernel, int rowKernelNumElements, float* vColKernel, int colKernelNumElements);
void ImageConvolutionGaussianKernel(float* mO, float* mI, float* mTmp, int numRows, int numCols, float gaussianStd, int stdToRadiusFactor);
void ImageConvolutionGaussianKernelOpt(float* mO, float* mI, float* mTmp, int numRows, int numCols, float gaussianStd, int stdToRadiusFactor);
