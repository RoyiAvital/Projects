#ifndef DLL_TEMPLATE
	#define DLL_TEMPLATE

	#ifdef _WIN32
		#ifdef EXPORT_FCNS
			#define EXPORTED_FUNCTION __declspec(dllexport)
		#else
			#define EXPORTED_FUNCTION __declspec(dllimport)
		#endif
	#else
		#define EXPORTED_FUNCTION
	#endif

	#ifdef  __cplusplus
		extern "C" {
	#endif

	EXPORTED_FUNCTION void ImageConvolution(float* mO, float* mI, int numRows, int numCols, float* mConvKernel, int kernelNumRows, int kernelNumCols);
	EXPORTED_FUNCTION void ImageConvolutionSeparableKernel(float* mO, float* mI, float* mTmp, int numRows, int numCols, float* vRowKernel, int rowKernelLength, float* vColKernel, int colKernelLength);
	EXPORTED_FUNCTION void ImageConvolutionGaussianKernel(float* mO, float* mI, float* mTmp, int numRows, int numCols, float gaussianStd, int stdToRadiusFactor);

	#ifdef  __cplusplus
		}
	#endif

#endif



