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

	EXPORTED_FUNCTION void CalcDistanceMatrixVanilla(float* mD, float* mA, float* mB, int vecDim, int numRowsA, int numRowsB);
	EXPORTED_FUNCTION void CalcDistanceMatrixSse(float* mD, float* mA, float* mB, int vecDim, int numRowsA, int numRowsB);
	EXPORTED_FUNCTION void CalcDistanceMatrixAvx(float* mD, float* mA, float* mB, int vecDim, int numRowsA, int numRowsB);
	EXPORTED_FUNCTION void CalcDistanceMatrixEigen(float* mD, float* mA, float* mB, int vecDim, int numRowsA, int numRowsB);
	EXPORTED_FUNCTION void CalcDistanceMatrixRefTime(float* mD, float* mA, float* mB, int vecDim, int numRowsA, int numRowsB);


#ifdef  __cplusplus
}
#endif

#endif



