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

	EXPORTED_FUNCTION void SineSleefSse(float* vO, float* vI, int numElements);
	EXPORTED_FUNCTION void SineSleefAvx(float* vO, float* vI, int numElements);
	EXPORTED_FUNCTION void SineSvmlSse(float* vO, float* vI, int numElements);
	EXPORTED_FUNCTION void SineSvmlAvx(float* vO, float* vI, int numElements);
	//EXPORTED_FUNCTION void CosineSleefSse(float* vO, float* vI, int numElements);
	//EXPORTED_FUNCTION void CosineSleefAvx(float* vO, float* vI, int numElements);
	//EXPORTED_FUNCTION void CosineSvmlSse(float* vO, float* vI, int numElements);
	//EXPORTED_FUNCTION void CosineSvmlAvx(float* vO, float* vI, int numElements);
	EXPORTED_FUNCTION void ExpSleefSse(float* vO, float* vI, int numElements);
	EXPORTED_FUNCTION void ExpSleefAvx(float* vO, float* vI, int numElements);
	EXPORTED_FUNCTION void ExpSvmlSse(float* vO, float* vI, int numElements);
	EXPORTED_FUNCTION void ExpSvmlAvx(float* vO, float* vI, int numElements);

#ifdef  __cplusplus
}
#endif

#endif



