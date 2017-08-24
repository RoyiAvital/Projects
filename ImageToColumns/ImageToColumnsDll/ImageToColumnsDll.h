#ifndef DLL_TEMPLATE
	#define DLL_TEMPLATE

	#define DATA_TYPE float
	// #define DATA_TYPE double

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

	EXPORTED_FUNCTION void ImageToColumns(DATA_TYPE* mO, DATA_TYPE* mI, int numRows, int numCols, int blockRadius);

	#ifdef  __cplusplus
		}
	#endif

#endif



