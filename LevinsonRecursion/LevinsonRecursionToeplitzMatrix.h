#define SSE_ALIGNMENT 16

#include <stdio.h>
#include <string.h>


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

EXPORTED_FUNCTION void LevinsonRecursionToeplitzMatrix(float *mT, float *vY, float *vX, int numRows);

#ifdef  __cplusplus
}
#endif

#endif

