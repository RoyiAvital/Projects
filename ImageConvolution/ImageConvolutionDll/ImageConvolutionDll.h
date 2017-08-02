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

	EXPORTED_FUNCTION void ApplyGammaCurveFa(float* mO, float* mI, int numRows, int numCols, int numColsPad, float gammaFactor);
	EXPORTED_FUNCTION void ApplyTonalMaskFa(float* mO, float* mI, float* mE, int numRows, int numCols, int numColsPad, float shadowsFactor, float midtonesFactor, float highlightsFactor);
	EXPORTED_FUNCTION void ApplyUnsharpMaskFa(float* mO, float* mI, float* mBlurredImage, int numRows, int numCols, int numColsPad, float ampLevel);
	EXPORTED_FUNCTION void ApplyUnsharpMaskPhotoshopFa(float* mO, float* mI, float* mBlurredImage, int numRows, int numCols, int numColsPad, float ampLevel, float thrLevel);
	EXPORTED_FUNCTION void BlendModeLinearBurnFa(float* mO, float* mF, float* mB, int numRows, int numCols, int numColsPad);
	EXPORTED_FUNCTION void BlendModeLinearDodgeFa(float* mO, float* mF, float* mB, int numRows, int numCols, int numColsPad);
	EXPORTED_FUNCTION void BlendModeOpacityFa(float* mO, float* mF, float* mB, int numRows, int numCols, int numColsPad, float opacityFctr);
	EXPORTED_FUNCTION void BlendModeOverlayFa(float* mO, float* mF, float* mB, int numRows, int numCols, int numColsPad);
	EXPORTED_FUNCTION void ConvertColorToGrayFa(float *mO, float *mR, float *mG, float *mB, int numRows, int numCols, int numColsPad, float *colCoeff);
	EXPORTED_FUNCTION void ExtractLChannelFa(float* mL, float* mR, float* mG, float* mB, int numRows, int numCols, int numColsPad, int sRgbMode, int whitePoint);
	EXPORTED_FUNCTION void ExtractYChannelFa(float* mY, float* mR, float* mG, float* mB, int numRows, int numCols, int numColsPad);
	EXPORTED_FUNCTION void HighPassFilterPsAuxFa(float* mO, float* mI, float* mB, int numRows, int numCols, int numColsPad);
	EXPORTED_FUNCTION void SetYChannelOutputFa(float* mR, float* mG, float* mB, float* mYIn, float* mYOut, int numRows, int numCols, int numColsPad, float opacityLevel);

#ifdef  __cplusplus
}
#endif

#endif



