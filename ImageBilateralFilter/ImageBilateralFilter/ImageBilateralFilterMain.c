#include "ImageBilateralFilterMain.h"
#include "ImageBilateralFilter.h"

int main(int argc, char *argv[]) {
	int numRows = 4000;
	int numCols = 4000;
	int numIter = 10;

	int paramK;
	float spatialStd, rangeStd;

	spatialStd	= 3.0f;
	rangeStd	= 10.0f / 255.0f;
	paramK		= 5;

	int unitTestFlag;

	if (argc == 1) {
		unitTestFlag = UNIT_TEST_IMAGE_BILATERAL_FILTER;
	}
	else if (argc > 8) {
		printf("Error: Number of Argument Is Larger Than 7");
		return 0;
	}
	else { // argc == 2 -> Single Input
		unitTestFlag	= atoi(argv[1]);
		numRows			= atoi(argv[2]);
		numCols			= atoi(argv[3]);
		numIter			= atoi(argv[4]);

		switch (unitTestFlag)
		{
		case UNIT_TEST_IMAGE_GAUSSIAN_BLUR:
			spatialStd	= atof(argv[6]);
			paramK		= atoi(argv[7]);
            break;
		case UNIT_TEST_IMAGE_BILATERAL_FILTER:
			spatialStd	= atof(argv[5]);
			rangeStd	= atof(argv[6]);
			paramK		= atoi(argv[7]);
			break;
		default:
			break;
		}
	}

	switch (unitTestFlag)
	{
	case UNIT_TEST_IMAGE_GAUSSIAN_BLUR:
		ImageConvolutionGaussianKernelUnitTest(numRows, numCols, numIter, spatialStd, paramK);
		break;
	case UNIT_TEST_IMAGE_BILATERAL_FILTER:		
		BilateralFilterFastCompressiveUnitTest(numRows, numCols, numIter, spatialStd, rangeStd, paramK);
		break;
	default:
		break;
	}


	return 0;
}

