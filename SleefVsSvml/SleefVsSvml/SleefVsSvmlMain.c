#include "SleefVsSvmlMain.h"

int main(int argc, char *argv[]) {
	// See https://docs.microsoft.com/en-us/cpp/c-language/arguments-to-main
	int numElements, numIter, negValFlag;
	int unitTestFlag;
	float maxVal;

	if (argc == 1) {
		 printf("\nError: Undefined Unit Test\n");
		 return 0;
	}
	else if (argc == 2) {
		unitTestFlag	= atoi(argv[1]);
		numElements		= 1e7;
		numIter			= 100;
		negValFlag		= OFF;
		maxVal			= 100;
		printf("\nUsing Default Settings\n");
	}
	else if (argc > 6) {
		printf("\nError: Number of Argument Is Larger Than 5\n");
		return 0;
	}
	else { // argc ==2 -> Single Input
		unitTestFlag	= atoi(argv[1]);
		numElements		= atoi(argv[2]);
		numIter			= atoi(argv[3]);
		negValFlag		= atoi(argv[4]);
		maxVal			= (float)(atof(argv[5]));
        printf("\nUsing Users Settings\n");
	}
    
    printf("Number of Elements - %d\n", numElements);
	printf("Number of Iterations - %d\n", numIter);
	printf("Negative Values Flag - %d\n", negValFlag);
	printf("Maximum Value - %f\n", maxVal);

	switch (unitTestFlag)
	{
	case UNIT_TEST_SINE:
		UnitTestSine(numElements, numIter, maxVal, negValFlag);
		break;
	case UNIT_TEST_EXP:
		UnitTestExp(numElements, numIter, maxVal, negValFlag);
		break;
	default:
		printf("\n");
		printf("Error: Invalid Test Number");
		printf("\n");
		break;
	}


	return 0;
}

