#include "LevinsonRecursionToeplitzMatrix.h"
#include "xmmintrin.h"

// --------------------------------------- CalcDistMatSse -------------------------------------- //
void LevinsonRecursionToeplitzMatrix(float *mT, float *vY, float *vX, int numRows)
{
	int ii, jj;
	float epsF, epsB, epsX, scalingFctr;
	//DECLARE_ALIGN float coordDiff[SSE_STRIDE];
	float *vF, *vFPrev, *vB, *vBPrev;

	vF		= (float*)_mm_malloc(numRows * sizeof(float), 16);
	vFPrev	= (float*)_mm_malloc(numRows * sizeof(float), 16);
	vB		= (float*)_mm_malloc(numRows * sizeof(float), 16);
	vBPrev	= (float*)_mm_malloc(numRows * sizeof(float), 16);
	//vX = (float*)_aligned_malloc(numRows * sizeof(float), SSE_ALIGNMENT);


	// Use faster initialization
	/*for (ii = 0; ii < numRows; ii++){
		vF[ii] = 0;
		vB[ii] = 0;
		vX[ii] = 0;
	}*/

	memset(vF, 0.0, numRows * sizeof(float));
	memset(vB, 0.0, numRows * sizeof(float));
	memset(vX, 0.0, numRows * sizeof(float));

	vF[0] = 1 / mT[0];
	vB[0] = 1 / mT[0];
	vX[0] = vY[0] / mT[0];

	for (ii = 1; ii < numRows; ii++){
		epsF = 0;
		epsB = 0;
		epsX = 0;
		for (jj = 0; jj < ii; jj++){
			epsF += mT[(jj * numRows) + ii] * vF[jj];
			epsB += mT[(jj + 1) * numRows] * vB[jj];
			epsX += mT[(jj * numRows) + ii] * vX[jj];
		}
		
		// Use faster MemCopy
		/*for (jj = 0; jj <= ii; jj++){
			vFPrev[jj] = vF[jj];
			vBPrev[jj] = vB[jj];
		}*/

		memcpy(vFPrev, vF, (ii + 1) * sizeof(float));
		memcpy(vBPrev, vB, (ii + 1) * sizeof(float));

		scalingFctr = 1 / (1 - epsF * epsB);

		vF[0] = scalingFctr * vF[0];
		for (jj = 1; jj <= ii; jj++){
			vF[jj] = (scalingFctr * vF[jj]) - ((epsF * scalingFctr) * vB[jj - 1]);

		}

		vB[0] = -(epsB * scalingFctr) * vFPrev[0];
		for (jj = 1; jj <= ii; jj++){
			vB[jj] = (scalingFctr * vBPrev[jj - 1]) - ((epsB * scalingFctr) * vFPrev[jj]);

		}

		for (jj = 0; jj <= ii; jj++){
			vX[jj] = vX[jj] + ((vY[ii] - epsX) * vB[jj]);
		}
	}

	_mm_free(vF);
	_mm_free(vB);
}
