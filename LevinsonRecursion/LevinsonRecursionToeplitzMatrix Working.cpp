#include "LevinsonRecursionToeplitzMatrix.h"

// --------------------------------------- CalcDistMatSse -------------------------------------- //
void LevinsonRecursionToeplitzMatrix(float *mT, float *vY, float *vX, int numRows)
{
	int ii, jj;
	float epsF, epsB, epsX, scalingFctr;
	//DECLARE_ALIGN float coordDiff[SSE_STRIDE];
	float *vF, *vFPrev, *vB, *vBPrev;

	vF		= (float*)_aligned_malloc(numRows * sizeof(float), SSE_ALIGNMENT);
	vFPrev	= (float*)_aligned_malloc(numRows * sizeof(float), SSE_ALIGNMENT);
	vB		= (float*)_aligned_malloc(numRows * sizeof(float), SSE_ALIGNMENT);
	vBPrev	= (float*)_aligned_malloc(numRows * sizeof(float), SSE_ALIGNMENT);
	//vX = (float*)_aligned_malloc(numRows * sizeof(float), SSE_ALIGNMENT);

	for (ii = 0; ii < numRows; ii++){
		vF[ii] = 0;
		vB[ii] = 0;
		vX[ii] = 0;
	}

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
		
		for (jj = 0; jj <= ii; jj++){
			vFPrev[jj] = vF[jj];
			vBPrev[jj] = vB[jj];
		}

		scalingFctr = 1 / (1 - epsF * epsB);

		for (jj = 0; jj <= ii; jj++){
			if (jj == 0){
				vF[0] = scalingFctr * vF[0];
			}
			else{
				vF[jj] = (scalingFctr * vF[jj]) - ((epsF * scalingFctr) * vB[jj - 1]);
			}
		}
		for (jj = 0; jj <= ii; jj++){
			if (jj == 0){
				vB[0] = -(epsB * scalingFctr) * vFPrev[0];
			}
			else{
				vB[jj] = (scalingFctr * vBPrev[jj - 1]) - ((epsB * scalingFctr) * vFPrev[jj]);
			}
		}
		for (jj = 0; jj <= ii; jj++){
			vX[jj] = vX[jj] + ((vY[ii] - epsX) * vB[jj]);
		}
	}

	_aligned_free(vF);
	_aligned_free(vB);
}
