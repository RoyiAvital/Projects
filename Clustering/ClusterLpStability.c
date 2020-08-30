/*
 * Clustering via LP Based Stabilities (`ClusterLpStability()`)
 * This module implements the Clustering via LP Based Stabilities algorithm.
 * 
 *
 * References:
 *	1. 	A
 * Remarks:
 *	1.	Matrices are column wise contiguous.
 * TODO:
 *	1.	C
 * Release Notes:
 *	-	1.0.000	30/08/2020	Royi Avital
 *		*	First release version.
 */

#define ELEMENT_TYPE_IDX 1

#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>

#define ELM_TYPE double
#include "../AuxiliaryFunctions/AuxFun.c"

struct structLpStability
{
	unsigned int numSamples;
	double* mD;
	double* mH;
	unsigned int numMedoids;
	unsigned int* vQ;
	unsigned int* vQs; // Sorted version of vQ
	unsigned int numNotMedoids;
	unsigned int* vP;
	unsigned int* vBuffer1; // _Distribute()
	unsigned int* vBuffer2; // _Distribute()
	unsigned int* vBuffer3; // _Distribute()
	unsigned int* vBuffer4;	// _Distribute()
	unsigned int* vBuffer5; // _CalcMargin()
	unsigned int* vBuffer6; // _CalcMargin
	double* vBuffer7; // _Distribute()
	double* vBuffer8; // _Distribute()
	double* vBuffer9;
	double* mBuffer1;
	double* mBuffer2;
};

//void _SetDiffQ(unsigned int *vP, unsigned int* vQ, unsigned int numMedoids, unsigned int numSamples)
//{
//	unsigned int ii, idxQ, idxP;
//	if (numMedoids == 0)
//	{
//		for (ii = 0; ii < numSamples; ii++)
//		{
//			vP[ii] = ii;
//		}
//		return;
//	}
//
//	// Better assume vQ is sorted
//	InsertionSort(vQ, numMedoids); // Optimized for paritally sorted array
//
//	for (ii = 0; ii < numSamples; ii++)
//	{
//		if (vQ[idxQ] == ii)
//		{
//			idxQ++;
//			idxQ = min(idxQ, numMedoids);
//		}
//		else
//		{
//			vP[idxP] = ii;
//			idxP++;
//		}
//	}
//}

int _CmpDouble(const void* a, const void* b)
{
	if (*(double*)a > * (double*)b)
	{
		return 1;
	}
	else if (*(double*)a < *(double*)b)
	{
		return -1;
	}
	else
	{
		return 0;
	}
}

void _SortArray(double* vA, double* vB, unsigned int numElements)
{
	memcpy(vA, vB, numElements * sizeof(double));
	qsort(vA, numElements, sizeof(double), _CmpDouble);
}

void _UpdateP(unsigned int* vP, unsigned int numNotMedoids, unsigned int qq)
{
	unsigned int ii;
	int shiftFlag;

	shiftFlag = 0;
	
	for (ii = 0; ii < numNotMedoids; ii++)
	{
		shiftFlag = shiftFlag || (vP[ii] == qq);
		if (shiftFlag)
		{
			vP[ii] = vP[ii + 1];
		}
	}
}

int _IsMember(unsigned int valIn, unsigned int* vA, unsigned int numElements)
{
	unsigned int ii;
	int isMember;

	isMember = FALSE;

	for (ii = 0; ii < numElements; ii++)
	{
		isMember = (valIn == vA[ii]);
		if (isMember)
		{
			break;
		}
	}

	return isMember;
}

double _CalcMargin( struct structLpStability* sLpStability, unsigned int qq )
{
	double valMargin, *mH, *mD, minVal, minValQ;
	unsigned int numSamples, * vP, pp, numNotMedoids, ii, jj, idxH, idxMin, idxMinQ, idxPq, idxPn, idxPnq, idxQq, idxQn, *vN, *vNq;

	valMargin = 0.0;

	numSamples = sLpStability->numSamples;
	mH = sLpStability->mH;
	mD = sLpStability->mD;
	vP = sLpStability->vP;
	numNotMedoids = sLpStability->numNotMedoids;
	vN = sLpStability->vBuffer5; // Prevent interaction with _Distribute()
	vNq = sLpStability->vBuffer6; // Prevent interaction with _Distribute()

	for (ii = 0; ii < numSamples; ii++)
	{
		idxH = ii; // Iterating over the row
		minVal = mH[ii];
		idxMin = 0;
		
		if (qq != 0)
		{
			minValQ = mH[ii]; // Without the qq column
			idxMinQ = 0; // Without the qq column
		}
		else
		{
			minValQ = mH[ii + numSamples]; // Without the qq column
			idxMinQ = 1; // Without the qq column
		}

		for (jj = 1; jj < numSamples; jj++)
		{
			idxH += numSamples;
			if (mH[idxH] < minVal)
			{
				minVal = mH[idxH];
				idxMin = jj;
			}
			if ((mH[idxH] < minValQ) && (jj != qq))
			{
				minValQ = mH[idxH];
				idxMinQ = jj;
			}
		}
		vN[ii] = idxMin;
		vNq[ii] = idxMinQ;
	}

	for (ii = 0; ii < numNotMedoids; ii++)
	{
		pp = vP[ii];

		idxPq	= (qq * numSamples) + pp;
		idxPn	= (vN[pp] * numSamples) + pp;
		idxPnq	= (vNq[pp] * numSamples) + pp;

		if (mH[idxPq] == mH[idxPn])
		{
			valMargin += mH[idxPnq] - mH[idxPq];
		}
		if (pp != qq)
		{
			valMargin -= (mH[idxPq] - max(mH[idxPn], mD[idxPq]));
		}
	}

	idxQq = (qq * numSamples) + qq;
	idxQn = (vN[qq] * numSamples) + qq;
	
	valMargin -= (mH[idxQq] - mH[idxQn]);

	return valMargin;

}

double _CalcDualObjective( struct structLpStability* sLpStability )
{
	double daulObjVal, * mH, minVal;
	unsigned int numSamples, ii, jj, idxH;

	numSamples = sLpStability->numSamples;
	mH = sLpStability->mH;

	daulObjVal = 0.0;

	for (ii = 0; ii < numSamples; ii++)
	{
		idxH = ii;
		minVal = mH[ii];
		for (jj = 1; jj < numSamples; jj++)
		{
			idxH += numSamples;
			minVal = min(minVal, mH[idxH]);
		}
		daulObjVal += minVal;
	}

	return daulObjVal;

}

void _Distribute(struct structLpStability* sLpStability)
{
	double *mH, *mD, minVal, *mHH, *mHT, *vNVal, * vH, valMarginQ;
	unsigned int numSamples, numElements, numMedoids, Nqc, qq, ii, jj, idxH, idxMin, *vQ, *vQc, *vN, *vLq, *vV, *vIdx, pp, idxPn, idxQ, idxD, cardVq, shiftVal, idxPq;

	numSamples = sLpStability->numSamples;
	numMedoids = sLpStability->numMedoids;
	mD = sLpStability->mD;
	mH = sLpStability->mH;
	vQ = sLpStability->vQ;
	vQc = sLpStability->vP;
	Nqc = sLpStability->numNotMedoids;
	vN = sLpStability->vBuffer1;
	vNVal = sLpStability->vBuffer7;
	vH = sLpStability->vBuffer8;
	vLq = sLpStability->vBuffer2;
	vV = sLpStability->vBuffer4;
	vIdx = sLpStability->vBuffer3;
	mHH = sLpStability->mBuffer1;
	mHT = sLpStability->mBuffer2;

	numElements = numSamples * numSamples;

	memcpy(mHH, mH, numElements * sizeof(double));
	memcpy(mHT, mH, numElements * sizeof(double));

	for (ii = 0; ii < numSamples; ii++)
	{
		idxH = ii;
		minVal = mH[ii];
		idxMin = 0;

		for (jj = 1; jj < numSamples; jj++)
		{
			idxH += numSamples;
			if (mH[idxH] < minVal)
			{
				minVal = mH[idxH];
				idxMin = jj;
			}
		}
		vN[ii] = idxMin;
		vNVal[ii] = minVal;
	}

	for (ii = 0; ii < Nqc; ii++)
	{
		pp = vQc[ii];
		idxPn = (vN[pp] * numSamples) + pp;
		mHH[idxPn] = 1e30;
	}

	for (ii = 0; ii < numSamples; ii++)
	{
		idxH = ii;
		minVal = mHH[ii];
		for (jj = 1; jj < numSamples; jj++)
		{
			idxH += numSamples;
			minVal = min(minVal, mHH[idxH]);
		}
		vH[ii] = minVal; // h_hat in MATLAB

		vLq[ii] = numSamples + 1;
	}

	idxQ = 0;
	for (ii = 0; ii < Nqc; ii++)
	{
		if (_IsMember(vN[vQc[ii]], vQ, numMedoids))
		{
			vLq[idxQ] = vQc[ii];
			idxQ++;
		}
		vIdx[ii] = FALSE;
	}

	for (ii = 0; ii < Nqc; ii++)
	{
		vV[ii] = !(_IsMember(vQc[ii], vLq, numMedoids));
	}

	for (ii = 0; ii < Nqc; ii++)
	{
		qq = vQc[ii];
		valMarginQ = _CalcMargin(sLpStability, qq);

		cardVq = 0;
		shiftVal = 1;
		for (jj = 0; jj < Nqc; jj++)
		{
			idxD = (qq * numSamples) + vQc[jj];
			vIdx[jj] = vV[jj] && (vNVal[vQc[jj]] >= mD[idxD]);
			if (vIdx[jj])
			{
				cardVq++;
				if (vQc[jj] == qq)
				{
					shiftVal = 0;
				}
			}
			
		}
		cardVq += shiftVal;
		for (jj = 0; jj < Nqc; jj++)
		{
			pp = vQc[jj];
			idxPq = (qq * numSamples) + pp;
			if ((pp != qq) && (_IsMember(pp, vLq, Nqc) || vNVal[pp] < mD[idxPq]))
			{
				mHT[idxPq] = max(vNVal[pp], mD[idxPq]);
			}
			else if (mHT[idxPq] > vNVal[pp])
			{
				mHT[idxPq] = vNVal[pp] - (valMarginQ / (double)cardVq);
			}
			else if (mH[idxPq] == vNVal[pp])
			{
				mHT[idxPq] = vH[pp] - (valMarginQ / (double)cardVq);
			}
		}
	}

	memcpy(mH, mHT, numElements * sizeof(double));

}

void _ProjectMedoid( struct structLpStability * sLpStability, unsigned int qq )
{
	double *mH, *mD;
	unsigned int numNotMedoids, numSamples, ii, pp, *vP, idxPp, idxPq, idxQp, idxQq;

	numNotMedoids = sLpStability->numNotMedoids;
	vP = sLpStability->vP;
	mH = sLpStability->mH;
	mD = sLpStability->mD;
	numSamples = sLpStability->numSamples;

	// _SetDiffQ(vP, vQs, numMedoids, numSamples);

	for (ii = 0; ii < numNotMedoids; ii++)
	{
		pp = vP[ii];
		idxPp = (pp * numSamples) + pp;
		idxPq = (qq * numSamples) + pp;
		idxQp = (pp * numSamples) + qq;

		mH[idxPp] = mH[idxPp] + mH[idxQp] - mD[idxQp];
		mH[idxQp] = mD[idxQp];
		mH[idxPq] = mD[idxPq]; // Should be symmetric so might not be needed for the new idx

	}
	
	idxQq = (qq * numSamples) + qq;
	mH[idxQq] = mD[idxQq];
}

unsigned int _SearchStablePoint( struct structLpStability* sLpStability )
{
	double valEps, dualObjVal, dualObjValPrev, maxMargin, currMargin;
	unsigned int *vQc, Nqc, numSamples, ii, idxMax;

	valEps = 1e-5;

	numSamples	= sLpStability->numSamples;
	vQc			= sLpStability->vP;
	Nqc			= sLpStability->numNotMedoids;

	maxMargin	= _CalcMargin(sLpStability, vQc[0]);
	idxMax		= 0;

	for (ii = 1; ii < Nqc; ii++)
	{
		currMargin = _CalcMargin(sLpStability, vQc[ii]);
		if (currMargin > maxMargin)
		{
			maxMargin	= currMargin;
			idxMax		= ii;
		}
	}

	dualObjVal		= _CalcDualObjective(sLpStability);
	dualObjValPrev	= 1e20;

	while ((maxMargin < 0.0) && ((fabs(dualObjVal - dualObjValPrev) / numSamples) > valEps))
	{
		dualObjValPrev = dualObjVal;
		_Distribute(sLpStability);
		dualObjVal = _CalcDualObjective(sLpStability);
		
		maxMargin = _CalcMargin(sLpStability, vQc[0]);
		idxMax = 0;

		for (ii = 1; ii < Nqc; ii++)
		{
			currMargin = _CalcMargin(sLpStability, vQc[ii]);
			if (currMargin > maxMargin)
			{
				maxMargin = currMargin;
				idxMax = ii;
			}
		}
	}

	return vQc[idxMax];
	
}

void ClusterLpStability( unsigned int* vMedoidIdx, unsigned int numMedoids[1], double* mD, unsigned int numSamples, double paramMu, unsigned int maxNumMedoids )
{
	unsigned int numElements, ii, qq, *vP;
	double medoidPenalty;
	
	struct structLpStability sLpStability;
	
	numElements = numSamples * numSamples;

	medoidPenalty = paramMu * CalcMedianVal(mD, numElements);
	// mexPrintf("medoidPenalty = %f\n", medoidPenalty);
	maxNumMedoids = min(numSamples, maxNumMedoids);

	sLpStability.numSamples = numSamples;
	sLpStability.mD = (double*)_mm_malloc(numElements * sizeof(double), SIMD_ALIGNMENT);
	sLpStability.mH = (double*)_mm_malloc(numElements * sizeof(double), SIMD_ALIGNMENT);
	sLpStability.numMedoids = 0;
	sLpStability.vQ = (unsigned int*)_mm_malloc(maxNumMedoids * sizeof(unsigned int), SIMD_ALIGNMENT);
	sLpStability.vQs = (unsigned int*)_mm_malloc(maxNumMedoids * sizeof(unsigned int), SIMD_ALIGNMENT); // Sorted version of vQ
	sLpStability.numNotMedoids = numSamples;
	sLpStability.vP = (unsigned int*)_mm_malloc(numSamples * sizeof(unsigned int), SIMD_ALIGNMENT);
	sLpStability.vBuffer1 = (unsigned int*)_mm_malloc(numSamples * sizeof(unsigned int), SIMD_ALIGNMENT);
	sLpStability.vBuffer2 = (unsigned int*)_mm_malloc(numSamples * sizeof(unsigned int), SIMD_ALIGNMENT);
	sLpStability.vBuffer3 = (unsigned int*)_mm_malloc(numSamples * sizeof(unsigned int), SIMD_ALIGNMENT);
	sLpStability.vBuffer4 = (unsigned int*)_mm_malloc(numSamples * sizeof(unsigned int), SIMD_ALIGNMENT);
	sLpStability.vBuffer5 = (unsigned int*)_mm_malloc(numSamples * sizeof(unsigned int), SIMD_ALIGNMENT);
	sLpStability.vBuffer6 = (unsigned int*)_mm_malloc(numSamples * sizeof(unsigned int), SIMD_ALIGNMENT);
	sLpStability.vBuffer7 = (double*)_mm_malloc(numSamples * sizeof(double), SIMD_ALIGNMENT);
	sLpStability.vBuffer8 = (double*)_mm_malloc(numSamples * sizeof(double), SIMD_ALIGNMENT);
	sLpStability.vBuffer9 = (double*)_mm_malloc(numSamples * sizeof(double), SIMD_ALIGNMENT);
	sLpStability.mBuffer1 = (double*)_mm_malloc(numElements * sizeof(double), SIMD_ALIGNMENT);
	sLpStability.mBuffer2 = (double*)_mm_malloc(numElements * sizeof(double), SIMD_ALIGNMENT);

	memcpy(sLpStability.mD, mD, numElements * sizeof(double));

	//_SortArray(sLpStability.mBuffer1, mD, numElements);
	//medoidPenalty = paramMu * ((sLpStability.mBuffer1[numElements / 2] + sLpStability.mBuffer1[(numElements / 2) - 1]) / 2.0);
	//mexPrintf("medoidPenalty = %f\n", medoidPenalty);

	for (ii = 0; ii < numElements; ii += (numSamples + 1))
	{
		sLpStability.mD[ii] += medoidPenalty;
	}

	memcpy(sLpStability.mH, sLpStability.mD, numElements * sizeof(double));

	vP = sLpStability.vP;
	for (ii = 0; ii < numSamples; ii++)
	{
		vP[ii] = ii;
	}

	qq = _SearchStablePoint(&sLpStability); // Index of medoid candidate

	while (_CalcMargin(&sLpStability, qq) >= 0)
	{
		sLpStability.vQ[sLpStability.numMedoids] = qq;
		sLpStability.vQs[sLpStability.numMedoids] = qq;
		sLpStability.numMedoids++;
		if (sLpStability.numMedoids == maxNumMedoids)
		{
			break;
		}
		sLpStability.numNotMedoids--;
		_UpdateP(sLpStability.vP, sLpStability.numNotMedoids, qq);
		_ProjectMedoid(&sLpStability, qq);
		qq = _SearchStablePoint(&sLpStability);
	}

	numMedoids[0] = sLpStability.numMedoids;
	for (ii = 0; ii < numMedoids[0]; ii++)
	{
		vMedoidIdx[ii] = sLpStability.vQ[ii];
	}

	_mm_free(sLpStability.mD);
	_mm_free(sLpStability.mH);
	_mm_free(sLpStability.vQ);
	_mm_free(sLpStability.vQs);
	_mm_free(sLpStability.vP);
	_mm_free(sLpStability.vBuffer1);
	_mm_free(sLpStability.vBuffer2);
	_mm_free(sLpStability.vBuffer3);
	_mm_free(sLpStability.vBuffer4);
	_mm_free(sLpStability.vBuffer5);
	_mm_free(sLpStability.vBuffer6);
	_mm_free(sLpStability.mBuffer1);

}