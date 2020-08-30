/*
 * Clustering via LP Based Stabilities (`ClusterLpStability()`) - MEX Wrapper
 * This module implements the Cholesky Decomposition. It uses upper boundary to the number of
 * elements to prevent recurring allocations. It is implemented with the preconditioned conjugate
 * gradient in mind hence the preconditioning step is implemented as well. The decomposed Matrix
 * is given by mA + shiftVal * diag(mA).
 * Input:
 *	- mD                -   Input Positive Definite Sparse Matrix.
 *							The matrix to decompose.
 *                          Structure: Sparse (CSC) Matrix (numRows * numCols).
 *                          Type: 'Double'.
 *                          Range: (-inf, inf).
 *	- paramMu			-   Discard Threshold.
 *							Values in the cholesky decomposition which are smaller than
 *							`discardThr` are zeroed. If `discardThr` = 0 then the full
 *							Cholesky Decomposition is applied (No discarding).
 *                          Structure: Scalar.
 *                          Type: 'Double'.
 *                          Range: [0, inf).
 * - maxNumMedoids      -   Discard Threshold.
 *							Values in the cholesky decomposition which are smaller than
 *							`discardThr` are zeroed.
 *                          Structure: Scalar.
 *                          Type: 'Double'.
 *                          Range: (0, inf).
 * Output:
 *	- vClusterIdx       -   Input Positive Definite Sparse Matrix.
 *							The matrix to decompose.
 *                          Structure: Sparse (CSC) Matrix (numRows * numCols).
 *                          Type: 'Double'.
 *                          Range: (-inf, inf).
 *  - vClusterIdx       -   Input Positive Definite Sparse Matrix.
 *							The matrix to decompose.
 *                          Structure: Sparse (CSC) Matrix (numRows * numCols).
 *                          Type: 'Double'.
 *                          Range: (-inf, inf).
 * References:
 *	1. 	https://github.com/pymatting/pymatting/blob/master/pymatting/preconditioner/ichol.py
 *	2.	Incomplete Cholesky Decomposition: https://en.wikipedia.org/wiki/Incomplete_Cholesky_factorization
 * Remarks:
 *	1.	The Sparse Matrices are given in Compressed Sparse Column (CSC) format.
 * TODO:
 *	1.	Add "Zero Fill" variant of the decomposition.
 * Release Notes:
 *	-	1.0.000	10/07/2020	Royi Avital
 *		*	First release version.
 */


#include <stdlib.h>

#include "mex.h"
#include "ClusterLpStability.c"

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

	mwSize numSamples;
	double *mD, paramMu, maxNumMedoids, * vMedoidIdxD;
	unsigned int ii, *vMedoidIdx, numMedoids[1];
	
	// Validate Input
	if ( nrhs != 3 )
	{
		mexErrMsgIdAndTxt("ClusterLpStability:ClusterLpStability:nrhs", "There must be 3 inputs: A real symmetric with non negative elements matrix of type double and 2 scalars of type double");
	}
	// Validating the 1st input
	if( !mxIsNumeric(prhs[0]) || mxIsEmpty(prhs[0]) || mxIsComplex(prhs[0]) || (mxGetClassID(prhs[0]) != mxDOUBLE_CLASS) || (mxGetM(prhs[0]) != mxGetN(prhs[0])) )
	{
		mexErrMsgIdAndTxt("ClusterLpStability:ClusterLpStability:nrhs", "Matrix mD: The 1st input must be a real symmetric with non negative elements matrix of type double");
	}

	// Validating the 2nd input
	if ( !mxIsNumeric(prhs[1]) || mxIsEmpty(prhs[1]) || mxIsComplex(prhs[1]) || (mxGetClassID(prhs[1]) != mxDOUBLE_CLASS) || (mxGetNumberOfElements(prhs[1]) != 1) )
	{
		mexErrMsgIdAndTxt("ClusterLpStability:ClusterLpStability:nrhs", "Parameter Mu: The 2nd input must be a real non negative scalar of type double");
	}

	// Validating the 3rd input
	if ( !mxIsNumeric(prhs[2]) || mxIsEmpty(prhs[2]) || mxIsComplex(prhs[2]) || (mxGetClassID(prhs[2]) != mxDOUBLE_CLASS) || (mxGetNumberOfElements(prhs[2]) != 1) )
	{
		mexErrMsgIdAndTxt("ClusterLpStability:ClusterLpStability:nrhs", "Maximum Number of Medoids: The 3rd input must be a real non negative integer scalar of type double");
	}

	// Validate Output
	if (nlhs != 1)
	{
		mexErrMsgIdAndTxt("ClusterLpStability:ClusterLpStability:nlhs", "There must be 2 output vairables");
	}
	
	// Input Matrix Dimensions -> Number of samples
	numSamples = mxGetM(prhs[0]);
	
	// Input Matrix Data
	mD = (double*)mxGetData(prhs[0]);

	paramMu			= (double)mxGetScalar(prhs[1]);
	maxNumMedoids	= (double)mxGetScalar(prhs[2]);

	if ( (paramMu < 0.0) || mxIsNaN(paramMu) || mxIsInf(paramMu) )
	{
		mexErrMsgIdAndTxt("ClusterLpStability:ClusterLpStability:nrhs", "Parameter Mu: The 2nd input must be a real non negative scalar of type double");
	}

	if ( (maxNumMedoids < 0.0) || (maxNumMedoids != round(maxNumMedoids)) || mxIsNaN(maxNumMedoids) || mxIsInf(maxNumMedoids) )
	{
		mexErrMsgIdAndTxt("ClusterLpStability:ClusterLpStability:nrhs", "Maximum Number of Medoids: The 3rd input must be a real non negative integer scalar of type double");
	}

	vMedoidIdx = (unsigned int*)_mm_malloc(numSamples * sizeof(unsigned int), SIMD_ALIGNMENT);

	ClusterLpStability(vMedoidIdx, numMedoids, mD, (unsigned int)numSamples, paramMu, (unsigned int)maxNumMedoids);

	plhs[0] = mxCreateNumericMatrix(numMedoids[0], 1, mxUINT32_CLASS, mxREAL);
	plhs[0] = mxCreateDoubleMatrix(numMedoids[0], 1, mxREAL);

	vMedoidIdxD = (double*)mxGetData(plhs[0]);

	for (ii = 0; ii < numMedoids[0]; ii++)
	{
		vMedoidIdxD[ii] = (double)vMedoidIdx[ii] + 1.0;
	}

	// mxSetUint32s(plhs[0], vMedoidIdx);

	_mm_free(vMedoidIdx);
	
	
}