#include "CalcDistanceMatrix.h"

#ifdef _USRDLL
#define EXPORT_FCNS
	#include "../CalcDistanceMatrix/CalcDistanceMatrixDll.h"
#endif // _USRDLL
#ifndef _USRDLL
// #include "mersenneTwister2002.c" // C File, can't be included in the .h file
#endif // !_USRDLL

// ---------------------------------- CalcDistanceMatrixVanilla --------------------------------- //
/*
Calculates the Distance Matrix between 2 sets of vectors. The output mD(i, j) = dist(mA(:, i), mB(:, j)).  
Input:
- mD			-	Distance Matrix.
					Structure: Matrix (numVectorsA x numVectorsB).
					Type : 'Single'.
					Range : (-inf, inf).
- mA			-	Vector Set A Matrix .
					Structure: Matrix (vecDim x numVectorsA).
					Type : 'Single'.
					Range : (-inf, inf).
- mB			-	Vector Set B Matrix .
					Structure: Matrix (vecDim x numVectorsB).
					Type : 'Single'.
					Range : (-inf, inf).
- vecDim		-	Vector Dimension.
					The number of elements in each vector (Column of a matrix).
					Structure: Scalar.
					Type : 'Int'.
					Range : {1, 2, ...}.
- numRowsA		-	Number of Rows in Matrix A.
					Basically number of Vectors in Set A.
					Since C is Row Major this is the number of rows as Vectors are contiguous along
					a Row.
					Structure: Scalar.
					Type : 'Int'.
					Range : {1, 2, ...}.
- numRowsB		-	Number of Rows in Matrix B.
					Basically number of Vectors in Set B.
					Since C is Row Major this is the number of rows as Vectors are contiguous along
					a Row.
					Structure: Scalar.
					Type : 'Int'.
					Range : {1, 2, ...}.
Reference:
 1.	See MATLAB's reference implementation.
Remarks:
 1.	The columns of a Matrix are contiguous in memory. Since C is Row Major it means that for 2D
	array the vectors are along a row.
 2.	This version should be written in a Compiler Friendly (Multi Threading & Vectorization) manner.
TODO :
 1.	C
Release Notes:
 -	1.0.000	20/04/2018	Royi Avital
	*   First release version.
*/
// ---------------------------------- CalcDistanceMatrixVanilla --------------------------------- //
#if defined(__GNUC__)
void CalcDistanceMatrixVanilla(float* mD, float* __restrict__ mA, float* __restrict__ mB, int vecDim, int numRowsA, int numRowsB)
#elif defined(__INTEL_COMPILER)
void CalcDistanceMatrixVanilla(float* mD, float* mA, float* mB, int vecDim, int numRowsA, int numRowsB)
#else
void CalcDistanceMatrixVanilla(float* mD, float* mA, float* mB, int vecDim, int numRowsA, int numRowsB)
#endif
{

	int ii, jj, kk;
	float distVal;
	float* ptrMatrixA;
	float* ptrMatrixB;
	float* ptrMatrixD;

// #if defined(__GNUC__) || defined(__INTEL_COMPILER) || !defined(_USRDLL)
#pragma omp parallel for private(ptrMatrixB, ptrMatrixD, jj, ptrMatrixA, distVal, kk)
// #endif
	for (ii = 0; ii < numRowsB; ii++) {
		ptrMatrixB = &mB[ii * vecDim];
		ptrMatrixD = &mD[ii * numRowsA];
		for (jj = 0; jj < numRowsA; jj++)
		{
			ptrMatrixA = &mA[jj * vecDim];
			distVal = 0.0f;
			for (kk = 0; kk < vecDim; kk++)
			{
				distVal += (ptrMatrixA[kk] - ptrMatrixB[kk]) * (ptrMatrixA[kk] - ptrMatrixB[kk]);
			}
			ptrMatrixD[jj] = distVal;
		}
	}


}


// Code to measure MATLAB in & Out Time
void CalcDistanceMatrixRefTime(float* mD, float* mA, float* mB, int vecDim, int numRowsA, int numRowsB)
{

	mD[1] = mA[1] + mB[1];


}

