#ifdef _USRDLL
#define EXPORT_FCNS
#include "../CalcDistanceMatrix/CalcDistanceMatrixDll.h"
#endif // _USRDLL
#include "CalcDistanceMatrix.h"
#ifndef _USRDLL
// #include "mersenneTwister2002.c" // C File, can't be included in the .h file
#endif // !_USRDLL

// ------------------------------------ CalcDistanceMatrixSse ----------------------------------- //
/*
Calculates the Distance Matrix between 2 sets of vectors. The output mD(i, j) = dist(mA(:, i), mB(:, j)).
This is an SSE Optimized version.
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
TODO :
1.	C
Release Notes:
-	1.0.000	20/04/2018	Royi Avital
*   First release version.
*/
// ------------------------------------ CalcDistanceMatrixSse ----------------------------------- //
void CalcDistanceMatrixSse(float* mD, float* mA, float* mB, int vecDim, int numRowsA, int numRowsB)
{

	int ii, jj, kk, vecDimSseFloatPack;
	float distVal;
	float* ptrMatrixA;
	float* ptrMatrixB;
	float* ptrMatrixD;

	__m128 diffValSse;
	__m128 distValSse;

	vecDimSseFloatPack = vecDim - (vecDim % SSE_STRIDE);

#pragma omp parallel for private(ptrMatrixB, ptrMatrixD, jj, ptrMatrixA, distValSse, distVal, kk, diffValSse)
	for (ii = 0; ii < numRowsB; ii++) {
		ptrMatrixB = &mB[ii * vecDim];
		ptrMatrixD = &mD[ii * numRowsA];
		for (jj = 0; jj < numRowsA; jj++)
		{
			ptrMatrixA = &mA[jj * vecDim];
			distValSse = _mm_setzero_ps();
			distVal = 0.0f;
			for (kk = 0; kk < vecDimSseFloatPack; kk += SSE_STRIDE)
			{
				diffValSse = _mm_sub_ps(_mm_loadu_ps(&ptrMatrixA[kk]), _mm_loadu_ps(&ptrMatrixB[kk]));
				distValSse = _mm_add_ps(distValSse, _mm_mul_ps(diffValSse, diffValSse));
			}
			distVal = HorizontalSumSse(distValSse);
			for (kk = vecDimSseFloatPack; kk < vecDim; kk++)
			{
				distVal += (ptrMatrixA[kk] - ptrMatrixB[kk]) * (ptrMatrixA[kk] - ptrMatrixB[kk]);
			}
			ptrMatrixD[jj] = distVal;
		}
	}


}


// ------------------------------------ CalcDistanceMatrixAvx ----------------------------------- //
/*
Calculates the Distance Matrix between 2 sets of vectors. The output mD(i, j) = dist(mA(:, i), mB(:, j)).
This is an AVX Optimized version.
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
TODO :
1.	C
Release Notes:
-	1.0.000	20/04/2018	Royi Avital
*   First release version.
*/
// ------------------------------------ CalcDistanceMatrixAvx ----------------------------------- //
void CalcDistanceMatrixAvx(float* mD, float* mA, float* mB, int vecDim, int numRowsA, int numRowsB)

{

	int ii, jj, kk, vecDimAvxFloatPack;
	float distVal;
	float* ptrMatrixA;
	float* ptrMatrixB;
	float* ptrMatrixD;

	__m256 diffValAvx;
	__m256 distValAvx;
	__m128 lowPack, highPack, shufReg, sumsReg;

	vecDimAvxFloatPack = vecDim - (vecDim % AVX_STRIDE);

#pragma omp parallel for private(ptrMatrixB, ptrMatrixD, jj, ptrMatrixA, distValAvx, distVal, kk, diffValAvx, lowPack, highPack, shufReg, sumsReg)
	for (ii = 0; ii < numRowsB; ii++) {
		ptrMatrixB = &mB[ii * vecDim];
		ptrMatrixD = &mD[ii * numRowsA];
		for (jj = 0; jj < numRowsA; jj++)
		{
			ptrMatrixA = &mA[jj * vecDim];
			distValAvx = _mm256_setzero_ps();
			distVal = 0.0f;
			for (kk = 0; kk < vecDimAvxFloatPack; kk += AVX_STRIDE)
			{
				diffValAvx = _mm256_sub_ps(_mm256_loadu_ps(&ptrMatrixA[kk]), _mm256_loadu_ps(&ptrMatrixB[kk]));
				distValAvx = _mm256_add_ps(distValAvx, _mm256_mul_ps(diffValAvx, diffValAvx));
			}

			// Works only with GCC (Only when writing functions explicilty)
			lowPack = _mm256_castps256_ps128(distValAvx);
			highPack = _mm256_extractf128_ps(distValAvx, 1); // high 128
			lowPack = _mm_add_ps(lowPack, highPack);     // add the low 128

			shufReg = _mm_movehdup_ps(lowPack);        // Broadcast elements 3,1 to 2,0
			sumsReg = _mm_add_ps(lowPack, shufReg);
			shufReg = _mm_movehl_ps(shufReg, sumsReg); // High Half -> Low Half
			sumsReg = _mm_add_ss(sumsReg, shufReg);
			distVal = _mm_cvtss_f32(sumsReg); // Result in the lower part of the SSE Register


			// distVal = HorizontalSumAvx(distValAvx); // For some reason doesn't work as EXE (Only as DLL in MATLAB)
			for (kk = vecDimAvxFloatPack; kk < vecDim; kk++)
			{
				distVal += (ptrMatrixA[kk] - ptrMatrixB[kk]) * (ptrMatrixA[kk] - ptrMatrixB[kk]);
			}
			ptrMatrixD[jj] = distVal;
		}
	}


}

