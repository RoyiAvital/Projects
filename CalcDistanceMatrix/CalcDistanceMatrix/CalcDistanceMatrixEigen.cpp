#ifdef _USRDLL
#define EXPORT_FCNS
#include "../CalcDistanceMatrix/CalcDistanceMatrixDll.h"
#endif // _USRDLL
#include "CalcDistanceMatrix.h"
#ifndef _USRDLL
// #include "mersenneTwister2002.c" // C File, can't be included in the .h file
#endif // !_USRDLL

// ----------------------------------- CalcDistanceMatrixEigen ---------------------------------- //
/*
Calculates the Distance Matrix between 2 sets of vectors. The output mD(i, j) = dist(mA(:, i), mB(:, j)).
This is an Eigen based version.
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
	array the vectors are along a row. Yet Eigen uses Column Major definition of Matrix.
2.	This version works between a vector and a matrix.
TODO :
1.	C
Release Notes:
-	1.0.000	20/04/2018	Royi Avital
*   First release version.
*/
// ----------------------------------- CalcDistanceMatrixEigen ---------------------------------- //
void CalcDistanceMatrixEigen(float* mD, float* mA, float* mB, int vecDim, int numRowsA, int numRowsB)
{

	EigenMatExt meA(mA, vecDim, numRowsA);
	EigenMatExt meB(mB, vecDim, numRowsB);

	EigenMatExt meD(mD, numRowsA, numRowsB);

	int rowIdx;

	for (rowIdx = 0; rowIdx < numRowsA; rowIdx++)
	{
		meD.row(rowIdx) = (meB.colwise() - meA.col(rowIdx)).matrix().colwise().squaredNorm(); // Broadcasting -> The '-' must be after the .colwise()
	}


}


// ----------------------------------- CalcDistanceMatrixEigen ---------------------------------- //
/*
Calculates the Distance Matrix between 2 sets of vectors. The output mD(i, j) = dist(mA(:, i), mB(:, j)).
This is an Eigen based version.
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
	array the vectors are along a row. Yet Eigen uses Column Major definition of Matrix.
 2.	This version uses the fact that in MATLAB Notation:
	mD = sum(mA .^ 2, 1).' - (2 * mA.' * mB) + sum(mB .^ 2, 1)
	Yet in practice it is not as fast (Maybe since MATLAB uses Intel MKL).
TODO :
 1.	C
Release Notes:
 -	1.0.000	20/04/2018	Royi Avital
	*   First release version.
*/
// ----------------------------------- CalcDistanceMatrixEigen ---------------------------------- //
void CalcDistanceMatrixEigenM(float* mD, float* mA, float* mB, int vecDim, int numRowsA, int numRowsB)
{

	EigenMatExt meA(mA, vecDim, numRowsA);
	EigenMatExt meB(mB, vecDim, numRowsB);

	EigenMatExt meD(mD, numRowsA, numRowsB);

	// meD = meA.colwise().squaredNorm().transpose() - (2 * meA.transpose() * meB) + meB.colwise().squaredNorm();
	// Not even close to be as fast as MATLAB (Intel MKL???)
	meD = ((-2 * meA.transpose() * meB).colwise() + meA.colwise().squaredNorm().transpose()).rowwise() + meB.colwise().squaredNorm();


}

