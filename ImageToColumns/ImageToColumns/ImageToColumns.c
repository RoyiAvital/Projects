#ifdef _USRDLL
	#define EXPORT_FCNS
	#include "../ImageToColumnsDll/ImageToColumnsDll.h"
#endif // _USRDLL

#ifndef _USRDLL
	#include "ImageToColumnsMain.h"
#endif // !_USRDLL

#include "ImageToColumns.h"
// -------------------------------------- ImageConvolution -------------------------------------- //
/*
Applies Convolution on an Image (2D Array) using given 2D Kernel.
Input:
- mO			-	Output Image.
					Structure: Image Array (Single Channel).
					Type : 'Single'.
					Range : [0, 1].
- mI			-	Input Image.
					Structure: Image Matrix (Single Channel).
					Type : 'Single'.
					Range : [0, 1].
- numRows		-	Number of Rows.
					Structure: Scalar.
					Type : 'Int'.
					Range : {1, 2, ...}.
- numCols		-	Number of Columns.
					Assumbed to be a factor of 4 (SSE Stride).
					Structure: Scalar.
					Type : 'Int'.
					Range : {1, 2, ...}.
- mConvKernel	-	Convolution Kernel.
					Structure: 2D Array.
					Type : 'Single'.
					Range : (-inf, inf).
- kernelNumRows	-	Kernel Number of Rows.
					Structure: Scalar.
					Type : 'Int'.
					Range : {1, 2, ...}.
- kernelNumCols -	Kernel Number of Columns.
					Assumbed to be a factor of 4 (SSE Stride).
					Structure: Scalar.
					Type : 'Int'.
					Range : {1, 2, ...}.
Reference:
1.	[].
Remarks:
1.	Actually applies Correlation and not Convolution.
2.	Using 'Replicate' / 'Nearest Neighbor' Boundary Conditions.
2.	It seems VS 2015 generates faster code than GCC/
TODO :
1.	Why the OpenMP DLL Generated in VS 2015 crashes MATLAB (GCC's DLL Works in MATLAB).
Release Notes:
-	1.0.000	04/08/2017	Royi Avital
*   First release version
*/
// -------------------------------------- ImageConvolution -------------------------------------- //
void ImageToColumns(DATA_TYPE* mO, DATA_TYPE* mI, int numRows, int numCols, int blockRadius)
{
	int blockSize, blockNumElements, colImageColIdx, colImageRowIdx;
	int ii, jj, kk, ll;

	blockSize			= (2 * blockRadius) + 1;
	blockNumElements	= blockSize * blockSize;
	colImageRowIdx		= -1;

// #pragma omp parallel for private(jj, colImageColIdx, kk, ll) firstprivate(blockRadius, numRows, numCols, colImageRowIdx)
	for (ii = blockRadius; ii < (numRows - blockRadius); ii++)
	{
		for (jj = blockRadius; jj < (numCols - blockRadius); jj++)
		{
			colImageRowIdx = colImageRowIdx + 1;
			colImageColIdx = -1;

			for (kk = -blockRadius; kk <= blockRadius; kk++)
			{
				for (ll = -blockRadius; ll <= blockRadius; ll++)
				{
					colImageColIdx = colImageColIdx + 1;
					// mO[(colImageRowIdx * blockNumElements) + colImageColIdx] = 1;
					mO[(colImageRowIdx * blockNumElements) + colImageColIdx] = mI[((ii + kk) * numCols) + jj + ll];
				}

			}
		}
	}


}
