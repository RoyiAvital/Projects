/*
 * Array Sorting
 * Sorts a contiguous array using Sorting Network and High Level sorting algorithm.
 *
 * References:
 *	1. 	Based on https://github.com/eggpi/sorting-networks-test by Guilherme P. Goncalves.
 * Remarks:
 *	1.	Supports `unisgned int`, `int`, `single` and `double` by changing `ELEMENT_TYPE_IDX`.
 * TODO:
 *	1.	Add option for Merge Sort in addition to Quick Sort.
 * Release Notes:
 *	-	1.0.000	13/07/2020	Royi Avital
 *		*	First release version.
 */

#ifndef ELEMENT_TYPE_IDX
#define ELEMENT_TYPE_IDX 1
#endif // !ELEMENT_TYPE_IDX


#if ELEMENT_TYPE_IDX == 1
#define ELEMENT_TYPE unsigned int
#elif ELEMENT_TYPE_IDX == 2
#define ELEMENT_TYPE int
#elif ELEMENT_TYPE_IDX == 3
#define ELEMENT_TYPE float
#elif ELEMENT_TYPE_IDX == 4
#define ELEMENT_TYPE double
#endif // ELEMENT_TYPE_IDX == 1


#define SWAP(vA, ii, jj) \
if ((vA)[ii] > (vA)[jj]) { \
	ELEMENT_TYPE tmpVal; \
	tmpVal = (vA)[ii]; \
	(vA)[ii] = (vA)[jj]; \
	(vA)[jj] = tmpVal; \
}

#define NSORT2(vA) \
SWAP((vA), 0, 1)

#define NSORT3(vA) \
SWAP((vA), 1, 2); SWAP((vA), 0, 2); SWAP((vA), 0, 1);

#define NSORT4(vA) \
SWAP((vA), 0, 1); SWAP((vA), 2, 3); SWAP((vA), 0, 2); \
SWAP((vA), 1, 3); SWAP((vA), 1, 2);

#define NSORT5(vA) \
SWAP((vA), 0, 1); SWAP((vA), 3, 4); SWAP((vA), 2, 4); \
SWAP((vA), 2, 3); SWAP((vA), 0, 3); SWAP((vA), 0, 2); \
SWAP((vA), 1, 4); SWAP((vA), 1, 3); SWAP((vA), 1, 2);

#define NSORT6(vA) \
SWAP((vA), 1, 2); SWAP((vA), 0, 2); SWAP((vA), 0, 1); \
SWAP((vA), 4, 5); SWAP((vA), 3, 5); SWAP((vA), 3, 4); \
SWAP((vA), 0, 3); SWAP((vA), 1, 4); SWAP((vA), 2, 5); \
SWAP((vA), 2, 4); SWAP((vA), 1, 3); SWAP((vA), 2, 3);

#define NSORT7(vA) \
SWAP((vA), 1, 2); SWAP((vA), 0, 2); SWAP((vA), 0, 1); \
SWAP((vA), 3, 4); SWAP((vA), 5, 6); SWAP((vA), 3, 5); \
SWAP((vA), 4, 6); SWAP((vA), 4, 5); SWAP((vA), 0, 4); \
SWAP((vA), 0, 3); SWAP((vA), 1, 5); SWAP((vA), 2, 6); \
SWAP((vA), 2, 5); SWAP((vA), 1, 3); SWAP((vA), 2, 4); \
SWAP((vA), 2, 3);

#define NSORT8(vA) \
SWAP((vA), 0, 1); SWAP((vA), 2, 3); SWAP((vA), 0, 2); \
SWAP((vA), 1, 3); SWAP((vA), 1, 2); SWAP((vA), 4, 5); \
SWAP((vA), 6, 7); SWAP((vA), 4, 6); SWAP((vA), 5, 7); \
SWAP((vA), 5, 6); SWAP((vA), 0, 4); SWAP((vA), 1, 5); \
SWAP((vA), 1, 4); SWAP((vA), 2, 6); SWAP((vA), 3, 7); \
SWAP((vA), 3, 6); SWAP((vA), 2, 4); SWAP((vA), 3, 5); \
SWAP((vA), 3, 4);


int ArrayPartition(ELEMENT_TYPE *vA, unsigned int s, unsigned int e) {
	unsigned int ii, jj, p;
	ELEMENT_TYPE tmpVal;

	/* Get pivot from median of three */
	p = (s + e)/2;

	if (vA[p] > vA[s]) {
		if (vA[e] > vA[s]) {
			/* p > e > s or e > p > s */
			p = (vA[p] > vA[e]) ? e : p;
		} else p = s;
	} else {
		if (vA[e] > vA[p]) {
			/* s > e > p or e > s > p */
			p = (vA[e] > vA[s]) ? s : e;
		}
	}

	/* Place pivot in the beginning */
	tmpVal = vA[p];
	vA[p] = vA[s];
	vA[s] = tmpVal;

	/* Partition <= | pivot | >= */
	ii = s;
	jj = e;
	while (1) {
		while (ii < e && vA[ii] <= vA[s]) ii++;
		while (jj > s && vA[jj] >= vA[s]) jj--;

		if (ii > jj) break;

		tmpVal = vA[ii];
		vA[ii] = vA[jj];
		vA[jj] = tmpVal;

		ii++;
		jj--;
	}

	/* Swap pivot back to its place */
	tmpVal = vA[--ii];
	vA[ii] = vA[s];
	vA[s] = tmpVal;

	return ii;
}

static void ArrayQuickSort(ELEMENT_TYPE *vA, unsigned int s, unsigned int e) {
	unsigned int p, size;

	size = e - s + 1;
	if (size <= 1) {
		return;
	}

	switch (size) {
		case 2:
			NSORT2(vA + s);
			break;
		case 3:
			NSORT3(vA + s);
			break;
		case 4:
			NSORT4(vA + s);
			break;
		case 5:
			NSORT5(vA + s);
			break;
		case 6:
			NSORT6(vA + s);
			break;
		case 7:
			NSORT7(vA + s);
			break;
		case 8:
			NSORT8(vA + s);
			break;
		default:
			p = ArrayPartition(vA, s, e);
			ArrayQuickSort(vA, s, p - 1);
			ArrayQuickSort(vA, p + 1, e);
	}

	return;
}

void ArrayInsertionSort(ELEMENT_TYPE A[], unsigned int numElements) {
	// Insertion Sort
	// Optimized for partially sorted array
	// https://opendsa-server.cs.vt.edu/ODSA/Books/Everything/html/SortOpt.html
	
	unsigned int ii, jj;
	ELEMENT_TYPE tmpVal;

	for (ii = 1; ii < numElements; ii++) { // Insert i'th record
		tmpVal = A[i];
		for (jj = ii; (jj > 0) && (tmpVal < A[jj - 1]); j--)
			A[jj] = A[jj - 1];
		A[jj] = tmpVal;
	}
}

