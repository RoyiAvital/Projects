# Projects
Various Small Projects on Various Subjects

## Image Convolution
Implementation of Image Convolution using Multi Threading and SIMD (Vectrorization) for speed.
The project generates DLL / DYLIB / SO to be used in MATLAB.

## Image to Columns
Optimized C based implementation of MATLAB's [`im2col()`](https://www.mathworks.com/help/images/ref/im2col.html) function.

## Levinson Recursion
Implementation of the [Levinson Recursion](https://en.wikipedia.org/wiki/Levinson_recursion) on MATLAB.  
This allows solving Linear System `A x = b` where `A` is a [Toeplitz Matrix](https://en.wikipedia.org/wiki/Toeplitz_matrix) in O(N^2) complexity instead of O(N^3) using classic solution.

## Optimization
Various projects which are Optimization (Numerical Methods) oriented such as:

 * Projection into Balls (Like Norm Balls, the Simplex Ball, etc...).
 * Interior Point () - Implementation in MATLAB.
 * Steepest Descent - Calculation of the Steepest Descent Direction for various norms ($ {L}_{1} $, $ {L}_{2} $, $ {L}_{\infty} $).

## Singular Value Decomposition
A Presentation about [Singular Value Decomposition](https://en.wikipedia.org/wiki/Singular_value_decomposition) (SVD).  
It covers the following:

 * Singular Value Theorem.
 * Applications
    * Order Reduction.
    * Solving Linear Equation System (Least Squares).
    * Total Least Squares.
    * [Principal Component Analysis](https://en.wikipedia.org/wiki/Principal_component_analysis) (PCA).

The PDF can be viewed directly - [SVD Presentation PDF](https://docs.google.com/viewer?url=https://github.com/RoyiAvital/Projects/raw/master/SingularValueDecomposition/SVD.pdf).
