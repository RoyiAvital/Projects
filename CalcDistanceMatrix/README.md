# Calculate Distance Matrix

## Introduction
This project tries to compare different methods to calculate the Distance Matrix between 2 set of vectors.  
The projects compares 4 different approaches:
 *  Using MATLAB.
 *  Using Vanilla C Implementation.  
    The C Implementation will compiled using 3 different compilers (MSVC 2017 15.6.6, Intel Compiler 18.1, GCC 7,3) with highest optimization levels.
 *  Using SSE / AVX Intrinsics Code.  
    Using SSE / AVX intrinsics code to be compleletly insensitive to compiler.
 *  Using Eigen.


## Implementation Notes

There are 2 matrices $ A \in \mathbb{R}^{k \times m} $, $ B \in \mathbb{R}^{k \times n} $ which are a set of vectors in $ \mathbb{R}^{k} $.  
The distance matrix $ D \in \mathbb{R}^{m \times n} $ is given by:

$$ {D}_{ij} = {\left\| {A}_{i} - {B}_{j} \right\|}_{2}^{2} $$

Where $ {A}_{i} $ is the $ i $ -th vector in the set $ A $ (Namely the $ i $ -th column os the matrix).  

In matrix form the first index is the one which is continious in memory.

The [Levinson Recursion](https://en.wikipedia.org/wiki/Levinson_recursion) can be used for solving Linear System of Equations `A x = y` in O(N^2) assuming `A` is [Toeplitz Matrix](https://en.wikipedia.org/wiki/Toeplitz_matrix).  
Classic solvers, which doesnt take advatnage of the structure of `A` have O(N^3) complexity.

The method is basically the [Levinson Durbin Recursion](https://www.mathworks.com/help/signal/ref/levinson.html) for estimating [Auto Regressive Model](https://en.wikipedia.org/wiki/Autoregressive_model) Parameters.  
Yet in MATLAB it is limited to the case matrix `A` is composed from `y`.

## Results

![Matrix Generation](https://raw.githubusercontent.com/RoyiAvital/Projects/master/LevinsonRecursion/Figure0001.png)


### System Configuration
 * CPU - Intel Core I7 6800K @ 3.4 [GHz].
 * Memory - 4 * 8 [GB] @ 2166 [MHz] - G.Skill F4 2800C-16-8GRK.
 * Mother Board - ASRock X99 Killer (BIOS Version P3.20).
 * MATLAB R2018a.
    * BLAS Version (`version -blas`) - `Intel(R) Math Kernel Library Version 2017.0.31 Product Build 20170606 for Intel(R) 64 architecture applications, CNR branch AVX2`.
    * LAPACK Version (`version -lapack`) - `Intel(R) Math Kernel Library Version 2017.0.31 Product Build 20170606 for Intel(R) 64 architecture applications, CNR branch AVX2 Linear Algebra PACKage Version 3.7.0`.
 * Windows 10 Professional 64 Bit (`Version 1709 (OS Build 16299.371)`).

## How to Run
Compile the `LevinsonRecursionToeplitzMatrix.cpp` file into a DLL using MSVC / GCC (See `CompileDllGcc.m` for reference).  
Run `LevinsonRecursionRunTime.m`.


## ToDo
 * Update the README file about the method.
 * Improve the C code (Better `h` file and more established code).