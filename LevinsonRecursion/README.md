# Levinson Recursion


## Introduction
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
 * MATLAB R2017a.
    * BLAS Version (`version -blas`) - `Intel(R) Math Kernel Library Version 11.3.1 Product Build 20151021 for Intel(R) 64 architecture applications, CNR branch AVX2`.
    * LAPACK Version (`version -lapack`) - `Intel(R) Math Kernel Library Version 11.3.1 Product Build 20151021 for Intel(R) 64 architecture applications, CNR branch AVX2; Linear Algebra PACKage Version 3.5.0`.
 * Windows 10 Professional 64 Bit.

## How to Run
Compile the `LevinsonRecursionToeplitzMatrix.cpp` file using MSVC / GCC (See `CompileDllGcc.m` for reference).


## ToDo
 * Update the README file about the method.
 * Improve the C code (Better `h` file and more established code).