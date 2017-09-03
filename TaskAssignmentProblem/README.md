# Task Assignment Problem
This problem is a variation of the [Knapsack Problem](https://en.wikipedia.org/wiki/Knapsack_problem).

## The Problem Statement
Given the following model.  
There is a set of $ n $ workers, each with a given capacity $ \left\{ {c}_{j} \right\}_{j = 1}^{n} $.  
There are $ m $ tasks where the same task can be done by many workers or none at all. Each task has different gain and cost for each worker.

 1. The gain of completion of the $ i $ -th task by the $ j $ -th worker is given by $ {g}_{ij} $.  
 2. The work of doing the $ i $ -th task by the $ j $ -th worker is given by $ {w}_{ij} $.  
 3. The flag whether the $ i $ -th task is assigned to the $ j $ -th worker is given by $ {x}_{ij} $. Where $ {x}_{ij} = 1 $ means the work is assigned and $ {x}_{ij} = 0 $ means it is not assigned.

The objective is to maximize the gain of the assigned tasks according to the following constraints:

 * Each task can be either not be assigned at all or be allocated to at least $ k $ workers and not more than $ l $ workers.
 * No worker will exceed its work capacity.

This can be formulated as:

$$
\begin{alignat}{3}
\text{maximize} & \quad & \sum_{i = 1}^{m} \sum_{j = 1}^{n} {x}_{i j} {g}_{i j} \\
\text{subject to} & \quad & \sum_{i = 1}^{m} {x}_{i j} {w}_{i j} \leq {c}_{j} & \quad \forall j \\
& \quad & \sum_{j} {x}_{i, j} \leq l & \quad & \forall i \in T \\
& \quad & \sum_{j} {x}_{i, j} \geq k & \quad & \forall i \in T \\
& \quad & {x}_{i, j} \in \left\{ 0, 1 \right\}
\end{alignat}
$$

Where $ T = \left\{ i \mid \exists j : {x}_{i, j} = 1 \right\} $, namely the set of tasks which are assigned.

In the above, $ i $ (Rows of the matrix) is the task index and $ j $ is the worker index.


## Literature Review
Optimized C based implementation of MATLAB's [`im2col()`](https://www.mathworks.com/help/images/ref/im2col.html) function.

## Solutions
Implementation of the [Levinson Recursion](https://en.wikipedia.org/wiki/Levinson_recursion) on MATLAB.  
This allows solving Linear System `A x = b` where `A` is a [Toeplitz Matrix](https://en.wikipedia.org/wiki/Toeplitz_matrix) in O(N^2) complexity instead of O(N^3) using classic solution.

## References
 * [Knapsack Problem - Wikipedia](https://en.wikipedia.org/wiki/Knapsack_problem). 
 * [List of Knapsack Problem](https://en.wikipedia.org/wiki/List_of_knapsack_problems).
 * [David Pisinger's Optimization Codes](http://www.diku.dk/~pisinger/codes.html)  
 Code for many variations of the Knapsack Problem.
 * fd