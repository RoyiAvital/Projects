# ${L}_{1}$ Regularized Least Squares - Solvers Analysis
This project analyze the performance and implementation of solvers which solves the ${L}_{1}$ Regularized Least Squares problem:

$$ \arg \min_{x} \frac{1}{2} {\left\| A x - b \right\|}_{2}^{2} + \lambda {\left\| x \right\|}_{1} $$

## Solvers

### Sub Gradient Method (SGM)
TBC

### Proximal Gradient Method (PGM)
TBC

### Smoothing (Huber)
TBC

### Alternating Direction Method of Multipliers (ADMM)
TBC

### Iterative Reweighted Least Squares (IRLS)
TBC

### Coordinate Descent (CD)
TBC

## Results

![](Figure0001.png)

## TO DO List
 *  Check the performance of Proximal Gradient Method with Line Search.
 *  Add the [PANOC Algorithm (Proximal Averaged Newton Type Method for Optimal Control)](https://arxiv.org/abs/1709.06487).
 *	Add Newton Based solver (See https://github.com/cvxgrp/l1_ls).
