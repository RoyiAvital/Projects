# Prox Functions

This repository is dedicated to have functions which apply the $ \operatorname{Prox} \left( \cdot \right) $ operator.

The $ \operatorname{Prox} \left( \cdot \right) $ is given by:

$$ \operatorname{Prox}_{\lambda g \left( \cdot \right)} \left( y \right) = \arg \min_{x} \frac{1}{2} {\left\| x - y \right\|}_{2}^{2} + \lambda g \left( x \right) $$

The funcitons in this peoject solve the problem for various functions $ g \left( \cdot \right) $.

In thi project we'll assume $ x \in \mathbb{R}^{n} $,

## Total Variation

[Total Variation](https://en.wikipedia.org/wiki/Total_variation) (TV) norm is given by:

$$ {\left\| x \right\|}_{TV} = \sum_{i = 2}^{n} \left| {x}_{i} - {x}_{i - 1} \right| $$

This yields a formulation as given by:

$$ \arg \min_{ x \in \mathbb{R}^{n} } \frac{1}{2} {\left\| A x - b \right\|}_{2}^{2} + \lambda {\left\| D x \right\|}_{1} $$

Where $ D \in \mathbb{R}^{ \left( n - 1 \right) \times n } $ is the Finite Differences operator given by:

$$ \begin{bmatrix}
-1 & 1 & 0 & \ldots & 0 & 0 \\ 
0 & -1 & 1 & 0 & \ldots & 0 \\ 
\vdots & \vdots & \ddots &  & \vdots & \vdots \\ 
0 & \ldots &  & 0 & -1 & 1
\end{bmatrix} $$

## TO DO List
 *  Look for others methods.
