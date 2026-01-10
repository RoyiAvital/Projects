<!-- https://visitor-badge.glitch.me/ -->
<!-- ![Visitors](https://visitor-badge.glitch.me/badge?page_id=RoyiAvital.StackExchangeCodes) -->
<!-- https://hits.seeyoufarm.com/ -->
[![Visitors](https://hits.seeyoufarm.com/api/count/incr/badge.svg?url=https%3A%2F%2Fgithub.com%2FRoyiAvital%2FStackExchangeCodes&count_bg=%2379C83D&title_bg=%23555555&icon=&icon_color=%23E7E7E7&title=Visitors+%28Daily+%2F+Total%29&edge_flat=false)](https://github.com/RoyiAvital/StackExchangeCodes)
[![Visitors](https://api.visitorbadge.io/api/combined?path=https%3A%2F%2Fgithub.com%2FRoyiAvital%2FStackExchangeCodes&labelColor=%23f47373&countColor=%23555555&style=plastic)](https://github.com/RoyiAvital/StackExchangeCodes) <!-- https://www.visitorbadge.io -->
[![DOI](https://zenodo.org/badge/44807437.svg)](https://zenodo.org/badge/latestdoi/44807437)
<a href="https://liberapay.com/Royi/donate"><img alt="Donate using Liberapay" src="https://liberapay.com/assets/widgets/donate.svg"></a>
[![ko-fi](https://ko-fi.com/img/githubbutton_sm.svg)](https://ko-fi.com/K3K8BLS2B)

# Curve Smoothing

Implementation of various 2D and 3D curve smoothing algorithms.

## Spline Like Optimization Model

Let $\left\{ {y}_{i} \in \mathbb{R} \right\}_{i = 1}^{n}$ a set of points on a uniformly sampled grid.  
The optimization problem generates a smoothed path  $\left\{ {x}_{i} \in \mathbb{R} \right\}_{i = 1}^{n}$ with the following properties:

 - Smooth curve by applying regularization on the derivative of the path.
 - Equality on a given set of reference indices $\mathcal{I}$: $\forall i \in \mathcal{I}: \; {x}_{i} = {y}_{i}$.
 - Monotonicity of the path within a segment between 2 reference points.  
   That is, $k, l \in \mathcal{I}, k < l, {y}_{k} \leq {y}_{l} \Rightarrow {y}_{l} \leq {y}_{l + 1} \leq \ldots \leq {y}_{k}$. 

The problem can be formulated into the following Quadratic Problem:

$$
\begin{alignat*}{3}
\arg \min_{ \boldsymbol{x} } & \quad & \frac{1}{2} \left\| \boldsymbol{x} - \boldsymbol{y} \right\|_{2}^{2} + \frac{\lambda}{2} {\left\| \boldsymbol{D}^{m} \boldsymbol{x} \right\|}_{2}^{2}  \\
\text{subject to} & \quad & {x}_{i} = {y}_{i} \; \forall i \in \mathcal{I} \\
& \quad & \boldsymbol{A} \boldsymbol{x} \leq \boldsymbol{0} \\
\end{alignat*}
$$

Where the matrix $\boldsymbol{A}$ forces the monotonic property of each segment.


This problem is a [Quadratic Programming](https://en.wikipedia.org/wiki/Quadratic_programming) (QP) as the Hessian of the problem is a single Symmetric Positive Definite (SPD) Matrix $\boldsymbol{I} + \lambda \boldsymbol{D}^{\top} \boldsymbol{D}$.  
Because the Hessian of the objective function is positive definite (Due to the $\boldsymbol{I}$ term), the problem is strictly convex and has a unique global minimum.

### ADMM Solve in ProxQP Style


The `ProxQP` ([PROXQP: an Efficient and Versatile Quadratic Programming Solver for Real Time Robotics Applications and Beyond](https://inria.hal.science/hal-04198663)) and `PIQP` ([PIQP: A Proximal Interior-Point Quadratic Programming Solver](https://arxiv.org/abs/2304.00290)) are based on the _Proximal Augmented Lagrangian Method_. This is slightly more robust than standard ADMM because it uses proximal regularization to ensure the sub problems are always well conditioned.  
This _extra_ feature is not needed in the above case as the problem is strictly convex to begin with.

#### Define the Quadratic Form

With $\boldsymbol{P} = \boldsymbol{I} + \lambda \boldsymbol{D}^{\top} \boldsymbol{D}$ the Quadratic Objective becomes:

$$ \frac{1}{2} \boldsymbol{x}^{\top} \boldsymbol{P} \boldsymbol{x} + \boldsymbol{q}^{\top} \boldsymbol{x} $$


#### The Augmented Lagrangian

The inequality $\boldsymbol{A} \boldsymbol{x} \leq \boldsymbol{0}$ is converted to equality constraint by introducing a slack variable $\boldsymbol{s} \geq \boldsymbol{0}$ such that $\boldsymbol{A} \boldsymbol{x} + \boldsymbol{s} = \boldsymbol{0}$.  
Using the matrix $\boldsymbol{E}$ one also represent the equality constraint as linear equation which yields the _Augmented Lagrangian_:

$$ \mathcal{L}_{\rho} \left( \boldsymbol{x}, \boldsymbol{s}, \boldsymbol{\mu}, \boldsymbol{\nu} \right) = \frac{1}{2} \boldsymbol{x}^{\top} \boldsymbol{P} \boldsymbol{x} + \boldsymbol{q}^{\top} \boldsymbol{x} + \frac{\rho}{2} {\left\| \boldsymbol{A} \boldsymbol{x} + \boldsymbol{s} + {\rho}^{-1} \boldsymbol{\mu} \right\|}_{2}^{2} + \frac{\rho}{2} {\left\| \boldsymbol{E} \boldsymbol{x} - \boldsymbol{d} + {\rho}^{-1} \boldsymbol{\nu} \right\|}_{2}^{2} $$

#### The Iterative Scheme

Like `ProxQP` the solution applies the primal update by solving a regularized KKT system.  

##### Primal Update

Minimizing the _Augmented Lagrangian_ with respect to $\boldsymbol{x}$ by solving a _Linear System_:

$$ \boldsymbol{x}^{\left( k + 1 \right)} = \boldsymbol{K}^{-1} \boldsymbol{r}^{\left( k \right)} $$

With $\boldsymbol{K} = \boldsymbol{P} + \rho \boldsymbol{A}^{\top} \boldsymbol{A} + \rho \boldsymbol{E}^{\top} \boldsymbol{E}$ and $\boldsymbol{r}^{\left( k \right)} = - \boldsymbol{q} - \boldsymbol{A}^{\top} \left( \rho \boldsymbol{s}^{\left( k \right)} + \boldsymbol{\mu}^{\left( k \right)} \right) + \boldsymbol{E}^{\top} \left( \rho \boldsymbol{d} - \boldsymbol{\nu} \right)$.

Since the matrix $\boldsymbol{K}$ is SPD and independent of the iteration its [Cholesky Decomposition](https://en.wikipedia.org/wiki/Cholesky_decomposition) can be pre calculated.


##### Slack Update

Minimizing $\mathcal{L}_{\rho}$ with respect to $\boldsymbol{s}$ subject to $\boldsymbol{s} \geq \boldsymbol{0}$:

$$ \boldsymbol{s}^{\left( k + 1 \right)} = \max \left(0, - \boldsymbol{A} \boldsymbol{x}^{\left( k + 1 \right)} - {\rho}^{-1} \boldsymbol{\mu}^{\left( k \right)} \right) $$

This is the _Prox_ part, basically a projection.


##### Dual Updates

$$
\begin{align*}
\boldsymbol{\mu}^{\left( k + 1 \right)} & = \boldsymbol{\mu}^{\left( k \right)} + \rho \left( \boldsymbol{A} \boldsymbol{x}^{\left( k + 1 \right)} + \boldsymbol{s}^{\left( k + 1 \right)} \right) \\
\boldsymbol{\nu}^{\left( k + 1 \right)} & = \boldsymbol{\nu}^{\left( k \right)} + \rho \left( \boldsymbol{E} \boldsymbol{x}^{\left( k + 1 \right)} - \boldsymbol{d} \right)
\end{align*}
$$

#### Difference from the ProxQP Algorithm

The `ProxQP` adds another term to the _Augmented Lagrangian_: $\frac{\eta}{2} {\left\| \boldsymbol{x} - \boldsymbol{x}^{\left( k \right)} \right\|}^{2}$.   
It yields a different matrix $\boldsymbol{K} = \boldsymbol{P} + \rho \boldsymbol{A}^{\top} \boldsymbol{A} + \rho \boldsymbol{E}^{\top} \boldsymbol{E} + \eta \boldsymbol{I}$ which is an SPD even if $\boldsymbol{Q}$ is Symmetric Semi Positive Definite (SPSD).  
Yet, in the case above it is unnecessary. 

#### Convergence Tests

Define 2 types of residuals:

 * Primal Residual: ${r}^{\left( k \right)} = {\left\| \boldsymbol{A} \boldsymbol{z}^{\left( k \right)} + \boldsymbol{s}^{\left( k \right)} \right\|}_{\infty}$.
 * Dual Residual: ${s}^{\left( k \right)} = {\left\| \rho \boldsymbol{A}^{\top} \left( \boldsymbol{s}^{\left( k \right)} - \boldsymbol{s}^{\left( k - 1 \right)} \right) \right\|}_{\infty}$.

The use of the ${L}^{\infty}$ Norm decouples the residual value from the dimensions of the vectors.

The convergence test is given by: ${r}^{\left( k \right)} \leq {\varepsilon}_{\text{Abs}} \land {s}^{\left( k \right)} \leq {\varepsilon}_{\text{Abs}}$ for a threshold ${\varepsilon}_{\text{Abs}} > 0$.

> [!TIP]
> The tests are computational intensive. Hence one might apply them once in a few iterations.

> [!NOTE]
> Some use additional relative test in addition to the absolute test.

> [!CAUTION]
> This test the _Primal_ / _Dual_ _Residual_ yet not the _Primal Dual Gap_.  
> In the `ProxQP` and `PIQP` papers they also show the gap calculation.

#### Adaptive ρ

Faster convergence happens when the ratio between the _Primal Residual_ and _Dual Residual_ is close to $1$.  
In order to achieve this the penalty parameter $\rho$ should be adjusted.  

The logic of the update:

 - If the _Primal Residual_ is much larger, increase the value of $\rho$. Make solution obey the constraints.
 - If the _Dual Residual_ is much larger, decrease the value of $\rho$. Make exploration of better solutions easier.

Copying the logic in `OSQP` ([OSQP: An Operator Splitting Solver for Quadratic Programs](https://arxiv.org/abs/1711.08013)):

$${\rho}^{\left( k + 1 \right)} \gets {\rho}^{\left( k \right)} \sqrt{\frac{{r}^{\left( k \right)}}{{s}^{\left( k \right)}}}$$

Once the stiffness parameter is updated, the matrix $\boldsymbol{K}$ and its decomposition should also be updated.  
One heuristic to test if one should update the parameter if $\frac{{r}^{\left( k \right)}}{{s}^{\left( k \right)}} \geq \tau \lor \frac{{s}^{\left( k \right)}}{{r}^{\left( k \right)}} \geq \tau$.

> [!TIP]
> The calculation are computational intensive. Hence one might apply them once in a few iterations.
