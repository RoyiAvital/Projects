# Total Variation (TV) Regularized Least Squares - Solvers Analysis
This project analyze the performance and implementation of solvers which solves the ${L}_{1}$ Regularized Least Squares problem:

$$ \arg \min_{ x \in \mathbb{R}^{n} } \frac{1}{2} {\left\| A x - b \right\|}_{2}^{2} + \lambda {\left\| x \right\|}_{TV} $$

Where the [Total Variation](https://en.wikipedia.org/wiki/Total_variation) (TV) norm is given by:

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

Usually the solvers will solve for any Matrix $ D $ unless stated otherwise.

For clarity we'll define:

$$ F \left( x \right) = \frac{1}{2} {\left\| A x - b \right\|}_{2}^{2} + \lambda {\left\| D x \right\|}_{1}, \; f \left( x \right) = \frac{1}{2} {\left\| A x - b \right\|}_{2}^{2}, \; g \left( x \right) = {\left\| D x \right\|}_{1} $$

## Solvers

### Sub Gradient Method (SGM)

The sub graduent is given by:

$$ \frac{\partial F \left( x \right ) }{\partial x} = {A}^{T} \left( A x - b \right ) + \lambda {D}^{T} \operatorname{sign} \left( D x \right ) $$

### The Chambolle Method (Chambolle Pock)

$$\begin{aligned}
\arg \min_{ x \in \mathbb{R}^{n} } \frac{1}{2} {\left\| A x - b \right\|}_{2}^{2} + \lambda {\left\| D x \right\|}_{1} & \overset{1}{=} \arg \min_{ x \in \mathbb{R}^{n} } \max_{ {\left\| p \right\|}_{\infty} \leq 1 } \frac{1}{2} {\left\| A x - b \right\|}_{2}^{2} + \lambda {p}^{T} D x \\
& \overset{2}{\Rightarrow} \arg \max_{ {\left\| p \right\|}_{\infty} \leq 1 } \min_{ x \in \mathbb{R}^{n} } \frac{1}{2} {\left\| A x - b \right\|}_{2}^{2} + \lambda {p}^{T} D x \\
& \overset{3}{=} \arg \max_{ {\left\| p \right\|}_{\infty} \leq 1 } \frac{1}{2} {\left\| A {x}^{\ast} - b \right\|}_{2}^{2} + \lambda {p}^{T} D {x}^{\ast} \\
& \overset{4}{=} \arg \min_{ {\left\| p \right\|}_{\infty} \leq 1 } -\frac{1}{2} {\left\| A {x}^{\ast} - b \right\|}_{2}^{2} - \lambda {p}^{T} D {x}^{\ast}
\end{aligned}$$

Where:

 1. As the Dual Norm (Support Function) of $ {\left\| x \right\|}_{1} = \sup \left\{ {p}^{T} x  \mid {\left\| p \right\|}_{\infty} \leq 1 \right\} $.
 2. By the [Min Max Theorem](https://en.wikipedia.org/wiki/Minimax_theorem) (The objective is Convex in 𝑥 and Concave in 𝑝) one could switch the order of the Maximum and Minimum.
 3. Setting the optimal $ {x}^{\ast} = {\left( {A}^{T} A \right)}^{-1} \left( {A}^{T} b - \lambda {D}^{T} p \right) $.
 4. Switching Maximization to Minimization by negation of the function.

Now, the function $ \frac{1}{2} {\left\| A {x}^{\ast} - b \right\|}_{2}^{2} + \lambda {p}^{T} D {x}^{\ast} $ is concave and smooth in $ p $ and constraint $ {\left\| p \right\|}_{\infty} \leq 1 $ is convex set. Hence the problem can be solved using Projected Gradient Descent.  
The gradient is given by:

$$ {\nabla}_{p} \left( \frac{1}{2} {\left\| A {x}^{\ast} - b \right\|}_{2}^{2} + \lambda {p}^{T} D {x}^{\ast} \right) = \lambda D {\left( {A}^{T} A \right)}^{-1} {A}^{T} b - {\lambda}^{2} D {\left( {A}^{T} A \right)}^{-1} {D}^{T} p $$

For clarity we'll negate the function (And the derivative) to work with Convex Function over Convex Set in a minimization problem.

So the step is given by:

$$ {p}^{k + 1} = \mathbf{P}_{{\left\| p \right\|}_{\infty} \leq 1} \left( {p}^{k} - {t}_{k} \left( {\lambda}^{2} D {\left( {A}^{T} A \right)}^{-1} {D}^{T} p - \lambda D {\left( {A}^{T} A \right)}^{-1} {A}^{T} b \right) \right) $$

Where $ {t}_{k} $ is the step size of the $ k $ -th step and $ \mathbf{P}_{{\left\| p \right\|}_{\infty} \leq 1} \left( \cdot \right) $ is the projection onto the set which is given by:

$$ \mathbf{P}_{{\left\| p \right\|}_{\infty} \leq 1} \left( x \right)_{i} = \frac{ {x}_{i} }{ \max \left\{ 1, \left| {x}_{i} \right| \right\} } $$

### Alternating Direction Method of Multipliers (ADMM)

The problem is rewritten as:

$$\begin{aligned}
\arg \min_{x, z} & \quad \frac{1}{2} {\left\| A x - b \right\|}_{2}^{2} + \lambda {\left\| z \right\|}_{1} \\
\text{subject to} & \quad D x = z
\end{aligned}$$

The Augmented Lagrangian is given by:

$$ {L}_{\rho} \left( x, z, u \right) = \frac{1}{2} {\left\| A x - b \right\|}_{2}^{2} + \lambda {\left\| z \right\|}_{1} + \rho {u}^{T} \left( D x - z \right) + \frac{\rho}{2} {\left\| D x - z \right\|}_{2}^{2} $$

The steps are:

 1. $ {x}^{k + 1} = \arg \min_{x} \frac{1}{2} {\left\| A x - b \right\|}_{2}^{2} + \rho {u}^{T} \left( D x - z \right) + \frac{\rho}{2} {\left\| D x - z \right\|}_{2}^{2} $ which is given by $ {x}^{k + 1} = {\left( {A}^{T} A + \rho {D}^{T} D \right)}^{-1} \left( {A}^{T} b + \rho {D}^{T} \left( {z}^{k} - {u}^{k} \right) \right) $.
 2. $ {z}^{k + 1} = \arg \min_{z} \lambda {\left\| z \right\|}_{1} + \rho {u}^{T} \left( D x - z \right) + \frac{\rho}{2} {\left\| D x - z \right\|}_{2}^{2} $. Since minimizing $ \rho {u}^{T} \left( D x - z \right) + \frac{\rho}{2} {\left\| D x - z \right\|}_{2}^{2} $ with respect to $ z $ is equivalent to minimizing $ \frac{\rho}{2} {\left\| D x - z + u \right\|}_{2}^{2} $ (Open the expression and the only difference is $ {\left\| u \right\|}_{2}^{2} $ which has no effect) then $ {z}^{k + 1} = \arg \min_{z} \lambda {\left\| z \right\|}_{1} + \frac{\rho}{2} {\left\| D x + u - z \right\|}_{2}^{2} $ which in the ADMM iteration means $ {z}^{k + 1} = \mathbf{S}_{ \frac{\lambda}{\rho} } \left( D {x}^{k + 1} + {u}^{k} \right) $ where $ \mathbf{S}_{\lambda} \left( \cdot \right) $ is the Soft Threshold operator with parameter $ \lambda $.
 3. $ {u}^{k + 1} = {u}^{k} + D {x}^{ k + 1} - {z}^{k + 1} $.

## Results

![](.\Figure0001.png)

## TO DO List
 *  Look for others methods.
 *  [Total Variation Denoising (An MM Algorithm)](https://eeweb.engineering.nyu.edu/iselesni/lecture_notes/TVDmm/).  
    Work of [Ivan W. Selesnick](https://eeweb.engineering.nyu.edu/iselesni/).
 *  
