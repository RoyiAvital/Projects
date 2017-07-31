# Exam A (27/07/2017) Question 003

## The Problem

\begin{align*}
    \text{minimize}     & \quad & \left\| x \right\|^{2} \\
    \text{subject to}   & \quad & x - a \leq 0 \\
    \text{}             & \quad & \boldsymbol{1}^{T} x = b
\end{align*}

## Solution

For the solution it will be assumed, for simplicity, that $ {a}_{1} \leq {a}_{2} \leq \ldots \leq {a}_{n} $.

The Lagrangian of the problem is given by:

$$ L \left( x, \lambda, \nu \right) = {x}^{T} x + {\lambda}^{T} \left( x - a \right) + \nu \left( \boldsymbol{1}^{T} x - b \right) $$

The KKT System is given by:

\begin{align*}
    \nabla L \left( x, \lambda, \nu \right) = 2x + \lambda + \nu \boldsymbol{1} & = 0 & \text{(1)} \\
    {\lambda}_{i} \left( {x}_{i} - {a}_{i} \right) & = 0 & \text{(2)} \\
    \boldsymbol{1}^{T} x & = b & \text{(3)} \\
    \lambda & \succeq 0 & \text{(4)}
\end{align*}

Define the sets $ A = \left\{ i \mid {x}_{i} \neq {a}_{i} \right\} $ and $ B = \left\{ i \mid {x}_{i} = {a}_{i} \right\} $.  
Looking at $ A $ yields $ \forall j \in A, \: {\lambda}_{j} = 0 $ which means:

\begin{align*}
    2 {x}_{A} + \nu \boldsymbol{1} & = 0 & \text{} \\
    & \Rightarrow 2 \boldsymbol{1}^{T} {x}_{A} + \nu \boldsymbol{1}^{T} \boldsymbol{1} = 0 & \text{Summing both ends} \\
    & = 2 \left( b - \boldsymbol{1}^{T} {x}_{B} \right) + \nu \left| A \right| & \text{By (3)} \\
    & \Rightarrow \nu = -2 \frac{b - \boldsymbol{1}^{T} {x}_{B}}{\left| A \right|} & \text{} \\
    & \Rightarrow {{x}_{A}}_{i} = \frac{b - \boldsymbol{1}^{T} {x}_{B}}{\left| A \right|}
\end{align*}

Where $ {x}_{A} $ is a vector composed of elements which are in set $ A $ by their indices.  
The intuition is that for any $ {x}_{i} $ such that $ i \in A $ the value is the mean of the residual $ b - \boldsymbol{1}^{T} {x}_{B} $.

Looking at $ B $ suggests $ \forall i \in B, \: {x}_{i} = {a}_{i} $ which yields:

\begin{align*}
    2 {a}_{i} + {\lambda}_{i} - 2 \frac{b - \boldsymbol{1}^{T} {x}_{B}}{\left| A \right|} & = 0 & \text{By (1)} \\
    & \Rightarrow {\lambda}_{i} = 2 \left( \frac{b - \boldsymbol{1}^{T} {x}_{B}}{\left| A \right|} - {a}_{i} \right)
\end{align*}

From constraint (4) the above is valid only in the case:
$$ i : {\lambda}_{i} = 2 \left( \frac{b - \boldsymbol{1}^{T} {x}_{B}}{\left| A \right|} - {a}_{i} \right) \geq 0 $$

Since, for the case with no constraints of $ a $ the solution would be given by the mean, namely $ x = \frac{b}{n} \boldsymbol{1} $ the logic is to increase $ \left| A \right| $.

Since $ {a}_{i} $ are sorted in ascending order one could search for the following:

 * Start with $ A = \left\{ 1, 2, \ldots, n \right\} $ and $ B = \emptyset $.
 * Do Loop $ i = 1, 2, \ldots n $:
     * If $ \left( \frac{b - \boldsymbol{1}^{T} {x}_{B}}{\left| A \right|} \leq {a}_{ {A}_{i} } \right) $ break the loop.
     * Move $ {A}_{i} $ to $ B $.
 * Calculate $ x $ as above.

Remarks:

 * If at the end $ A = \emptyset $ the problem isn't feasible as feasibility requires $ \boldsymbol{1}^{T} a \geq b $.
 * The solution is non negative, namely $ x \geq 0 $ assuming $ a \geq 0 $.


 

