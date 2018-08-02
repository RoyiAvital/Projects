# Derivation of TDOA CRLB
The derivation of the [Cramer Rao Lower Bound][1] for TDOA Estimation.

## Definition

 *	Set of  $ M $ sensors located at $ \left[ \boldsymbol{q}_{1}, \boldsymbol{q}_{2}, \ldots, \boldsymbol{q}_{M} \right] $ where $ \boldsymbol{q}_{i} \in \mathbb{R}^{d} $. 
 *	An Emitter at location $ \boldsymbol{p} \in \mathbb{R}^{d} $.
 *	Wave propagation time given by $ v $.
 *  The distance function $ {d}_{i} \left( \boldsymbol{p} \right) = {\left\| \boldsymbol{p} - \boldsymbol{q}_{i} \right\|}_{2} $.
 *  The set of sensors pairs $ I = \left\{ \left( i, j \right) \mid 1 \leq i < j \leq M \right\} $ where $ \left| I \right| = N = \frac{M \left( M - 1 \right)}{2} $.
 *  The cross distance function $ {d}_{ij} \left( \boldsymbol{p} \right) = {d}_{i} \left( \boldsymbol{p} \right) - {d}_{j} \left( \boldsymbol{p} \right) $.
 *	The cross distance vector from the emitter to each sensor $ \boldsymbol{d} \in \mathbb{R}^{N} $ given by $ \boldsymbol{d} = {\left[ {d}_{ij} \left( p \right), \ldots \right]}_{ \left(i, j \right) \in I }^{T} $.
 * 	The measured time difference vector $ \boldsymbol{ \tau } \in \mathbb{R}^{N} $ given by $ \boldsymbol{ \tau } = \frac{1}{v} \boldsymbol{d} + \boldsymbol{n} $ where the noise vector $ \boldsymbol{n} \in \mathbb{R}^{N} $ is given by $ \boldsymbol{n} = {\left[ {n}_{ij}, \ldots \right]}_{ \left(i, j \right) \in I }^{T} $.

## Derivation

Under the assumption the noise vector $ \boldsymbol{n} $  is Gaussian with zero mean and the covariance matrix $ {\sigma}^{2} I $ the measurement vector $ \boldsymbol{ \tau } \sim \mathcal{N} \left( \frac{1}{v} \boldsymbol{d} \left( \boldsymbol{p} \right), {\sigma}^{2} I \right) $.

The [Fischer Information Matrix][2] for Multivariate Normal Distribution $ \boldsymbol{x} \sim \mathcal{N} \left( \boldsymbol{\mu} \left( \boldsymbol{\theta} \right), C \right) $ is given by:

$$ {I}_{m, k} \left( \boldsymbol{\theta} \right) = \frac{ \partial \boldsymbol{\mu}^{T} \left( \boldsymbol{\theta} \right) }{ \partial \boldsymbol{\theta}_{m} } {C}^{-1} \frac{ \partial \boldsymbol{\mu}^ \left( \boldsymbol{\theta} \right) }{ \partial \boldsymbol{\theta}_{k} } $$

In the case above it becomes:

$$\begin{align*}
I \left( \boldsymbol{\theta} \right) & = \frac{1}{v} \frac{\partial \boldsymbol{d}^{T} \left( \boldsymbol{p} \right) }{d \boldsymbol{p}} \frac{1}{ {\sigma}^{2} } \frac{1}{v} \frac{\partial \boldsymbol{d} \left( \boldsymbol{p} \right) }{d \boldsymbol{p}} \\
& = \frac{1}{ { \left( v \sigma \right) }^{2} } {\left[ \frac{ \boldsymbol{p} - \boldsymbol{q}_{i} }{ \left\| \boldsymbol{p} - \boldsymbol{q}_{i} \right\| } - \frac{ \boldsymbol{p} - \boldsymbol{q}_{j} }{ \left\| \boldsymbol{p} - \boldsymbol{q}_{j} \right\| }, \ldots \right]}_{ \left(i, j \right) \in I } {\left[ \frac{ \boldsymbol{p} - \boldsymbol{q}_{i} }{ \left\| \boldsymbol{p} - \boldsymbol{q}_{i} \right\| } - \frac{ \boldsymbol{p} - \boldsymbol{q}_{j} }{ \left\| \boldsymbol{p} - \boldsymbol{q}_{j} \right\| }, \ldots \right]}_{ \left(i, j \right) \in I }^{T} \\
& = \frac{1}{ { \left( v \sigma \right) }^{2} } G {G}^{T}
\end{align*}$$

The matrix $ G \in \mathbb{R}^{d \times N} $ is basically the paris difference between unit vectors pointing from each sensor to the emitter.

The CRLB is given by the Inverse of the [Fischer Information Matrix][2]:

$$ \operatorname{CRLB} \left( \boldsymbol{p} \right) \geq {I}^{-1} \left( \boldsymbol{\theta} \right) = { \left( v \sigma \right) }^{2} { \left( G {G}^{T} \right) }^{-1} $$

## Remarks

 *  In practice the noise variance is a function of the location $ \boldsymbol{p} $ as well.




  [1]: 	https://en.wikipedia.org/wiki/Cram%C3%A9r%E2%80%93Rao_bound "Cramer Rao Lower Bound"
  [2]:  https://en.wikipedia.org/wiki/Fisher_information_matrix "Fischer Information Matrix"
