# Estimate Harmonic Signals Frequency

Estimating the parameters of _Harmonic Signal_.

## Problem 

Given a signal in the form $ y \left[ n \right] = A \sin \left[ 2 \pi f n + \phi \right] + w \left[ n \right] $ where:

 * $ f $ is the normalized frequency, namely $ f = \frac{ {f}_{0} }{ {f}_{s} } $ where $ {f}_{0} $ is the signal frequency and $ {f}_{s} $ is the sampling frequency.   
 * The signal is sampled above [_Nyquist Frequency_](https://en.wikipedia.org/wiki/Nyquist_frequency) hence $ 0 < f < 0.5 $.
 * The phase is given in the range $ \phi \in \left[0, 2 \pi \right] $.
 * The noise $ w \left[ n \right] $ is [_Additive White Gaussian Noise_](https://en.wikipedia.org/wiki/Additive_white_Gaussian_noise) (AWGN) with variance of $ {\sigma}^{2} $.

The parameter estimation problem is to estimate the parameters $ \boldsymbol{\theta} = {\left[ A, f, \phi \right]}^{T} $.

## The Cramer Rao Lower Bound (CRLB)

The [Fisher Information Matrix](https://en.wikipedia.org/wiki/Fisher_information) is given by[^001]:

$$ I \left[ \boldsymbol{\theta} \right]_{i, j} = \frac{1}{ {\sigma}^{2} } \left[ \frac{\partial \boldsymbol{x}_{\boldsymbol{\theta}} }{\partial \boldsymbol{\theta}_{i}} \right] {\left[ \frac{\partial \boldsymbol{x}_{\boldsymbol{\theta}} }{\partial \boldsymbol{\theta}_{j}} \right]}^{T} = \frac{1}{ {\sigma}^{2} } \sum_{n = 0}^{N - 1} \frac{\partial x \left[ n ; \boldsymbol{\theta} \right] }{\partial \boldsymbol{\theta}_{i}} \frac{\partial x \left[ n ; \boldsymbol{\theta} \right] }{\partial \boldsymbol{\theta}_{j}} $$

Which yields:

$$ I \left[ \boldsymbol{\theta} \right] \approx \begin{bmatrix}
\frac{N}{2 {\sigma}^{2}} & 0 & 0 \\ 
0 & \frac{ {A}^{2} }{2 {\sigma}^{2}} \sum_{n = 0}^{N - 1} {n}^{2} & \frac{ {A}^{2} }{2 {\sigma}^{2}} \sum_{n = 0}^{N - 1} n \\ 
0 & \frac{ {A}^{2} }{2 {\sigma}^{2}} \sum_{n = 0}^{N - 1} n & \frac{N {A}^{2} }{2 {\sigma}^{2}}
\end{bmatrix} $$

From the block diagonal form of the matrix one can conclude that the estimation of the amplitude is independent of the estimation of the frequency and phase which are dependent on each other.

The specific [Cramer Rao Lower Bound](https://en.wikipedia.org/wiki/Cram%C3%A9r%E2%80%93Rao_bound) (CRLB) lower bound for inverse is given by the inverse matrix which yields:

 * $ \operatorname{var} \left( \hat{A} \right) \geq \frac{2 {\sigma}^{2}}{N} $.
 * $ \operatorname{var} \left( \hat{f} \right) \geq \frac{12}{ {\left( 2 \pi \right)}^{2} \eta N \left( {N}^{2} - 1 \right) } $.
 * $ \operatorname{var} \left( \hat{\phi} \right) \geq \frac{2 \left( 2 N - 1 \right)}{ \eta N \left( N + 1 \right) } $.

Where $ \eta = \frac{ {A}^{2} }{ 2 {\sigma}^{2} } $ is the SNR.  

One can observe that the the ability to estimate the frequency is inversely proportional to $ {N}^{3} $ which means the longer we observe the signal the better we can estimate its frequency.
This is aligns with the intuition that the sampling interval ($ \frac{N}{ {T}_{s} } $) determines the resolution of the DFT.

## Estimators

### Efficient Maximum Likelihood Estimator  

The [Maximum Likelihood Estimator](https://en.wikipedia.org/wiki/Maximum_likelihood_estimation) for the problem is given by[^002]:

$$ \arg \min_{A, f, \phi} \frac{1}{2} {\left\| A \boldsymbol{s} - \boldsymbol{y} \right\|}_{2}^{2} $$  

Where $ \boldsymbol{s} = \begin{bmatrix} A \sin \left[ 2 \pi f 0 + \phi \right] \\ A \sin \left[ 2 \pi f 1 + \phi \right] \\ A \sin \left[ 2 \pi f 2 + \phi \right] \\ \vdots \\ A \sin \left[ 2 \pi f (N - 1) + \phi \right] \end{bmatrix} $.  
One should pay attention that the we can always choose to start with $ n = 0 $ as the phase can compensate for this.

This section describes an efficient[^003] method to implement the estimator.  
The idea is transform this 3D problem into a 1D problem where we can employ many methods to solve.

By using a [trigonometric identity](https://en.wikipedia.org/wiki/List_of_trigonometric_identities) we can _linearize_ the estimation of the phase:

$$ A \sin \left[ 2 \pi f n + \phi \right] = A \cos \left( \phi \right) \sin \left[ 2 \pi f n \right] + A \sin \left( \phi \right) \cos \left[ 2 \pi f n \right] $$

Let's define the matrix $ \boldsymbol{X} \left( f \right) $ and the vector $ \boldsymbol{a} \left( \phi \right) $:

$$ \boldsymbol{X} \left( f \right) = \begin{bmatrix} \sin \left[ 2 \pi f 0 \right] & \cos \left[ 2 \pi f 0 \right] \\ \sin \left[ 2 \pi f 1 \right] & \cos \left[ 2 \pi f 1 \right] \\ \sin \left[ 2 \pi f 2 \right] & \cos \left[ 2 \pi f 2 \right] \\ \vdots & \vdots \\ \sin \left[ 2 \pi f \left( N - 1 \right) \right] & \cos \left[ 2 \pi f \left( N - 1 \right) \right] \end{bmatrix}, \; \boldsymbol{a} \left( A, \phi \right) = \begin{bmatrix} A \cos \left( \phi \right) \\ A \sin \left( \phi \right) \end{bmatrix} $$

Then the minimization problem becomes:

$$ \arg \min_{A, f, \phi} \frac{1}{2} {\left\| \boldsymbol{X} \boldsymbol{a} - \boldsymbol{y} \right\|}_{2}^{2} $$  

Yet, the estimator for $ \boldsymbol{a} $ is given by $ \hat{\boldsymbol{a}} = \boldsymbol{X}^{\dagger} \boldsymbol{y} $.  
Integrating this into the minimization problem yields[^004]:

$$
\begin{align*}
\arg \min_{A, f, \phi} \frac{1}{2} {\left\| \boldsymbol{X} \boldsymbol{a} - \boldsymbol{y} \right\|}_{2}^{2} & = \arg \min_{f} \frac{1}{2} {\left\| \boldsymbol{X} \boldsymbol{X}^{\dagger} \boldsymbol{y} - \boldsymbol{y} \right\|}_{2}^{2} \\
& = \arg \min_{f} \frac{1}{2} {\left\| \boldsymbol{X} {\left( \boldsymbol{X}^{T} \boldsymbol{X} \right)}^{-1} \boldsymbol{X}^{T} \boldsymbol{y} - \boldsymbol{y} \right\|}_{2}^{2} \\
& = \arg \min_{f} \frac{1}{2} {\left\| \left( \boldsymbol{X} {\left( \boldsymbol{X}^{T} \boldsymbol{X} \right)}^{-1} \boldsymbol{X}^{T} - \boldsymbol{I} \right) \boldsymbol{y} \right\|}_{2}^{2}
\end{align*}
$$  

> TODO: This can be simplified into something like $ \boldsymbol{y}^{T} \boldsymbol{X} {\left( \boldsymbol{X}^{T} \boldsymbol{X} \right)}^{-1} \boldsymbol{X}^{T} \boldsymbol{y} $.

In our case, the matrix $ \left(\boldsymbol{X}^{T} \boldsymbol{X}\right) \in \mathbb{R}^{2 \times 2} $, which has a simple closed form.  
This makes the search for the optimal frequency quite simple even over a grid.  
A good initialization could be made by the frequency bin which maximizes the the magnitude of the DFT of the data.   

Once the frequency is found, the vector $ \boldsymbol{a} $ can be estimated by trigonometric identity by using a polar coordinates.

An alternative to a grid search can be with the initialization by the DFT and using _Brent_ method, which require no derivative, to find the minimum.  
In case of a need for the derivative, it is given by (Using the chain rule):

$$ \frac{\mathrm{d} g }{\mathrm{d} f} = \left \langle \frac{\mathrm{d} g }{\mathrm{d} \boldsymbol{X}}, \frac{\mathrm{d} \boldsymbol{X} }{\mathrm{d} f} \right \rangle $$

Calculating $ \frac{\mathrm{d} \boldsymbol{X} }{\mathrm{d} f} $ is straight forward as:

$$ \frac{\mathrm{d} \boldsymbol{X} }{\mathrm{d} f} = \begin{bmatrix} 2 \pi 0 \cos \left[ 2 \pi f 0 \right] & -2 \pi 0 \sin \left[ 2 \pi f 0 \right] \\ 2 \pi 1 \cos \left[ 2 \pi f 1 \right] & -2 \pi 1 \sin \left[ 2 \pi f 1 \right] \\ 2 \pi 2 \cos \left[ 2 \pi f 2 \right] & -2 \pi 2 \sin \left[ 2 \pi f 2 \right] \\ \vdots & \vdots \\ 2 \pi \left( N - 1 \right) \cos \left[ 2 \pi f \left( N - 1 \right) \right] & -2 \pi \left( N - 1 \right) \sin \left[ 2 \pi f \left( N - 1 \right) \right] \end{bmatrix} $$

The other derivative:

$$ \frac{\mathrm{d} g }{\mathrm{d} \boldsymbol{X}} = \boldsymbol{T}_{3} \boldsymbol{T}_{2} - \boldsymbol{T}_{1} \boldsymbol{T}_{4} + \boldsymbol{X} \boldsymbol{T}_{0} \boldsymbol{X}^{T} \boldsymbol{T}_{3} \boldsymbol{T}_{2}  + \boldsymbol{y} \boldsymbol{T}_{4} $$

Where

 * $ \boldsymbol{T}_{0} = {\left( \boldsymbol{X}^{T} \boldsymbol{X} \right)}^{-1} $.
 * $ \boldsymbol{T}_{1} = \boldsymbol{X} \boldsymbol{T}_{0} \boldsymbol{X}^{T} \boldsymbol{y} $.
 * $ \boldsymbol{T}_{2} = \boldsymbol{y}^{T} \boldsymbol{X} \boldsymbol{T}_{0} $.
 * $ \boldsymbol{T}_{3} = \boldsymbol{T}_{1} - \boldsymbol{y} $.
 * $ \boldsymbol{T}_{4} = \left( \boldsymbol{T}_{2} \boldsymbol{X}^{T} - \boldsymbol{y}^{T} \right) \boldsymbol{X} \boldsymbol{T}_{0} $.


**Remark**: In case of low SNR the DFT maximum might not be in the correct interval. Hence doing a grid search and then finding the optimum will be better.


#### Simpler Notation

The minimization argument can be farther simplified:

$$
\begin{aligned}
\arg \min_{f} \frac{1}{2} {\left\| \boldsymbol{X} {\left( \boldsymbol{X}^{T} \boldsymbol{X} \right)}^{-1} \boldsymbol{X}^{T} \boldsymbol{y} - \boldsymbol{y} \right\|}_{2}^{2} & = \arg \min_{f} {\left( \boldsymbol{X} {\left( \boldsymbol{X}^{T} \boldsymbol{X} \right)}^{-1} \boldsymbol{X}^{T} \boldsymbol{y} - \boldsymbol{y} \right)}^{T} \left( \boldsymbol{X} {\left( \boldsymbol{X}^{T} \boldsymbol{X} \right)}^{-1} \boldsymbol{X}^{T} \boldsymbol{y} - \boldsymbol{y} \right) && \text{} \\
& = \arg \min_{f} \left( \boldsymbol{y}^{T} \boldsymbol{X} {\left( \boldsymbol{X}^{T} \boldsymbol{X} \right)}^{-1} \boldsymbol{X}^{T} - \boldsymbol{y}^{T} \right) \left( \boldsymbol{X} {\left( \boldsymbol{X}^{T} \boldsymbol{X} \right)}^{-1} \boldsymbol{X}^{T} \boldsymbol{y} - \boldsymbol{y} \right) && \text{} \\
& = \arg \min_{f} \boldsymbol{y}^{T} \boldsymbol{X} {\left( \boldsymbol{X}^{T} \boldsymbol{X} \right)}^{-1} \boldsymbol{X}^{T} \boldsymbol{X} {\left( \boldsymbol{X}^{T} \boldsymbol{X} \right)}^{-1} \boldsymbol{X}^{T} \boldsymbol{y} - 2 \boldsymbol{y}^{T} \boldsymbol{X} {\left( \boldsymbol{X}^{T} \boldsymbol{X} \right)}^{-1} \boldsymbol{X}^{T} \boldsymbol{y} + \boldsymbol{y}^{T} \boldsymbol{y} && \text{} \\
& = \arg \min_{f} \boldsymbol{y}^{T} \boldsymbol{X} {\left( \boldsymbol{X}^{T} \boldsymbol{X} \right)}^{-1} \boldsymbol{X}^{T} \boldsymbol{y} - 2 \boldsymbol{y}^{T} \boldsymbol{X} {\left( \boldsymbol{X}^{T} \boldsymbol{X} \right)}^{-1} \boldsymbol{X}^{T} \boldsymbol{y} && \text{} \\
& = \arg \max_{f} \boldsymbol{y}^{T} \boldsymbol{X} {\left( \boldsymbol{X}^{T} \boldsymbol{X} \right)}^{-1} \boldsymbol{X}^{T} \boldsymbol{y} && \text{} \\
&& \blacksquare
\end{aligned}
$$

With the derivative:

$$ \frac{\mathrm{d} g }{\mathrm{d} \boldsymbol{X}} = 2 \left( \boldsymbol{y} \boldsymbol{t}_{1} - \boldsymbol{X} \boldsymbol{T}_{0} \boldsymbol{X}^{T} \boldsymbol{y} \boldsymbol{t}_{1} \right) $$

Where

 * $ \boldsymbol{T}_{0} = {\left( \boldsymbol{X}^{T} \boldsymbol{X} \right)}^{-1} $.
 * $ \boldsymbol{t}_{1} = \boldsymbol{y}^{T} \boldsymbol{X} \boldsymbol{T}_{0} $.



## To Do List

 * 	Add monotonic smoothness to smooth the variance of the estimators.
 
## Resources

 *	[Isotonic Regression Code MATLAB](https://www.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/47196/versions/1/previews/improve_JP/toolbox_imp_JP/lsqisotonic.m/index.html).
 *	[Pool Adjacent Violators Python](https://gist.github.com/fabianp/3081831).
 *	[Isotonic Regression Code MATLAB](https://www.mathworks.com/matlabcentral/fileexchange/64933).
 *	[Isotone Optimization in R: Pool Adjacent Violators Algorithm (PAVA) and Active Set Methods](https://cran.r-project.org/web/packages/isotone/vignettes/isotone.pdf).
 *	[Isotonic Regression](https://en.wikipedia.org/wiki/Isotonic_regression).
 *	[Stat 8054 Lecture Notes: Isotonic Regression](https://www.stat.umn.edu/geyer/8054/notes/isotonic.pdf).
 *	[A smoothed Monotonic Regression via L2 Regularization](https://link.springer.com/article/10.1007/s10115-018-1201-2).
 *  [Julien Arzi - Frequency Estimation](http://www.tsdconseil.fr/log/scriptscilab/festim/index-en.html).

## License
Free for non commercial use.


  [^001]: See Steven Kay - Fundamentals of Statistical Signal Processing, Volume I: Estimation Theory (Example `3.14` page 56).
  [^002]: Since the noise is AWGN the maximum likelihood estimator is given by the Least Squares estimator.
  [^003]: A more generalized form of it (For sum of sines) was derived in a project in the RADAR world.
  [^004]: The equality in the derivation assumes optimal estimation of the phase. As there is no efficient estimator for the phase this is only an approximation in practice.

