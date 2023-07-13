# Stack Exchange Question Question 87582

## Question - How to Solve a Composition of Convolutions from Regularized Least Squares Model in Frequency Domain

Assume we need to solve the model:

$$ \arg \min_{\boldsymbol{x}} \frac{1}{2} {\left\| \boldsymbol{h} \ast \boldsymbol{x} - \boldsymbol{y} \right\|}_{2}^{2} + \frac{\lambda}{2} {\left\| \boldsymbol{g} \ast \boldsymbol{x} \right\|}_{2}^{2} $$

Where $ \ast $ is the convolution operator.

How can we solve this in the frequency domain?  
Assuming different convolution models (`full`, `same` and `valid`).

The question is derived from [Solving Linear Equation of Discrete Convolution Kernels Using Black Box Model for the Convolution][1].


## Answer

In order to solve the problem let's formulate it in its matrix form:

$$ \arg \min_{\boldsymbol{x}} {\left( {H}^{T} H + \lambda {G}^{T} G \right)} \boldsymbol{x} = {H}^{T} \boldsymbol{y} $$

Where $ H $ and $ G $ are the matrix form of the convolutions.  
Since in the frequency domain all we can apply is a circular / cyclic / periodic convolution we need to create a form of the convolutions using a circulant matrix.

We can do as following:

$$ H = E {C}_{\boldsymbol{h}} P $$

Where:

 * $ P $ is the padding matrix, it adds elements to the argument.  
 * $ {C}_{\boldsymbol{h}} $ is a circulant matrix build on the convolution kernel of $ H $.
 * $ E $ is the extracting matrix. It extract what's needed from the cyclic convolution.

For instance, assume input $ \boldsymbol{x} \in \mathbb{R}^{5} $ and a kernel $ \boldsymbol{h} \in \mathbb{R}^{3} = {\left[1, 2, 3 \right]}^{T} $.  
If we want to apply a `valid` convolution using the above, then:

$$ E = \begin{bmatrix} 0 & 0 & 1 & 0 & 0 \\ 0 & 0 & 0 & 1 & 0 \\ 0 & 0 & 0 & 0 & 1 \end{bmatrix}, \; {C}_{\boldsymbol{h}} = \begin{bmatrix} 1 & 0 & 0 & 3 & 2 \\ 2 & 1 & 0 & 0 & 3 \\ 3 & 2 & 1 & 0 & 0 \\ 0 & 3 & 2 & 1 & 0 \\ 0 & 0 & 3 & 2 & 1 \end{bmatrix}, \; P = {I}_{5} $$

So for any $ \boldsymbol{x} \in \mathbb{R}^{5} $ the matrix operation $ E {C}_{\boldsymbol{h}} P \boldsymbol{x} $ will be equivalent to $ \boldsymbol{h} \ast \boldsymbol{x} $ with `valid` type of convolution. 

> * One can this in my implementation for [Replicate MATLAB's `conv2()` in Frequency Domain][2].
> * The `valid` / `same` types of convolution is not commutative. In all the above we assume the case where the kernel $ \boldsymbol{h} $ is rolling on the signal $ \boldsymbol{x} $. In order to make sense for `vlaid` we also assume the length of the signal is bigger than the length of the kernel.

Since $ {C}_{\boldsymbol{h}} $ is a circulant matrix is can be diagonalized using the DFT Matrix:

$$ {C}_{\boldsymbol{h}} = {F}^{H} {D}_{\boldsymbol{h}} F $$

Where:

 * $ {D}_{\boldsymbol{h}} $ (Diagonal Matrix) has the DFT transform of $ \boldsymbol{h} $ along its diagonal (Padded in zeros to the number of rows of $ P $).
 * $ {F} $ (Unitary Matrix) is the DFT Matrix.

Let's analyze $ {H}^{T} H $ from above:

$$ {H}^{T} H = {\left( E {C}_{\boldsymbol{h}} P \right)}^{T} \left( E {C}_{\boldsymbol{h}} P \right) = {\left( E {F}^{H} {D}_{\boldsymbol{h}} F P \right)}^{H} \left( E {F}^{H} {D}_{\boldsymbol{h}} F P \right) = {P}^{T} {F}^{H} {D}_{\boldsymbol{h}}^{H} F {E}^{T} E {F}^{H} {D}_{\boldsymbol{h}} F P $$

As can be seen above, the only way to have a closed form solution in the Frequency Domain is when $ {E}^{T} E = I $. This happens when the convolution type is `full` only.

Once that happens, we indeed get:

$$ {H}^{T} H = {P}^{T} {F}^{H} {D}_{\boldsymbol{h}}^{H} {D}_{\boldsymbol{h}} F P $$

Which basically means:

 - Pad the vector to length `length(vH) + length(vX) - 1` (With zeros).
 - Transform `vX` into DFT. Namely the 2 steps are equal to DFT with zero padding.
 - Multiply by the squared magnitude of the padded transform of the kernel.
 - Go back into the original domain (Time / Spatial, etc..).
 - Unpad the result.

For the composition, we get:

$$ {H}^{T} H + \lambda {G}^{T} G = {P}^{T} {F}^{H} \left( {D}_{\boldsymbol{h}}^{H} {D}_{\boldsymbol{h}} + \lambda {D}_{\boldsymbol{g}}^{H} {D}_{\boldsymbol{g}} \right) F P $$

So the whole process is the same.  
From here, solving the equation, in the case of `full` or `periodic` is trivial.

Unless I am missing a trick to by pass $ {E}^{T} E $, it can not be done for the other types of convolution.  
This is unfortunate because classic deconvolution (Not based on Deep Learning) is best when we assume `valid` convolution.  
Namely where `vX` is larger than the data `vY`. Then when we estimate `vX` we basically only look on its center and not edges (Estimated by far less data).



  [1]: https://dsp.stackexchange.com/questions/87542
  [2]: https://dsp.stackexchange.com/questions/74803