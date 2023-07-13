This is not an easy question, not because ideas are complex but since the details matter a lot.

As I wrote in my answer in the question [Solving Inverse Problem Using Black Box Implementation of the Convolution Kernel][1], the issue is defining how to handle the boundaries (See [The Different Solutions for Filter Coefficients Estimation for Periodic Convolution and Full Convolution][2]). This is super important for images.

### Matrix Case

In all cases, for this equation:

$$ {\left( {H}^{T} H + \lambda {G}^{T} G \right)} \boldsymbol{x} = \boldsymbol{y} \Leftrightarrow A \boldsymbol{x} = \boldsymbol{y} $$

Can be solved using sparse solvers taking advantage of the sparsity pattern of the matrices and the matrix being symmetric positive definite.  

This will work for any variant of the convolution.  
Though we have 2 options, direct solvers and iterative ones.  

The iterative ones allow us taking advantage of the convolution black box as they usually require applying the operator $ A $. Many iterative solvers require $ {A}^{T} $ but if they do, it means they don't model it as a symmetric operator. Hence it means they are not optimal. So stick with those assuming symmetric operator, hence require $ A $ only.

Now the question, how to implement $ A $ efficiently?

### Periodic / Cyclic Boundary Conditions

If we assume that the $ \ast $ (Convolution) is applied using cyclic boundary convolution we can apply everything in frequency domain.  
I wrote many answer about it, you may look at them:

 * [Applying Convolution in Frequency Domain by Element Wise Multiplication on Time Domain][3].
 * [Dealing with the Cyclic Boundary Conditions of Frequency Domain Convolution in Order to Apply Linear Convolution][4].
 * [Circular Convolution Matrix of $ {H}^{T} H $][5].
 * [The Matrix Form of a 2D Circular Convolution][6].
 * [Least Squares Solution Using the DFT vs Wiener Hopf Equations][7].
 * [Estimate Filter Coefficients from the Result of Linear Convolution with a Known Signal][8].
 * [Strategy / Method for Implementation of the Fastest 1D Linear Convolution / Correlation][9].
 * [Show Equivalence Between Multiplication in Time Domain to Convolution in Frequency Domain][10].

### The Convolution of Type `full` / `same` / `valid`

For the other convolution types (I Use the lingo of MATLAB's convolution functions, such as [`conv2()`][11]) things are trickier.

I created the matrix $ {H}^{T} H $ for the different convolution shapes from the same 3 coefficients kernel applied on 5 samples signal:

[![enter image description here][12]][12]

As one can see, while all cases generates a `5x5` matrix, yet only the `full` case generates a [Toeplitz Matrix][13].  
Hence, only this case can be represented by a convolution kernel with some convolution type.

#### General Solution

We can always apply each matrix by itself. Something like:

$$ {\left( {H}^{T} H + \lambda {G}^{T} G \right)} \boldsymbol{x} = {H}^{T} \left( H \boldsymbol{x} \right) + \lambda {G}^{H} \left( G \boldsymbol{x} \right) $$

So $ H \boldsymbol{x} $ / $ G \boldsymbol{x} $ is just applying the convolution, what about $ {H}^{T} \boldsymbol{z} $ / $ {G}^{T} \boldsymbol{z} $?
Well, we can have a recipe per case, where the adjoint is done by the flipped kernel:

 * Convolution Type `full` -> Correlation type `valid`. 
 * Convolution Type `same` -> Correlation type `same`. 
 * Convolution Type `valid` -> Correlation type `full`. 

Now we can apply all of the above one by one.

#### Optimized Solution for `full`

As we can see, in case our convolution is of type `full` the normal matrix $ {H}^{T} H $ is a Toeplitz Matrix, so we can generate the composite operation using a single kernel convolution.

### Results

I implemented all cases above and using Conjugate Gradient I got the following results for 1D:

[![enter image description here][14]][14]

As can be seen, we can get the same results using the iterative method.  
This method only applies convolution operations and doesn't require the matrix form.

### Some Other Resources

 * [Estimating Convolution Input Under the Assumption of Sparsity and Constant Non Zero Values Using Compressive Sensing Approach][15].


  [1]: https://dsp.stackexchange.com/questions/87500
  [2]: https://dsp.stackexchange.com/questions/87406
  [3]: https://dsp.stackexchange.com/questions/80687
  [4]: https://dsp.stackexchange.com/questions/78918
  [5]: https://dsp.stackexchange.com/questions/56709
  [6]: https://dsp.stackexchange.com/questions/81949
  [7]: https://dsp.stackexchange.com/questions/87326
  [8]: https://dsp.stackexchange.com/questions/45879
  [9]: https://dsp.stackexchange.com/questions/52760
  [10]: https://dsp.stackexchange.com/questions/64870
  [11]: https://www.mathworks.com/help/matlab/ref/conv2.html
  [12]: https://i.stack.imgur.com/nztwj.png
  [13]: https://en.wikipedia.org/wiki/Toeplitz_matrix
  [14]: https://i.stack.imgur.com/Z8BtP.png
  [15]: https://dsp.stackexchange.com/questions/66778