# Convex Set Projection
The orthogonal projection onto various convex sets.

## Cicruclant Matrix
Project into tyhe closest Circulant Matrix set given by $\mathcal{C}_{n}$:

$$ \arg \min_{X \in \mathcal{C}_{n}} \frac{1}{2} {\left\| X - Y \right\|}_{F} $$

See my solutions at [Projection onto the Set of Circulant Matrices](https://math.stackexchange.com/questions/2778195).  

**MATLAB Note**  
If `mF` is the DFT Matrix and `mF'` is the Inverse DFT Matrix then:
 *  Multiplication from Left - `mF * mY` is equivalent to `fft(mY, [], 1)`.  
    Namley applying DFT on the columns of `mY`.  
    This also folds for inverse DFT, namely `mF' * mY` is equivalent of `ifft(mY, [], 1)`.
 *  Multiplication from Right - `mY * mF` is equivalent to `fft(mY, [], 2)`.  
    Namley applying DFT on the rows of `mY`.  
    This also folds for inverse DFT, namely `mY * mF'` is equivalent of `ifft(mY, [], 2)`.
 *  Compositions - `mF * mY * mF'` is equivelent of `ifft(fft(mY, [], 1), [], 2)` and `mF' * mY * mF` is equivalent of `fft(ifft(mY, [], 1), [], 2)`.  
    The order could be switched (Doing first the left operator then the right).
 *  Pay attention that in MATLAB `fft()` isn't normalized by `sqrt(n)` and `ifft()` is normalized by `n`.  
    Hence the compositions above will work yet the first 2 points will require normalization.

## TO DO List
 *  Semi Positive Definite Matrix.
 *  Solution Set of a Linear System of equations.
