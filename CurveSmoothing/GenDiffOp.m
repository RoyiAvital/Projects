function [ mD ] = GenDiffOp( diffPow, numSamples )
% ----------------------------------------------------------------------------------------------- %
% [ mD ] = GenDiffOp( diffPow, numSamples )
%   Generates the m -th derivative finite differences operator in a from of
%   a matrix. The matrix is designed to operate on vectors. The algorithm
%   uses the Central Finite Difference Coefficients. 
% Input:
%   - diffPow       -   The Derivative Order.
%                       Sets the order of the derivative.
%                       Structure: Scalar.
%                       Type: 'Single' / 'Double'.
%                       Range: {1, 2, ..., 6}.
%   - numSamples    -   Number of Samples.
%                       The number of samples of the vector the operator
%                       should be designed to
%                       Structure: Scalar.
%                       Type: 'Single' / 'Double'.
%                       Range: {1, 2, ...}.
% Output:
%   - mD            -   Finite Differences Operator.
%                       The finite differences operator in a form of a
%                       sparse matrix.
%                       Structure: Matrix (Sparse).
%                       Type: 'Single' / 'Double'.
%                       Range: (-inf, inf).
% References
%   1.  https://en.wikipedia.org/wiki/Finite_difference_coefficient.
% Remarks:
%   1.  Using the Central variant of the coefficients.
% TODO:
%   1.  U.
% Release Notes:
%   -   1.0.000     09/01/2026  Royi Avital
%       *   First release version.
% ----------------------------------------------------------------------------------------------- %

arguments(Input)
    diffPow (1, 1) {mustBeNumeric, mustBeFinite, mustBeInteger, mustBePositive}
    numSamples (1, 1) {mustBeNumeric, mustBeFinite, mustBeInteger, mustBePositive}
end

arguments(Output)
    mD (:, :) {mustBeNumeric, mustBeFinite, mustBeSparse}
end

vC = GenerateCentralDiffsCoefs(diffPow);
numCoeff = length(vC);
coeffRadius = ((numCoeff - 1) / 2);

% Generate the 1st Derivative Diff Matrix
mD = spdiags(vC, -coeffRadius:coeffRadius, numSamples, numSamples);
mD = mD((1 + coeffRadius):(end - coeffRadius), :);

% for ii = 2:diffPow
%     mD = mD' * mD;
% end


end


function [ vC ] = GenerateCentralDiffsCoefs( diffPow )

arguments(Input)
    diffPow (1, 1) {mustBeNumeric, mustBeFinite, mustBeInteger, mustBePositive, mustBeInRange(diffPow, 1, 6)}
end

arguments(Output)
    vC (1, :) {mustBeNumeric, mustBeFinite}
end

switch(diffPow)
    case(1)
        vC = [-0.5, 0, 0.5];
    case(2)
        vC = [1, -2, 1];
    case(3)
        vC = [-0.5, 1, 0, -1, 0.5];
    case(4)
        vC = [1, -4, 6, 4, 1];
    case(5)
        vC = [-0.5, 2, -2.5, 0, 2.5, -2, 0.5];
    case(6)
        vC = [1 -6, 15, -20, 15, -6, 1];
end


end


function [ vC ] = GenerateCentralDiffsCoefs_( diffPow, accLvl, zeroThr )
% ----------------------------------------------------------------------------------------------- %
% [ vC ] = GenerateCentralDiffsCoefs_( diffPow, accLvl, zeroThr )
%   Calculates the central finite differences coefficients given the order
%   of the derivative and the accuracy level.
%
% Input:
%   - diffPow       -   The Derivative Order.
%                       Sets the order of the derivative.
%                       Structure: Scalar.
%                       Type: 'Single' / 'Double'.
%                       Range: {1, 2, ..., 8}.
%   - accLvl        -   The Accuracy Level.
%                       Sets the accuracy level of the filter.
%                       Structure: Scalar.
%                       Type: 'Single' / 'Double'.
%                       Range: {2, 4, ..., 16}.
%   - zeroThr       -   Zeroing Threshold.
%                       Sets absolute values which are smaller than
%                       threshold to zero.
%                       Structure: Scalar.
%                       Type: 'Single' / 'Double'.
%                       Range: (0, inf).
% Output:
%   - vC            -   Coefficients Vector.
%                       The central differences coefficients.
%                       Structure: Vector (Column).
%                       Type: 'Single' / 'Double'.
%                       Range: (-inf, inf).
% References
%   1.  https://en.wikipedia.org/wiki/Finite_difference_coefficient.
%       See the section about the calculation of the coefficients.
% Remarks:
%   1.  Assumes the derivative power is at most 8 and the accuracy is at
%       most 16.
%   2.  The `accLvl` must be even number.
%   3.  There are 2p + 1 or 2 * floor((m + 1) / 2) - 1 + n central
%       coefficients where m is the degree and n is the accuracy.
%       The coefficients are given by the solution of the linear system A * c = b :
%       [   1        1      .. 1 ..   1      1    ][a_-p]   [ 0 ]
%       [  -p      -p+1     .. 0 ..  p-1     p    ][  : ]   [ : ]
%       [ (-p)^2  (-p+1)^2  .. 0 .. (p-1)^2  p^2  ][  : ] = [ m!]
%       [   :        :         :      :       :   ][  : ]   [ : ]
%       [ (-p)^2p (-p+1)^2p .. 0 .. (p-1)^2p p^2p ][a_+p]   [ 0 ]
%
%       where the only non zero value on the R.H.S is in the (m+1)-th row.
%
% TODO:
%   1.  U.
% Release Notes:
%   -   1.0.000     09/01/2026  Royi Avital
%       *   First release version.
% ----------------------------------------------------------------------------------------------- %

arguments(Input)
    diffPow (1, 1) {mustBeNumeric, mustBeFinite, mustBeInteger, mustBePositive}
    accLvl (1, 1) {mustBeNumeric, mustBeFinite, mustBeInteger, mustBePositive} = 2
    zeroThr (1, 1) {mustBeNumeric, mustBeFinite, mustBeNonnegative} = sqrt(eps())
end

arguments(Output)
    vC (:, 1) {mustBeNumeric, mustBeFinite}
end

numCoeff =  2 * floor((diffPow + 1) / 2) - 1 + accLvl;
valP = (numCoeff - 1) / 2;


% Solve system mA * vA = vB
mA = power(-valP:valP, (0:(2 * valP))'); %<! Vandermonde matrix

vB     = zeros(2 * valP + 1, 1);
ii     = diffPow + 1;
vB(ii) = factorial(diffPow);
 
vC = mA \ vB;
% Round small numbers to zero
vC(abs(vC) < zeroThr) = 0;
vC = vC - sum(vC); %<! Ensure zero sum

end

