function [ vX, mX ] = SolveLsL1Irls( mA, vB, paramLambda, numIterations )
% ----------------------------------------------------------------------------------------------- %
%[ vX, mX ] = SolveLsL1Irls( mA, vB, paramLambda, numIterations )
% Solve L1 Regularized Least Squares Using Iteratively Reweighted Least Squares (IRLS) Method.
% Input:
%   - mA                -   Input Matirx.
%                           The model matrix.
%                           Structure: Matrix (m X n).
%                           Type: 'Single' / 'Double'.
%                           Range: (-inf, inf).
%   - vB                -   input Vector.
%                           The model known data.
%                           Structure: Vector (m X 1).
%                           Type: 'Single' / 'Double'.
%                           Range: (-inf, inf).
%   - paramLambda       -   Parameter Lambda.
%                           The L1 Regularization parameter.
%                           Structure: Scalar.
%                           Type: 'Single' / 'Double'.
%                           Range: (0, inf).
%   - numIterations     -   Number of Iterations.
%                           Number of iterations of the algorithm.
%                           Structure: Scalar.
%                           Type: 'Single' / 'Double'.
%                           Range {1, 2, ...}.
% Output:
%   - vX                -   Output Vector.
%                           Structure: Vector (n X 1).
%                           Type: 'Single' / 'Double'.
%                           Range: (-inf, inf).
% References
%   2.  IRLS Wikipedia - https://en.wikipedia.org/wiki/Iteratively_reweighted_least_squares.
%   1.  L1 IRLS (My Derivation) - https://stats.stackexchange.com/a/299381/6244.
% Remarks:
%   1.  Doesn't work well for large values of paramLambda.
% Known Issues:
%   1.  Doesn't converge to a "good" solution.
% TODO:
%   1.  P
% Release Notes:
%   -   1.0.000     23/08/2017
%       *   First realease version.
% ----------------------------------------------------------------------------------------------- %

mAA = mA.' * mA;
vAb = mA.' * vB;
vX  = pinv(mA) * vB; %<! Dealing with "Fat Matrix"

thrBaseVal = 1e-6;
thrMinVal  = 1e-12;

mX = zeros([size(vX, 1), numIterations]);
mX(:, 1) = vX;

thrVal = thrBaseVal;
mW = eye(size(vX, 1));

for ii = 1:numIterations
    
    vX = (mAA + (2 * paramLambda * mW)) \ vAb;
    
    mW = diag(1 ./ (abs(vX) + thrVal));
    
    mX(:, ii) = vX;
    thrVal = thrVal / 2;
    thrVal = max(thrVal, thrMinVal);
    
end


end

