function [ vX, mX ] = SolveLsL1Irls( mA, vB, paramLambda, numIterations )
% ----------------------------------------------------------------------------------------------- %
%[ vX, mX ] = SolveLsL1Irls( mA, vB, paramLambda, numIterations )
% Solve L1 Regularized Least Squares Using Iteratively Reweighted Least Squares (IRLS) Method.
% Input:
%   - mA                -   Input Matrix.
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
%   1.  D
% Known Issues:
%   1.  D
% TODO:
%   1.  P
% Release Notes:
%   -   1.1.000     13/05/2015
%       *   Fixed the 2 factor bug.
%       *   Updating only the Diagonal of mAD for more efficiency.
%   -   1.0.000     23/08/2017
%       *   First release version.
% ----------------------------------------------------------------------------------------------- %

mAA = mA.' * mA;
mAD = mAA;
vAb = mA.' * vB;
vX  = pinv(mA) * vB; %<! Dealing with "Fat Matrix"

thrBaseVal = 1e-6;
thrMinVal  = 1e-12;

mX = zeros([size(vX, 1), numIterations]);
mX(:, 1) = vX;

thrVal = thrBaseVal;
vW = ones([size(vX, 1), 1]);

vDiagIdx = linspace(1, size(vX, 1) * size(vX, 1), size(vX, 1));
vDiagIdx = vDiagIdx(:);

for ii = 1:numIterations
    
    mAD(vDiagIdx) = mAA(vDiagIdx) + (paramLambda * vW);
    vX = mAD \ vAb;
    
    vW = 1 ./ (abs(vX) + thrVal);
    
    mX(:, ii) = vX;
    thrVal = thrVal / 2;
    thrVal = max(thrVal, thrMinVal);
    
end


end

