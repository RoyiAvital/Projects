function [ vX, mX ] = SolveLsL0Omp( mA, vB, paramLambda, numIterations, tolVal )
% ----------------------------------------------------------------------------------------------- %
%[ vX, mX ] = SolveLsL0Mp( mA, vB, paramLambda, numIterations )
% Solve L0 Regularized Least Squares Using Matching Pursuit (MP) Method.
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
%                           The L0 Regularization parameter.
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
%   1.  Wikipedia MP - https://en.wikipedia.org/wiki/Matching_pursuit.
% Remarks:
%   1.  U.
% Known Issues:
%   1.  A
% TODO:
%   1.  Pre Process 'mA' by normalizing its columns.
% Release Notes:
%   -   1.0.000     04/12/2017
%       *   First realease version.
% ----------------------------------------------------------------------------------------------- %

numRows = size(mA, 1);
numCols = size(mA, 2);

% Since all atoms are used and the solutions is optimal, number of
% iterations can not exceed number of atoms as the resdiaul is at its
% minimum.
numIterations = min(numIterations, numCols);

vActiveIdx  = false([numCols, 1]);
vR          = vB;
vX          = zeros([numCols, 1]);
vCostFun    = zeros([numIterations, 1]);
vResNorm    = zeros([numIterations, 1]);
mX          = zeros([numCols, numIterations]);

mAA = mA ./ sum(mA .* mA);

for ii = 1:numIterations
    
    vResNorm(ii) = norm(vR);
    
    [~, activeIdx] = max(abs(mA.' * vR));
    vActiveIdx(activeIdx) = true([1, 1]);
    
    vX(vActiveIdx) = mA(:, vActiveIdx) \ vB;
    vR = vB - (mA(:, vActiveIdx) * vX(vActiveIdx));
    
    mX(:, ii) = vX;
    vCostFun(ii) = (0.5 * sum(((mA * vX) - vB) .^ 2)) + (paramLambda * sum(abs(vX) > tolVal));
    
end

[~, minCostValIdx] = min(vCostFun);
vX = mX(:, minCostValIdx);


end

