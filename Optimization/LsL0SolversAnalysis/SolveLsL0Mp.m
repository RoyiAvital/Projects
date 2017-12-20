function [ vX, mX ] = SolveLsL0Mp( mA, vB, paramLambda, numIterations, tolVal )
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

vActiveIdx  = false([numCols, 1]);
vR          = vB;
vX          = zeros([numCols, 1]);
vCostFun    = zeros([numIterations, 1]);
vResNorm    = zeros([numIterations, 1]);
mX          = zeros([numCols, numIterations]);

for ii = 1:numIterations
    
    vResNorm(ii) = norm(vR);
    
    % Maximum Correlation minimizes the L2 of the Error
    [~, activeIdx] = max(abs(mA.' * vR));
    % Which index minimizes the L2 Squared Error given best coefficient
    % [~, activeIdx] = min(sum((mA .* ((vR.' * mA) ./ sum(mA .^ 2)) - vR) .^ 2));
    
    vActiveIdx(activeIdx) = true([1, 1]);
    
    % Solve the equation with respect to the residual
    coeffVal = (mA(:, activeIdx).' * vR) ./ (mA(:, activeIdx).' * mA(:, activeIdx));
    % Update the coefficient according to the current solution
    vX(activeIdx) = vX(activeIdx) + coeffVal;
    
    % Remove the current reduction in error
    vR = vR - (mA(:, activeIdx) * coeffVal);
    % Equivalent
    % vR = vB - (mA * vX);
    
    mX(:, ii) = vX;
    vCostFun(ii) = (0.5 * sum(((mA * vX) - vB) .^ 2)) + (paramLambda * sum(abs(vX) > tolVal));
    
end

[~, minCostValIdx] = min(vCostFun);
vX = mX(:, minCostValIdx);


end

