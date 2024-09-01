function [ vX, mX ] = SolveLsL1Huber( mA, vB, paramLambda, numIterations )
% ----------------------------------------------------------------------------------------------- %
%[ vX, mX ] = SolveLsL1Huber( mA, vB, paramLambda, numIterations )
% Solve L1 Regularized Least Squares Using Smoothing (Huber Loss) Method.
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
%   1.  Huber Loss Wikipedia - https://en.wikipedia.org/wiki/Huber_loss.
% Remarks:
%   1.  As the smoothness term approaches zero the Huber Loss better
%       approximate the L1 Norm. Yet the lower the value the harder to
%       solve hence "Warm Start" is used.
% Known Issues:
%   1.  D.
% TODO:
%   1.  Add line search (Backtracking).
% Release Notes:
%   -   1.0.000     25/08/2017
%       *   First release version.
% ----------------------------------------------------------------------------------------------- %

mAA = mA.' * mA;
vAb = mA.' * vB;
vX  = pinv(mA) * vB; %<! Dealing with "Fat Matrix"

lipConst = norm(mA, 2) ^ 2;%<! Lipschitz Constant;

paramMuBase     = 0.005; %<! Smoothness term in Huber Loss
stepSizeBase    = 1 / lipConst;

mX(:, 1) = vX;

for ii = 2:numIterations
    
    paramMu     = paramMuBase / log2(ii);
    stepSize    = stepSizeBase / log2(ii);
    
    vG = (mAA * vX) - vAb + (paramLambda * HuberLossGrad(vX, paramMu));
    vX = vX - (stepSize * vG);
    
    mX(:, ii) = vX;
    
end


end


function [ vG ] = HuberLossGrad( vX, paramMu )

vG = ((abs(vX) <= paramMu) .* (vX ./ paramMu)) + ((abs(vX) > paramMu) .* sign(vX));


end

