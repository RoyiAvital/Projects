function [ vX, mX ] = SolveLsL1ProxAccel( mA, vB, lambdaFctr, numIterations )
% ----------------------------------------------------------------------------------------------- %
%[ vF ] = BayesianDft( vX, numFreqBins, varX, varN, numIterations )
% High Resolution DFT using Bayesian Estimation of the DFT coefficients.
% Input:
%   - vX                -   Input Vector.
%                           Structure: Vector (Column Vector).
%                           Type: 'Single' / 'Double'.
%                           Range: (-inf, inf).
%   - numFreqBins       -   Number of Frequency Bins.
%                           The number of Frequency Bins in the Frequency
%                           Domain.
%                           Structure: Scalar.
%                           Type: 'Single' / 'Double'.
%                           Range: {1, 2, ...}.
%   - varX              -   Variance of Signal Model.
%                           The variance of each peak in the frequency
%                           domain assuming Normal Distribution.
%                           Structure: Scalar.
%                           Type: 'Single' / 'Double'.
%                           Range: (0, inf).
%   - varN              -   Variance of Noise.
%                           The variance of the noise.
%                           Structure: Scalar.
%                           Type: 'Single' / 'Double'.
%                           Range: (0, inf).
%   - numIterations     -   Number of Iterations.
%                           Number of iterations of the algorithm.
%                           Structure: Scalar.
%                           Type: 'Single' / 'Double'.
%                           Range {1, 2, ...}.
% Output:
%   - vF                -   High Resolution DFT.
%                           High Resolution DFT generated according to the
%                           Bayesian Model.
%                           Structure: Vector (Column Vector).
%                           Type: 'Single' / 'Double'.
%                           Range: (-inf, inf).
% References
%   1.  C
% Remarks:
%   1.  T
% TODO:
%   1.  Use Levinson's Recursion to solve `vB = ((lambdaFctr * mI) + mF * mQ * mF') \ vX;`.
%   2.  Add "Stopping Condition".
% Release Notes:
%   -   1.0.000     07/11/2016
%       *   First realease version.
% ----------------------------------------------------------------------------------------------- %

paramAlphaBase = 0.00045;

mAA = mA.' * mA;
vAb = mA.' * vB;
% vX  = mAA \ vAb;
vX  = pinv(mA) * vB; %<! Dealing with "Fat Matrix"

paramAlphaBase = 4 / (1.1 * (norm(mA, 2) ^ 2));
paramAlphaBase = 1 / (2 * (norm(mA, 2) ^ 2));

mX = zeros([size(vX, 1), (numIterations + 1)]);
mX(:, 1) = vX;

vY = vX;
tK = 1;

for ii = 1:numIterations
    vYGrad      = (mAA * vY) - vAb;
    % paramAlpha  = paramAlphaBase / sqrt(ii);
    paramAlpha  = paramAlphaBase;
    
    vXPrev      = vX;
    vX          = ProxL1(vY - (paramAlpha * vYGrad), paramAlpha * lambdaFctr);
    
    tKPrev          = tK;
    tK              = (1 + sqrt(1 + (4 * tKPrev))) / 2;
    fistaStepSize   = (tKPrev - 1) / tK;
    
    vY = vX + (fistaStepSize * (vX - vXPrev));
    
    mX(:, (ii + 1)) = vX;
end


end


function [ vX ] = ProxL1( vX, lambdaFactor )

% Soft Thresholding
vX = max(vX - lambdaFactor, 0) + min(vX + lambdaFactor, 0);

end

