function [ vX, mX ] = SolveLsL1Admm( mA, vB, paramLambda, numIterations )
% ----------------------------------------------------------------------------------------------- %
%[ vX, mX ] = SolveLsL1Admm( mA, vB, lambdaFctr, numIterations )
% Solve L1 Regularized Least Squares Using Alternating Direction Method of Multipliers (ADMM) Method.
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
%   1.  Wikipedia ADMM - https://en.wikipedia.org/wiki/Augmented_Lagrangian_method#Alternating_direction_method_of_multipliers.
% Remarks:
%   1.  Using vanilla ADMM with no optimization of the parameter or
%       smoothing.
% Known Issues:
%   1.  A
% TODO:
%   1.  Pre calculate decomposition of the Linear System.
% Release Notes:
%   -   1.0.000     23/08/2017
%       *   First realease version.
% ----------------------------------------------------------------------------------------------- %

mAA = mA.' * mA;
vAb = mA.' * vB;
vX  = pinv(mA) * vB; %<! Dealing with "Fat Matrix"

paramRho    = 5;
mI          = eye(size(vX, 1));

vZ = vX;
vU = vX;

mX = zeros([size(vX, 1), numIterations]);
mX(:, 1) = vX;

for ii = 2:numIterations
    
    vX = (mAA + (paramRho * mI)) \ (vAb + (paramRho * vZ) - vU);
    vZ = ProxL1(vX + (vU / paramRho), paramLambda / paramRho);
    vU = vU + (paramRho * (vX - vZ));
    
    mX(:, ii) = vX;
    
end


end


function [ vX ] = ProxL1( vX, lambdaFactor )

% Soft Thresholding
vX = max(vX - lambdaFactor, 0) + min(vX + lambdaFactor, 0);


end

