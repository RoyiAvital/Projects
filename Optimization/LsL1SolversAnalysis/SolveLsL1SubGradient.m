function [ vX, mX ] = SolveLsL1SubGradient( mA, vB, paramLambda, numIterations )
% ----------------------------------------------------------------------------------------------- %
%[ vX, mX ] = SolveLsL1SubGradient( mA, vB, paramLambda, numIterations )
% Solve L1 Regularized Least Squares Using Sub Gradient (PGM) Method.
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
%   1.  Wikipedia Sub Gradient Method - https://en.wikipedia.org/wiki/Subgradient_method.
% Remarks:
%   1.  Using 4 options for Step Size Policy. 
% Known Issues:
%   1.  A
% TODO:
%   1.  B
% Release Notes:
%   -   1.0.000     23/08/2017
%       *   First release version.
% ----------------------------------------------------------------------------------------------- %

STEP_SIZE_POLICY_CONST              = 1;
STEP_SIZE_POLICY_SCALED             = 2;
STEP_SIZE_POLICY_NON_SUMMABLE       = 3;
STEP_SIZE_POLICY_SQUARE_SUMMABLE    = 4;


mAA = mA.' * mA;
vAb = mA.' * vB;
vX  = pinv(mA) * vB; %<! Dealing with "Fat Matrix"

stepSizeLipschitz   = 1 / (2 * (norm(mA, 2) ^ 2));
stepSizePolicy      = STEP_SIZE_POLICY_CONST;

mX = zeros([size(vX, 1), numIterations]);
mX(:, 1) = vX;

for ii = 2:numIterations
    
    vG = (mAA * vX) - vAb + (paramLambda * sign(vX));
    
    switch(stepSizePolicy)
        case(STEP_SIZE_POLICY_CONST)
            stepSize = stepSizeLipschitz;
        case(STEP_SIZE_POLICY_SCALED)
            stepSize = stepSizeLipschitz / norm(vG);
        case(STEP_SIZE_POLICY_NON_SUMMABLE)
            stepSize = 1 / sqrt(ii);
        case(STEP_SIZE_POLICY_SQUARE_SUMMABLE)
            stepSize = 1 / ii;
    end
    
    vX = (vX - stepSize * vG);
    
    
    mX(:, ii) = vX;
    
end


end

