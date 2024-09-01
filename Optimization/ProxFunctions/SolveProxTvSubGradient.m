function [ vX ] = SolveProxTvSubGradient( vY, mD, paramLambda, numIterations )
% ----------------------------------------------------------------------------------------------- %
%[ vX ] = SolveProxTvSubGradient( vY, mD, paramLambda, numIterations )
% Solves the Prox of the Total Variation (TV) Norm using Sub Gradient
% Method (SGM). Basically solves the problem given by:
% $$ \arg \min_{ x \in \mathbb{R}^{n} } \frac{1}{2} {\left\| x - y \right|}_{2}^{2} + \lambda {\left\| D x \right\|}_{1} $$
% Input:
%   - vY                -   input Vector.
%                           The model known data.
%                           Structure: Vector (m X 1).
%                           Type: 'Single' / 'Double'.
%                           Range: (-inf, inf).
%   - mD                -   Finite Differences Operator Matrix.
%                           Matrix which implements the Forward Finite
%                           Differences operator.
%                           Structure: Matrix ((n - 1) X n).
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
%   2.  While in TV form the matrix 'mD' is indeed the Finite Differences
%       operator this method will solve the problem for any matrix 'mD'.
% Known Issues:
%   1.  A
% TODO:
%   1.  B
% Release Notes:
%   -   1.0.000     26/11/2019  Royi Avital
%       *   First release version.
% ----------------------------------------------------------------------------------------------- %

STEP_SIZE_POLICY_CONST              = 1;
STEP_SIZE_POLICY_SCALED             = 2;
STEP_SIZE_POLICY_NON_SUMMABLE       = 3;
STEP_SIZE_POLICY_SQUARE_SUMMABLE    = 4;
STEP_SIZE_POLICY_NON_SUMMABLE_LIP   = 5;


vX = vY;

stepSizeLipschitz   = 1 / 2;
stepSizePolicy      = STEP_SIZE_POLICY_NON_SUMMABLE_LIP;

for ii = 2:numIterations
    
    vG = vX - vb + (paramLambda * mD.' * sign(mD * vX));
    
    switch(stepSizePolicy)
        case(STEP_SIZE_POLICY_CONST)
            stepSize = stepSizeLipschitz;
        case(STEP_SIZE_POLICY_SCALED)
            stepSize = stepSizeLipschitz / norm(vG);
        case(STEP_SIZE_POLICY_NON_SUMMABLE)
            stepSize = 1 / sqrt(ii);
        case(STEP_SIZE_POLICY_SQUARE_SUMMABLE)
            stepSize = 1 / ii;
        case(STEP_SIZE_POLICY_NON_SUMMABLE_LIP)
            stepSize = (1 / sqrt(ii)) * stepSizeLipschitz;
    end
    
    vX = (vX - stepSize * vG);
    
end


end

