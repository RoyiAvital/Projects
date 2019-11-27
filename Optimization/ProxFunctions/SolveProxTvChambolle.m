function [ vX ] = SolveProxTvChambolle( vY, mD, paramLambda, numIterations )
% ----------------------------------------------------------------------------------------------- %
%[ vX ] = SolveLsTvChambolle( vX, mA, vB, mD, paramLambda, numIterations )
% Solves the Prox of the Total Variation (TV) Norm using Chambolle's
% Method. Basically solves the problem given by:
% $$ \arg \min_{ x \in \mathbb{R}^{n} } \frac{1}{2} {\left\| x - y \right|}_{2}^{2} + \lambda {\left\| D x \right\|}_{1} $$
% Input:
%   - vY                -   input Vector.
%                           The model known data.
%                           Structure: Vector (m X 1).
%                           Type: 'Single' / 'Double'.
%                           Range: (-inf, inf).
%   - mD                -   Finite Differences Operator Matrix.
%                           Matrix which implemets the Forward Finite
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
%   1.  A
% Remarks:
%   1.  B
% Known Issues:
%   1.  A
% TODO:
%   1.  B
% Release Notes:
%   -   1.0.000     27/11/2019  Royi Avital
%       *   First realease version.
% ----------------------------------------------------------------------------------------------- %

STEP_SIZE_POLICY_CONST              = 1;
STEP_SIZE_POLICY_SCALED             = 2;
STEP_SIZE_POLICY_NON_SUMMABLE       = 3;
STEP_SIZE_POLICY_SQUARE_SUMMABLE    = 4;
STEP_SIZE_POLICY_NON_SUMMABLE_LIP   = 5;

mDD = mD * mD.';
vDy = mD * vY;
mDDLambda = paramLambda * paramLambda * mDD;
vDbLambda = paramLambda * vDy;


% Assuming the Norm (Maximum Eigen Value) of mD is 1
stepSizeLipschitz   = 1 / (8 * paramLambda * paramLambda);
stepSizePolicy      = STEP_SIZE_POLICY_NON_SUMMABLE_LIP;

vX = vY;
vP = (mDDLambda \ (mD * (vX - vY))) / paramLambda;


for ii = 2:numIterations
    vG = mDDLambda * vP - vDbLambda; %<! Gradient of vP
    
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
    
    vP = vP - (stepSize * vG);
    vP = vP ./ (max(1, abs(vP)));
    
    vX = vY - (paramLambda * mD.' * vP);
    
end


end

