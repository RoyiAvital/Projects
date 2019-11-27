function [ vX, mX ] = SolveLsTvChambolle( vX, mA, vB, mD, paramLambda, numIterations )
% ----------------------------------------------------------------------------------------------- %
%[ vX, mX ] = SolveLsTvChambolle( vX, mA, vB, mD, paramLambda, numIterations )
% Solve Total Variation (TV) Regularized Least Squares Using Chambolle's 
% Method. Basically solves the problem given by:
% $$ \arg \min_{ x \in \mathbb{R}^{n} } \frac{1}{2} {\left\| A x - b \right|}_{2}^{2} + \lambda {\left\| D x \right\|}_{1} $$
% Input:
%   - vX                -   input Vector.
%                           Initialization of the iterative process.
%                           Structure: Vector (n X 1).
%                           Type: 'Single' / 'Double'.
%                           Range: (-inf, inf).
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
%   -   1.0.000     27/11/2019  Royi Avital
%       *   First realease version.
% ----------------------------------------------------------------------------------------------- %

STEP_SIZE_POLICY_CONST              = 1;
STEP_SIZE_POLICY_SCALED             = 2;
STEP_SIZE_POLICY_NON_SUMMABLE       = 3;
STEP_SIZE_POLICY_SQUARE_SUMMABLE    = 4;
STEP_SIZE_POLICY_NON_SUMMABLE_LIP   = 5;

mDD = mD * mD.';
% vDB = mD * vB;
vAB = mA.' * vB;
mAAInv = inv(mA.' * mA); %<! Useful Constant
vC = paramLambda * mD * mAAInv * vAB; %<! Useful Constant
mC = (paramLambda * paramLambda) * mD * mAAInv * mD.'; %<! Useful Constant

if(isempty(vX))
    vX  = pinv(mA) * vB; %<! Dealing with "Fat Matrix"
end

% Assuming the Norm (Maximum Eigen Value) of mD is 1
stepSizeLipschitz   = 1 / (2 * (paramLambda * norm(mA, 2) ^ 2));
stepSizePolicy      = STEP_SIZE_POLICY_NON_SUMMABLE_LIP;

vP = (mDD \ (mD * mA.' * ((mA * vX) - vB))) / paramLambda;

mX = zeros([size(vX, 1), numIterations]);
mX(:, 1) = vX;

for ii = 2:numIterations
    vG = (mC * vP) - vC; %<! Gradient of vP
    
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
    
    % This is the by Lipschitz Constant
    stepSize = 1 / (8 * paramLambda * paramLambda);
    
    vP = vP - (stepSize * vG);
    vP = vP ./ (max(1, abs(vP)));
    
    vX = mAAInv * (vAB - paramLambda * mD.' * vP);
    
    mX(:, ii) = vX;
end


end

