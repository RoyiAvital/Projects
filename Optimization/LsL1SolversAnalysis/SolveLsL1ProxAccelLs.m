function [ vX, mX ] = SolveLsL1ProxAccelLs( mA, vB, lambdaFctr, numIterations )
% ----------------------------------------------------------------------------------------------- %
%[ vX, mX ] = SolveLsL1ProxAccelLs( mA, vB, paramLambda, numIterations )
% Solve L1 Regularized Least Squares Using Accelerated Proximal Gradient (PGM) Method with Line Search.
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
%   1.  Wikipedia PGM - https://en.wikipedia.org/wiki/Proximal_gradient_method.
%   2.  Wikipedia Fast Gradient Methods - https://en.wikipedia.org/wiki/Gradient_descent#Fast_gradient_methods.
% Remarks:
%   1.  Using Smooth / Aggressive FISTA step size (The smooth is the
%       original step size in the FISTA article).
%   2.  The line search is looking for the Lipschitz Constant of the
%       function from the condition it must obey.
% Known Issues:
%   1.  A
% TODO:
%   1.  B
% Release Notes:
%   -   1.0.000     23/08/2017
%       *   First release version.
% ----------------------------------------------------------------------------------------------- %

FISTA_STEP__METHOD_SMOOTH       = 1; %<! Monotonic, Slower
FISTA_STEP__METHOD_AGGRESSIVE   = 2; %<! Non Monotonic, Faster

hObjFun = @(vX) (0.5 * sum(((mA * vX) - vB) .^ 2)) + (lambdaFctr * sum(abs(vX)));

mAA = mA.' * mA;
vAb = mA.' * vB;
vX  = pinv(mA) * vB; %<! Dealing with "Fat Matrix"

lipConst    = 1; %<! Lipschitz Constant
% deltaThr - Slackness for the deltaVal (Given the true Lipschitz Constant
% should be positive) requiring bigger than small positive number.
deltaThr    = 1e-3;
paramBeta   = 1.05; %<! Back Tracking parameter
fistaStepMode = FISTA_STEP__METHOD_AGGRESSIVE;

mX = zeros([size(vX, 1), (numIterations + 1)]);
mX(:, 1) = vX;

vY = vX;
tK = 1;

for ii = 1:numIterations
    
    objVal  = hObjFun(vX);
    
    vYGrad      = (mAA * vY) - vAb;
    vG          = (mAA * vX) - vAb;
    
    vXPrev      = vX;
    deltaVal    = 1;
    
    while(deltaVal > deltaThr)
        
        stepSize    = 1 / lipConst;
        vX          = ProxL1(vY - (stepSize * vYGrad), stepSize * lambdaFctr);
        % For the true value of Lipschitz Constant deltaVal should be
        % negative (Assuming no numerical issues).
        deltaVal    =  hObjFun(vX) - objVal - (vG.' * (vX - vXPrev)) - ((lipConst / 2) * sum((vX - vXPrev) .^ 2));
        
        lipConst = lipConst * paramBeta;
        
    end
    
    lipConst = lipConst / paramBeta;
    
    switch(fistaStepMode)
        case(FISTA_STEP__METHOD_SMOOTH)
            tKPrev          = tK;
            tK              = (1 + sqrt(1 + (4 * tKPrev))) / 2;
            fistaStepSize   = (tKPrev - 1) / tK;
        case(FISTA_STEP__METHOD_AGGRESSIVE)
            fistaStepSize = (ii - 1) / (ii + 2);
    end
    
    vY = vX + (fistaStepSize * (vX - vXPrev));
    
    mX(:, (ii + 1)) = vX;
    
end


end


function [ vX ] = ProxL1( vX, lambdaFactor )

% Soft Thresholding
vX = max(vX - lambdaFactor, 0) + min(vX + lambdaFactor, 0);

end

