function [ vX, mX ] = SolveLsL0Prox( mA, vB, paramLambda, numIterations )
% ----------------------------------------------------------------------------------------------- %
%[ vX, mX ] = SolveLsL0Prox( mA, vB, paramLambda, numIterations )
% Solve L0 Regularized Least Squares Using Proximal Gradient (PGM) Method.
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
%   2.  Bummer, thought it is new, Yet Michael Elad referred me to
%       "Iterative Thresholding for Sparse Approximations" by Thomas 
%       Blumensath and Mike E. Davies.
% Remarks:
%   1.  Using vanilla PGM.
% Known Issues:
%   1.  A
% TODO:
%   1.  B
% Release Notes:
%   -   1.0.000     30/03/2018
%       *   First release version.
% ----------------------------------------------------------------------------------------------- %

mAA = mA.' * mA;
vAb = mA.' * vB;
vX  = pinv(mA) * vB; %<! Dealing with "Fat Matrix"

stepSize = 1 / (2 * (norm(mA, 2) ^ 2));
% stepSize = 1 / sum(mA(:) .^ 2); %<! Faster to calculate, conservative (Hence slower)

mX = zeros([size(vX, 1), numIterations]);
mX(:, 1) = vX;

for ii = 2:numIterations
    
    vG = (mAA * vX) - vAb;
    vX = ProxL0(vX - (stepSize * vG), stepSize * paramLambda);
    
    mX(:, ii) = vX;
    
end


end


function [ vX ] = ProxL0( vX, lambdaFactor )
% See https://math.stackexchange.com/questions/2713427

% Soft Thresholding
vX( (0.5 * (vX .* vX)) <= lambdaFactor ) = 0;


end

