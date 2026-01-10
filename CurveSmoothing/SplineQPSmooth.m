function [ vX, isConv ] = SplineQPSmooth( vY, mD, paramLambda, vI, mA, sParams )
% ----------------------------------------------------------------------------------------------- %
% [ vX, isConv ] = SplineQPSmooth( vY, mD, paramLambda, vI, mA, sParams )
%   Applies 1D smoothing using Spline like model with reference points and
%   monotonicity.
% Input:
%   - vY            -   Input Vector.
%                       The set of samples to smooth
%                       Structure: Vector (Column).
%                       Type: 'Single' / 'Double'.
%                       Range: (-inf, inf).
%   - mD            -   Finite Differences Operator.
%                       Applies finite differences on a vector.
%                       Structure: Matrix (Dense / Sparse).
%                       Type: 'Single' / 'Double'.
%                       Range: (-inf, inf).
%   - paramLambda   -   Regularization (Smoothing) Factor.
%                       Sets the level of smoothing.
%                       Structure: Scalar.
%                       Type: 'Single' / 'Double'.
%                       Range: (0, inf).
%   - vI            -   Segments Indices Vector.
%                       Each vI(k) and vI(k + 1) define a segment which
%                       should be monotonic.
%                       Also defines points of equality: vX(vI(k)) = vY(vI(k)).
%                       Structure: Vector (Column).
%                       Type: 'Single' / 'Double'.
%                       Range: {1, 2, ..., }.
%   - mA            -   Monotonic Operator.
%                       Each row forces decreasing or increasing property
%                       on a sample.
%                       Structure: Matrix (Sparse).
%                       Type: 'Single' / 'Double'.
%                       Range: {-1, 0, 1}.
% Output:
%   - vG            -   Gradient Vector.
%                       The numerical approximation of the gradient of the
%                       Objective Function at the input point 'vX'.
%                       Structure: Vector (Column).
%                       Type: 'Single' / 'Double'.
%                       Range: (-inf, inf).
%   - isConv        -   Converging Flag.
%                       Sets the type of convergence of the solver.
%                       Structure: Scalar.
%                       Type: 'Logical'.
%                       Range: {false, true}.
% References
%   1.  https://scicomp.stackexchange.com/questions/45334.
% Remarks:
%   1.  B.
% TODO:
%   1.  Optimize the multiplication by `mE` as it is equivalent to setting
%       values.
% Release Notes:
%   -   1.0.000     09/01/2026  Royi Avital
%       *   First release version.
% ----------------------------------------------------------------------------------------------- %

arguments(Input)
    vY (:, 1) {mustBeNumeric, mustBeFinite, mustBeReal}
    mD (:, :) {mustBeNumeric, mustBeFinite, mustBeReal}
    paramLambda (1, 1) {mustBeNumeric, mustBeFinite, mustBeReal, mustBeNonnegative}
    vI (:, 1) {mustBeNumeric, mustBeFinite, mustBeReal, mustBeInteger}
    mA (:, :) {mustBeNumeric, mustBeFinite, mustBeReal}
    sParams.paramRho (1, 1) {mustBeNumeric, mustBeFinite, mustBeReal, mustBePositive} = 1.0
    sParams.numIter (1, 1) {mustBeNumeric, mustBeFinite, mustBeInteger, mustBePositive} = 5000
    sParams.epsAbs (1, 1) {mustBeNumeric, mustBeFinite, mustBeReal, mustBePositive} = 1e-5
    sParams.epsRel (1, 1) {mustBeNumeric, mustBeFinite, mustBeReal, mustBePositive} = 1e-5
    sParams.convInterval (1, 1) {mustBeNumeric, mustBeFinite, mustBeInteger, mustBePositive} = 25
    sParams.paramTau (1, 1) {mustBeNumeric, mustBeFinite, mustBeReal, mustBePositive} = 10
end

arguments(Output)
    vX (:, 1) {mustBeNumeric, mustBeFinite, mustBeReal}
    isConv (1, 1) {mustBeA(isConv, 'logical')}
end

numSamples = length(vY);
numEq      = length(vI);
numInEq    = size(mA, 1);

% Quadratic Terms
mQ = sparse(eye(numSamples)) + paramLambda * (mD' * mD);
vQ = -vY;

% Equality Constraints
mE = sparse(1:numEq, vI, 1, numEq, numSamples);
vD = vY(vI);

paramRho     = sParams.paramRho;
paramRhoInv  = inv(paramRho);
numIter      = sParams.numIter;
epsAbs       = sParams.epsAbs;
epsRel       = sParams.epsRel;
convInterval = sParams.convInterval;
paramTau     = sParams.paramTau;

% ADMM Variables
vX  = vY;                %<! Optimization variable
vS  = zeros(numInEq, 1); %<! Slack variable for inequality
vS1 = zeros(numInEq, 1); %<! Previous iteration buffer
vMu = zeros(numInEq, 1); %<! Dual variable for inequality
vNu = zeros(numEq, 1);   %<! Dual variable for equality

% Factorize the KKT System
% (Q + rho * A' * A + rho * E' * E) * z = r
mK = mQ + paramRho * (mA.' * mA) + paramRho * (mE.' * mE);
sK = decomposition(mK, 'chol', 'CheckCondition', false);

isConv     = false;
updatedRho = false;

for ii = 1:numIter
    vS1(:) = vS; %<! Previous iteration

    % Solve the Linear System
    vR = -vQ - mA.' * (paramRho * vS + vMu) + mE.' * (paramRho * vD - vNu); %<! Right hand vector
    vX = sK \ vR;

    % Proximal / Projection Step
    % s = -A * z with s >= 0
    vS = max(0, -(mA * vX + paramRhoInv * vMu));

    % Update Dual Variables
    vMu = vMu + paramRho * (mA * vX + vS);
    vNu = vNu + paramRho * (mE * vX - vD);

    % Check Convergence
    if mod(ii, convInterval) == 0
        primRes = norm(mA * vX + vS, 'inf');
        dualRes = norm(paramRho * mA' * (vS - vS1), 'inf');
        if ((primRes < epsAbs) && (dualRes < epsAbs))
            isConv = true;
            break;
        end

        % Adapt `paramRho`
        resRatio = primRes / dualRes;
        % fprintf('Primal Residual: %0.7f, Dual Residual: %0.7f\n', primRes, dualRes);
        % fprintf('Residual Ratio: %0.2f, ρ = %0.3f\n', resRatio, paramRho);
        if (resRatio > paramTau) || (inv(resRatio) > paramTau)
            updatedRho = true;
        else
            updatedRho = false;
        end
        if updatedRho
            paramRho = paramRho * sqrt(resRatio);
            paramRho = clip(paramRho, 1e-5, 1e5);
            paramRhoInv  = inv(paramRho);
            mK = mQ + paramRho * (mA.' * mA) + paramRho * (mE.' * mE);
            sK = decomposition(mK, 'chol', 'CheckCondition', false);
        end
    end

end

end

