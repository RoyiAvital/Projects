function [ vX ] = SmoothSavitzkyGolayMonotonic( vX, sOpt )
% ----------------------------------------------------------------------------------------------- %
%[ estFreq ] = EstimateHarmonicFrequency( vX, samplingFreq, estType )
% Estimates the frequency of a single Real / Complex Harmonic signal with
% arbitrary real amplitude and phase. It is assumed only a single tone
% exists.
% Input:
%   - vX                -   Input Samples.
%                           The data samples of a single tone harmonic
%                           signal.
%                           Structure: Vector (numSamples X 1).
%                           Type: 'Single' / 'Double'.
%                           Range: (-inf, inf).
%   - samplingFreq      -   Sampling Frequency.
%                           Structure: Scalar.
%                           Type: 'Single' / 'Double'.
%                           Range: (0, inf).
%   - estType           -   Estimation Type.
%                           When set to 1 the algorithm uses the unweighted
%                           version (Eq 17 in the paper). When set to 2 the
%                           algorithm used the weighted version (Eq 18 in
%                           the paper).
%                           Structure: Scalar.
%                           Type: 'Single' / 'Double'.
%                           Range {1, 2}.
% Output:
%   - estFreq           -   Estimated Frequency.
%                           The estimated frequency of the single tone.
%                           Structure: Scalar.
%                           Type: 'Single' / 'Double'.
%                           Range (0, inf).
% References
%   1.  
% Remarks:
%   1.  Very useful for performance curve. For instance estimator variance
%       as a function of the SNR which usually is non decreasing.
%   2.  Should be implemented as a solution of an optimization problem:
%       \arg \min_x || x - y ||_2^2 + \lambda \sum ( x_i - x_{i + 1} )
%       subject to x_i <= x_{i + 1} or x_i >= x_{i + 1}
% Known Issues:
%   1.  C
% TODO:
%   1.  Look at some literature such as "A smoothed Monotonic Regression via L2 Regularization" (https://link.springer.com/article/10.1007/s10115-018-1201-2).
%   2.  Implement Pool Adjacent Violators (PAV).
% Release Notes:
%   -   1.0.000     11/08/2021  Royi Avital
%       *   First release version.
% ----------------------------------------------------------------------------------------------- %

arguments
    vX (:, 1) {mustBeNumeric}
    sOpt.monoDir (1, 1) {mustBeNumeric, mustBeReal, mustBePositive, mustBeInteger, mustBeMember(sOpt.monoDir, [1, 2])} = 2
    sOpt.kernelRadius (1, 1) {mustBeNumeric, mustBeReal, mustBePositive, mustBeInteger} = 2
    sOpt.polyDeg (1, 1) {mustBeNumeric, mustBeReal, mustBePositive, mustBeInteger} = 2
    sOpt.numIterations (1, 1) {mustBeNumeric, mustBeReal, mustBePositive, mustBeInteger} = 10
end

FALSE   = 0;
TRUE    = 1;
OFF     = 0;
ON      = 1;

MONO_NON_DEC = 1;
MONO_NON_INC = 2; %<! Better in high SNR

monoDir         = sOpt.monoDir;
kernelRadius    = sOpt.kernelRadius;
polyDeg         = sOpt.polyDeg;
numIterations   = sOpt.numIterations;

numSamples = size(vX, 1);

% Generate the Diff Operator (1D Gradient) by Finite Differences
switch(monoDir)
    case(MONO_NON_DEC)
        mD = spdiags([ones(numSamples, 1), -ones(numSamples, 1)], [0, 1], numSamples - 1, numSamples);
    case(MONO_NON_INC)
        mD = spdiags([-ones(numSamples, 1), ones(numSamples, 1)], [0, 1], numSamples - 1, numSamples);
end

for ii = 1:numIterations
    % Smoothing
    vX = SmoothSavitzkyGolay(vX, 'kernelRadius', kernelRadius, 'polyDeg', polyDeg);
    % Monotonic Regression (Can be solved faster with Pool Adjacent Violators (PAV)).
    sOpt    = optimoptions('lsqlin', 'Display', 'off');
    vX      = lsqlin(speye(numSamples), vX, mD, zeros(numSamples - 1, 1), [], [], [], [], vX, sOpt);
end


end

