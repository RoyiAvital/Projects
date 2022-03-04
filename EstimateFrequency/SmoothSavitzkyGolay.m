function [ vX ] = SmoothSavitzkyGolay( vX, sOpt )
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
%   1.  Steven Kay - A Fast and Accurate Single Frequency Estimator (https://ieeexplore.ieee.org/document/45547).
%   2.  Signal Processing StackExchange Q76644 - Simple and Effective Method to Estimate the Frequency of a Single Sine Signal in White Noise (https://dsp.stackexchange.com/questions/76644).
% Remarks:
%   1.  Works only for a single tone signals of the model (One of):
%       x[n] = A * sin[2 * pi * f / fs * n + phi] + w[n]
%       x[n] = A * exp[-1i * (2 * pi * f / fs * n + phi)] + w[n]
%   2.  Real signals (Sine / Cosine) will require Hilbert Transform in
%       order to be transformed into their analytic form. As the core
%       algorithm assumes the model:
%       x[n] = A * exp[-1i * (2 * pi * f / fs * n + phi)] + w[n]
% Known Issues:
%   1.  C
% TODO:
%   1.  D
% Release Notes:
%   -   1.0.000     11/08/2021  Royi Avital
%       *   First release version.
% ----------------------------------------------------------------------------------------------- %

arguments
    vX (:, 1) {mustBeNumeric}
    sOpt.kernelRadius (1, 1) {mustBeNumeric, mustBeReal, mustBePositive, mustBeInteger} = 2
    sOpt.polyDeg (1, 1) {mustBeNumeric, mustBeReal, mustBePositive, mustBeInteger} = 2
end

FALSE   = 0;
TRUE    = 1;
OFF     = 0;
ON      = 1;

kernelRadius    = sOpt.kernelRadius;
polyDeg         = sOpt.polyDeg;

numSamples = size(vX, 1);

% Validte the number of samples used by filter isn't longer than the
% signal.
% (numSamples - (mod(numSamples, 2) == 0) = An odd number
kernelRadius = min(kernelRadius, ((numSamples - (mod(numSamples, 2) == 0)) - 1) / 2);
kernelLength = (2 * kernelRadius) + 1;

mJ = (-kernelRadius:kernelRadius).' .^ (0:polyDeg);
% mC = (mJ.' * mJ) \ mJ.'; %<! Wikipedia: The first row are the filter coefficients

[mQ, ~] = qr(mJ, 0);
vK = mQ * mQ(kernelRadius + 1, :).';

vX1 = mQ(1:kernelRadius, :) * mQ.' * vX(1:kernelLength); %<! Begining
vX2 = conv2(vX, vK, 'valid'); %<! Middle valid
vX3 = mQ((kernelRadius + 2):end, :) * mQ.' * vX(numSamples - kernelLength + 1:numSamples); %<! End

vX(:) = cat(1, vX1, vX2, vX3);


end

