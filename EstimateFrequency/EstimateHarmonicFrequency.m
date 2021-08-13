function [ estFreq ] = EstimateHarmonicFrequency( vX, samplingFreq, estType )
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
    samplingFreq (1, 1) {mustBeNumeric, mustBeReal, mustBePositive}
    estType (1, 1) {mustBeNumeric, mustBeReal, mustBePositive, mustBeInteger, mustBeMember(estType, [1, 2])} = 2
end

FALSE   = 0;
TRUE    = 1;
OFF     = 0;
ON      = 1;

EST_TYPE_1 = 1;
EST_TYPE_2 = 2; %<! Better in high SNR

numSamples = size(vX, 1);

if(isreal(vX))
    % Sine / Cosine Signal -> Use Hilbert Transform to generate the
    % Analytic Signal
    vX = hilbert(vX);
end

estFreq = 0;
switch(estType)
    case(EST_TYPE_1)
        for ii = 1:(numSamples - 1)
            estFreq = estFreq + angle(vX(ii)' * vX(ii + 1));
        end
        estFreq = estFreq / (2 * pi * (numSamples - 1));
    case(EST_TYPE_2)
        weightDen = (numSamples ^ 3) - numSamples;
        for ii = 1:(numSamples - 1)
            sampleWeight    = (6 * ii * (numSamples - ii)) / weightDen;
            estFreq         = estFreq + (sampleWeight * angle(vX(ii)' * vX(ii + 1)));
        end
        estFreq = estFreq / (2 * pi);
end

estFreq = samplingFreq * estFreq; %<! Moving from normalized frequency to absolute frequency


end

