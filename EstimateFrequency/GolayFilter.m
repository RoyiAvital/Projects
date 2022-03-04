
clear();
close('all');
load('DataSmoothing.mat');

% Parameters
polyDeg = 3;
kernelRadius = 2;

% Solution

mJ = (-kernelRadius:kernelRadius).' .^ (0:polyDeg);
mC = (mJ.' * mJ) \ mJ.'; %<! Wikipedia

% The first row are the filter coefficients
mC * 35

[mQ, ~] = qr(mJ, 0);
vK = mQ * mQ(kernelRadius + 1, :).'


mFreqErr = 10 * log10(mFreqErr);

figure();
plot(mFreqErr);

for ii = 1:size(mFreqErr, 2)
    % mFreqErr(:, ii) = filtfilt(vK, 1, mFreqErr(:, ii));
    mFreqErr(:, ii) = SmoothSavitzkyGolayMonotonic(mFreqErr(:, ii), 'kernelRadius', 15, 'numIterations', 35);
end

figure();
plot(mFreqErr);