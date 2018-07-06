% DFT Matrix Multiplication Analysis
% How to rewrite DFT Matrix Multiplication by 'fft()'.
% References:
%   1.  See https://github.com/RoyiAvital/Projects/tree/master/Optimization/ConvexSetProjection.
%   2.  See https://github.com/RoyiAvital/StackExchangeCodes/blob/master/SignalProcessing/Q19646/Q19646.ipynb.
% Remarks:
%   1.  sa
% TODO:
% 	1.  ds
% Release Notes
% - 1.0.000     05/07/2018  Royi Avital
%   *   First release.


%% General Parameters

% subStreamNumberDefault = 2116;

run('InitScript.m');
addpath(genpath('AuxiliaryFunctions'));

figureIdx           = 0;
figureCounterSpec   = '%04d';
generateFigures     = OFF;


%% Setting Parameters

numRows = 6;
numCols = 4;


%% Generate Data

mW = ([0:(numRows - 1)].' * [0:(numCols - 1)]);
mF = exp(((-2i * pi) / numRows) * mW) / sqrt(numRows); %<! DFT Matrix (Normalized)

vD = randn([numCols, 1]);
mD = diag(vD);

vE = randn([numRows, 1]);
mE = diag(vE);


%% Analysis
% Try implement mF * mD * mF' efficiently where mD is diagonal

mARef   = mF * mD * mF';
mA      = ifft(fft(mD, numRows, 1), numRows, 2);

norm(mA - mARef, 'fro')

% Try implement mF' * mE * mF efficiently where mE is diagonal
% We'll use 'mF' for the Inverse DFT Matrix -> Implementing mF * mE * mF'.
% If we go down in the number of elements we can use fft(x, numElements,
% dim) but we need to do FFT on all samples and take the interesting ones.

mF = exp(((2i * pi) / numRows) * mW.') / numRows; %<! Inverse DFT Matrix

mARef   = mF * mE * mF';
mA      = ifft(mE, [], 1);
mA      = fft(mA(1:numCols, :), [], 2) / numRows;
mA      = mA(:, 1:numCols);
% mA      = fft(ifft(mE, numCols, 2), numCols, 1);

norm(mA - mARef, 'fro')

% numRows / numCols

mFF = exp(((-2i * pi) / numRows) * mW); %<! Forward
mFB = exp(((2i * pi) / numRows) * mW.') / numRows; %<! Backward

vF = randn([numRows, 1]);
vXRef = ifft(vF);
vXRef = vXRef(1:numCols);
vX = mFB * vF;

norm(vX - vXRef)

% Since the mtrices are diagonal the result could be done on the first
% column / row and then the other could be just the result of shifted
% signal (Multiplication by Complex Exponential). Yet in practice it might
% be faster just neglecting the fact this is diagonal.



%% Restore Defaults

% set(0, 'DefaultFigureWindowStyle', 'normal');
% set(0, 'DefaultAxesLooseInset', defaultLoosInset);

