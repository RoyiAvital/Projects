function [ mDistMat ] = CalcDistanceMatrix001( mA, mB )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

numSamples001   = size(mA, 2);
numSamples002   = size(mB, 2);

mDistMat = zeros([numSamples001, numSamples002], 'single');

for ii = 1:numSamples001
    mDistMat(ii, :) = sum((mA(:, ii) - mB) .^ 2);
end


end

