function [ mDistMat ] = CalcDistanceMatrix002( mDataSamples001, mDataSamples002 )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

mDistMat = sum(mDataSamples001 .^ 2).' - (2 * mDataSamples001.' * mDataSamples002) + sum(mDataSamples002 .^ 2);


end

