function [ mDistMat ] = CalcDistanceMatrix002( mA, mB )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

mDistMat = sum(mA .^ 2).' - (2 * mA.' * mB) + sum(mB .^ 2);


end

