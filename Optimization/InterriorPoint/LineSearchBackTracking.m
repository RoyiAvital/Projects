function [ stepSize ] = LineSearchBackTracking( hObjFun, vX, vD, paramAlpha )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

stepSize = 1;

vY      = vX + (stepSize * vD);
objVal  = hObjFun(vX);

while(hObjFun(vY) > objVal)
    
    stepSize    = paramAlpha * stepSize;
    vY          = vX + (stepSize * vD);
    
end


end

