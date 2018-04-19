function [  ] = DrawLinesBetweenCurves( vY1Idx, vY1, vY2Idx, vY2, hAxes )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

FALSE   = 0;
TRUE    = 1;
OFF     = 0;
ON      = 1;

lineWidth = 1;
vLineColor = 0.35 * ones([1, 3]);

numLines = size(vY1Idx, 1);

for ii = 1:numLines
    vX = [vY1Idx(ii), vY2Idx(ii)];
    vY = [vY1(vY1Idx(ii)), vY2(vY2Idx(ii))];
    hLineSeries = line(hAxes, vX, vY);
    set(hLineSeries, 'LineWidth', lineWidth, 'Color', vLineColor);
end


end

