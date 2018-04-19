function [ hImageObj ] = DisplayMatrixValues( mInputData, hAxes )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

FALSE   = 0;
TRUE    = 1;
OFF     = 0;
ON      = 1;

formatSpec = '%1.2f';

numRows = size(mInputData, 1);
numCols = size(mInputData, 2);

hImageObj   = imagesc(hAxes, mInputData);

for ii = 1:numRows
    for jj = 1:numCols
        hTextObj = text(hAxes, jj, ii, num2str(mInputData(ii, jj), formatSpec), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    end
end


end

