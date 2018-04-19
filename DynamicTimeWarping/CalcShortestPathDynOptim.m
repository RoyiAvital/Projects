function [ shortDist, mShortDist, mShortPath ] = CalcShortestPathDynOptim( mDistMtx, numRows, numCols )
% ----------------------------------------------------------------------------------------------- %
% [ shortDist, mShortDist, mShortPath ] = CalcShortestPathDyn( mDistMtx, numRows, numCols )
% Cacluate the shortest path over a Distance Matrix from point (1, 1) to
% (m, n) using Dyanmic Programming.
% Input:
%   - mDistMtx      -   Disatnce Matrix.
%                       Structure: Matrix (numRows x numCols).
%                       Type: 'Single' / 'Double'.
%                       Range: [0, inf).
%   - numRows       -   Number of Rows.
%                       The number of rows of the Distance Matrix.
%                       Structure: Scalar.
%                       Type: 'Single' / 'Double'.
%                       Range: {1, 2, ...}.
%   - numCols       -   Number of Columns.
%                       The number of columns of the Distance Matrix.
%                       Structure: Scalar.
%                       Type: 'Single' / 'Double'.
%                       Range: {1, 2, ...}.
% Output:
%   - shortDist     -   The Shortest Distance Value.
%                       The minimum distance from (1, 1) to (m, n).
%                       Structure: Scalar.
%                       Type: 'Single' / 'Double'.
%                       Range: [0, 1].
%   - mShortDist    -   Dynamic Programming Disatnce Matrix.
%                       The evaluated distance to the (i, j) element
%                       following the Dynamic Programming calculation.
%                       Structure: Matrix (numRows x numCols).
%                       Type: 'Single' / 'Double'.
%                       Range: [0, inf).
%   - mShortPath    -   Dynamic Programming Direction Matrix.
%                       Specifies the path to (i, j) element from its valid
%                       neighbors according to the Dynamic Programming
%                       calculation.
%                       Structure: Matrix (numRows x numCols).
%                       Type: 'Single' / 'Double'.
%                       Range: [0, inf).
% References:
%   1.  Wikipedia Dynamic Time Warping - https://en.wikipedia.org/wiki/Dynamic_time_warping.
% Remarks:
%   1.  This implementation assumes the strating point of 2 "Walkers" was
%       the same.
% TODO:
%   1.  C
%   Release Notes:
%   -   1.0.000     09/05/2016  Royi Avital
%       *   First release version.
% ----------------------------------------------------------------------------------------------- %
 
FALSE   = 0;
TRUE    = 1;
 
OFF = 0;
ON  = 1;
 
mShortDist = zeros([numRows, numCols]);
mShortPath = zeros([numRows, numCols]);
 
mShortDist(1:numRows, 1) = cumsum(mDistMtx(1:numRows, 1));
mShortDist(1, 1:numCols) = mDistMtx(1, 1:numCols);
 
mShortPath(1:numRows, 1) = [0:(numRows - 1)].';

for jj = 2:numCols
    for ii = 2:numRows
        
        currLocDist = mDistMtx(ii, jj);
                
        % vFeasiblePathVal = [mShortDist(ii - 1, jj), ...
        %     mShortDist(ii - 1, jj - 1), ...
        %     mShortDist(ii, jj - 1)];
        % [minPathVal, minPathIdx] = min(vFeasiblePathVal);
       
        stepDist1 = mShortDist(ii - 1, jj);
        stepDist2 = mShortDist(ii - 1, jj - 1);
        stepDist3 = mShortDist(ii, jj - 1);
        minStepVal = min(min(stepDist1, stepDist2), stepDist3);
        switch(minStepVal)
            case(stepDist1)
                mShortPath(ii, jj) = ((jj - 1) * numRows) + (ii - 1);
            case(stepDist2)
                mShortPath(ii, jj) = (((jj - 1) - 1) * numRows) + (ii - 1);
            case(stepDist3)
                mShortPath(ii, jj) = (((jj - 1) - 1) * numRows) + ii;
        end
        
        mShortDist(ii, jj) = currLocDist + minStepVal;
    end
end
 
shortDist = min(mShortDist(numRows, :));
 
 
end

