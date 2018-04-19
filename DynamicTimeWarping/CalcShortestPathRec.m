function [ shortDist ] = CalcShortestPathRec( mDistMtx, numRows, numCols )
% ----------------------------------------------------------------------------------------------- %
% [ shortDist, mShortDist, mShortPath ] = CalcShortestPathRec( mDistMtx, numRows, numCols )
% Cacluate the shortest path over a Distance Matrix from point (1, 1) to
% (m, n) using Recursive approach.
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

shortDist = sum(mDistMtx(:));

for ii = 0:(numCols - 1)
    currShortDist = CalcShortestPath(mDistMtx, numRows, numCols - ii);
    if(currShortDist < shortDist)
        shortDist = currShortDist;
    end
end


end

function [ shortDist ] = CalcShortestPath( mDistMtx, numRows, numCols )

if(numRows == 1)
    shortDist = min(mDistMtx(1, 1:numCols));
    return;
end

if(numCols == 1)
    shortDist = sum(mDistMtx(1:numRows, numCols));
    return;
end


vFeasiblePathVal = [CalcShortestPath(mDistMtx, numRows - 1, numCols), ...
    CalcShortestPath(mDistMtx, numRows - 1, numCols - 1), ...
    CalcShortestPath(mDistMtx, numRows, numCols - 1)];

shortDist = mDistMtx(numRows, numCols) + min(vFeasiblePathVal);


end

