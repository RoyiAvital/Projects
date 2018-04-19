function [ mShortPathImg ] = ExtractShortestPath( mShortDist, mShortPath )
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

numRows = size(mShortDist, 1);
numCols = size(mShortDist, 2);

mShortPathImg               = zeros([numRows, numCols]);
[shortDist, shortDistIdx]   = min(mShortDist(numRows, :));

currRow = numRows;
currCol = shortDistIdx;

mShortPathImg(currRow, currCol) = 1;

prevEntry = mShortPath(currRow, currCol);


while(prevEntry > 0)
    [currRow, currCol] = ind2sub([numRows, numCols], prevEntry);
    
    mShortPathImg(currRow, currCol) = 1;
    
    prevEntry = mShortPath(currRow, currCol);
end


end

