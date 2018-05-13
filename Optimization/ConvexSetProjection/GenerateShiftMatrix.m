function [ mPi ] = GenerateShiftMatrix( matDim, shiftDeg )
% ----------------------------------------------------------------------------------------------- %
% [ mPi ] = GenerateShiftMatrix( matDim, shiftDeg )
%   Generates Shifting Matrices which shifts a vector by shiftDeg.
% Input:
%   - matDim            -   Matrix Dimension.
%                           Structure: Scalar.
%                           Type: 'Single' / 'Double'.
%                           Range: {1, 2, ...}.
%   - shiftDeg          -   Shifting Degree.
%                           Number of shits of the matrix where shiftDeg =
%                           0 means the idnentity matrix.
%                           Structure: Scalar.
%                           Type: 'Single' / 'Double'.
%                           Range: {0, 1, 2, ..., matDim - 1}.
% Output:
%   - mX                -   Circulant Matrix.
%                           Structure: Matrix.
%                           Type: 'Single' / 'Double'.
%                           Range: (-inf, inf).
% References:
%   1.  StackOverflow (Q2778195) - https://math.stackexchange.com/a/2778260/33.
% Remarks:
%   1.  Prefixes:
%       -   'm' - Matrix.
%       -   'v' - Vector.
%   2.  Pay attention that GenerateShiftMatrix( matDim, 5 ) =
%       GenerateShiftMatrix( matDim, 1 ) ^ 5. Namely shifting mPi can be
%       done by multiplication of it by itself.
% TODO:
%   1.  S
%   Release Notes:
%   -   1.0.000     12/05/2018  Royi Avital
%       *   First release version.
% ----------------------------------------------------------------------------------------------- %

FALSE   = 0;
TRUE    = 1;

OFF = 0;
ON  = 1;

% The shifting matrix is permutation of the rows of I.
mPi = eye(matDim);
mPi = [mPi(shiftDeg + 1:end, :); mPi(1:shiftDeg, :)];


end

