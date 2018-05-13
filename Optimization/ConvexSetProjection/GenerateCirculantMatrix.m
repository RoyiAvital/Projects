function [ mX ] = GenerateCirculantMatrix( vX )
% ----------------------------------------------------------------------------------------------- %
% [ mX ] = GenerateCirculantMatrix( vX )
%   Generates Circulant Matrix from its generating first row.
% Input:
%   - vX                -   Input Vector.
%                           Structure: Vector.
%                           Type: 'Single' / 'Double'.
%                           Range: (-inf, inf).
% Output:
%   - mX                -   Circulant Matrix.
%                           Circulant Matrix where its first row is vX and
%                           rest of the rows are shifted versions of vX.
%                           Structure: Matrix.
%                           Type: 'Single' / 'Double'.
%                           Range: (-inf, inf).
% References:
%   1.  StackOverflow (Q23174345) - https://stackoverflow.com/questions/23174345.
% Remarks:
%   1.  Prefixes:
%       -   'm' - Matrix.
%       -   'v' - Vector.
%   2.  T
% TODO:
%   1.  Support Complex Matrices (Remove the 'symmetric' property).
%   Release Notes:
%   -   1.0.000     12/05/2018  Royi Avital
%       *   First release version.
% ----------------------------------------------------------------------------------------------- %

FALSE   = 0;
TRUE    = 1;

OFF = 0;
ON  = 1;

mX = ifft( fft( diag( fft(vX) ), [], 1 ), [], 2 ); %<! This is F * diag(fft(vX)) * F'

% mX = ifft( fft( diag( fft(vX) ), [], 1 ), [], 2, 'symmetric' );

% mX = real(ifft( fft( diag( fft(vX) ), [], 1 ), [], 2 ));
% mX = ifft( fft( diag( fft(real(vX)) ), [], 1 ), [], 2, 'symmetric' ); %<! More efficinet than one line above


end

