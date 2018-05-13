function [ mX ] = ProjectCirculantMatrixSet( mY, domainFlag )
% ----------------------------------------------------------------------------------------------- %
% [ mX ] = ProjectCirculantMatrixSet( mY )
%   Projects arbitrary matrix into the convex set of Circulant Matrices to
%   yield the closest Circulant Matrix to mY in the Frobenius Norm sense.
% Input:
%   - mY                -   Input Matrix.
%                           Structure: Matrix.
%                           Type: 'Single' / 'Double'.
%                           Range: (-inf, inf).
%   - domainFlag        -   Domain Flag.
%                           Sets whether the domain is real or complex.
%                           Structure: Scalar.
%                           Type: 'Single' / 'Double'.
%                           Range: {1, 2}.
% Output:
%   - mX                -   Circulant Matrix.
%                           The closest Circulant Matrix to mY in the
%                           Frobenius Norm meaning.
%                           Structure: Matrix.
%                           Type: 'Single' / 'Double'.
%                           Range: (-inf, inf).
% References:
%   1.  StackOverflow (2778195) - https://math.stackexchange.com/questions/2778195.
% Remarks:
%   1.  Prefixes:
%       -   'm' - Matrix.
%       -   'v' - Vector.
%   2.  T
% TODO:
%   1.  Accelerate calculation using FFT instead of the DFT Matrix.
%   Release Notes:
%   -   1.1.001     13/05/2018  Royi Avital
%       *   Accelerated using FFT instead of Matrix Multiplication.
%   -   1.1.000     13/05/2018  Royi Avital
%       *   Added support for Real / Complex output.
%   -   1.0.000     12/05/2018  Royi Avital
%       *   First release version.
% ----------------------------------------------------------------------------------------------- %

FALSE   = 0;
TRUE    = 1;

OFF = 0;
ON  = 1;

DOMAIN_FLAG_REAL    = 1;
DOMAIN_FLAG_COMPLEX = 2;

if(exist('domainFlag', 'var') == FALSE)
    % Default
    domainFlag = DOMAIN_FLAG_REAL;
end

% inputDim = size(mY, 1); %<! Matrix is Square Matrix

% mF = dftmtx(inputDim) / sqrt(inputDim); %<! Unitary DFT Matrix
% vY = diag(mF * mY * mF');
vY = diag(ifft( fft( mY, [], 1 ), [], 2 ));

% mX = mF' * diag(vY) * mF; %<! If mY is real will be very close to pure real
mX = fft( ifft( diag( vY ), [], 1 ), [], 2 ); %<! This is F' * diag(vY) * F

if(domainFlag == DOMAIN_FLAG_REAL)
    mX = real(mX);
end


end

