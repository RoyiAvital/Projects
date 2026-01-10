function [ mM ] = GenMonoOp( vI, vY )
% ----------------------------------------------------------------------------------------------- %
% [ mM ] = GenMonoOp( vI, vY )
%   Generates the Matrix Operator which represents a piece wise monotonic 
%   property by linear inequality: M * x <= 0.
% Input:
%   - vI            -   Segments Indices Vector.
%                       Each vI(k) and vI(k + 1) define a segment which
%                       should be monotonic.
%                       Structure: Vector (Column).
%                       Type: 'Single' / 'Double'.
%                       Range: {1, 2, ..., }.
%   - vY            -   Values Vector.
%                       Values to set the monotonicity property per
%                       segment, namely decreasing or increasing.
%                       Structure: Vector (Column).
%                       Type: 'Single' / 'Double'.
%                       Range: (-inf, inf).
% Output:
%   - mM            -   Monotonic Operator.
%                       Each row forces decreasing or increasing property
%                       on a sample.
%                       Structure: Matrix (Sparse).
%                       Type: 'Single' / 'Double'.
%                       Range: {-1, 0, 1}.
% References
%   1.  A
% Remarks:
%   1.  It is assumed `vI` is sorted.
% TODO:
%   1.  C.
% Release Notes:
%   -   1.0.000     09/01/2026  Royi Avital
%       *   First release version.
% ----------------------------------------------------------------------------------------------- %

arguments(Input)
    vI (:, 1) {mustBeNumeric, mustBeFinite, mustBeInteger, mustBePositive}
    vY (:, 1) {mustBeNumeric, mustBeFinite, mustBeReal}
end

arguments(Output)
    mM (:, :) {mustBeNumeric, mustBeFinite, mustBeSparse}
end

numRefPts  = length(vI);
numSamples = length(vY);

numSegments = numRefPts - 1;

% Default Monotonic Non Decreasing
% Forcing the monotonicity by the relationship with the next sample
mD = ones(numSamples, 2); %<! Current sample
mD(:, 2) = -1; %<! Next sample

for ii = 1:numSegments
    startIdx = vI(ii);
    endIdx   = vI(ii + 1);

    if vY(startIdx) <= vY(endIdx)
        % Monotonic Non Decreasing
        valSign = 1;
    else
        % Monotonic Non Increasing
        valSign = -1;
    end

    for jj = startIdx:(endIdx - 1)
        mD(jj, :) = valSign * mD(jj, :);
    end

end

% MATLAB's `spdiags()` takes, for positive k diagonal, from the `k` sample.
% Shifting to "align" the coefficients in `mD`.
mD(:, 2) = circshift(mD(:, 2), 1);

mM = spdiags(mD, [0, 1], numSamples, numSamples);
mM = mM(vI(1):(vI(end) - 1), :);


end

