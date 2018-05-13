% Project onto the Circulant Matrix Set (Real) Analysis
% 
% References:
%   1.  Projection onto the Set of Circulant Matrices - https://math.stackexchange.com/questions/2778195.
% Remarks:
%   1.  sa
% TODO:
% 	1.  ds
% Release Notes
% - 1.0.000     12/05/2018
%   *   First release.


%% General Parameters

run('InitScript.m');

figureIdx           = 0;
figureCounterSpec   = '%04d';

generateFigures = ON;


%% Parameters

matDim = 5;


%% Generate Data

mY = randn([matDim, matDim]);
mPi = eye(matDim);
mPi = [mPi(2:end, :); mPi(1, :)]; %<! Forward Shift Matrix

mF = dftmtx(matDim) / sqrt(matDim);


%% Solution by Reference

vC = zeros([matDim, 1]);
mC = zeros([matDim, matDim]);

for ii = 1:matDim
    mPi = GenerateShiftMatrix(matDim, ii - 1);
    % vC(ii) = trace(mY * (mPi ^ (ii - 1)).') / matDim;
    % vC(ii) = trace(mY.' * (mPi ^ (ii - 1))) / matDim;
    % vC(ii) = trace(mY.' * mPi) / matDim;
    
    vC(ii) = (mY(:).' * mPi(:)) / matDim;
    
    % mC = mC + (vC(ii) * (mPi ^ (ii - 1)));
    mC = mC + (vC(ii) * mPi);
end

mX = ProjectCirculantMatrixSet(mY);

% mX = mF' * diag(diag(mF * mY * mF')) * mF

disp(['The Reference Frobenius Norm - ', num2str(norm(mX - mY, 'fro'))]);
disp(['The Suggested Frobenius Norm - ', num2str(norm(mC - mY, 'fro'))]);


%% Restore Defaults

% set(0, 'DefaultFigureWindowStyle', 'normal');
% set(0, 'DefaultAxesLooseInset', defaultLoosInset);

