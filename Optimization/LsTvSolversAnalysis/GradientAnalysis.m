% Gradient Analysis
% References:
%   1.  aa
% Remarks:
%   1.  
% Known Issues:
%   1.  IRLS Solver doesn't work.
% TODO:
% 	1.  Add Coordinate Descent based algorithm.
% Release Notes
% - 1.0.000     27/11/2019
%   *   First release.


%% General Parameters

subStreamNumberDefault = 0;

run('InitScript.m');

figureIdx           = 0;
figureCounterSpec   = '%04d';

generateFigures = ON;

DIFF_MODE_FORWARD   = 1;
DIFF_MODE_BACKWARD  = 2;
DIFF_MODE_CENTRAL   = 3;
DIFF_MODE_COMPLEX   = 4;


%% Simulation Parameters

numRows = 10;
numCols = 5;

paramLambda = 0.1;

diffMode    = DIFF_MODE_CENTRAL;
epsVal      = 1e-7;


%% Generate Data

mA = randn(numRows, numCols);
% mA = eye(numRows);
vB = randn(numRows, 1);

vP = randn(numCols - 1, 1);

% Generate the Diff Operator (1D Gradient) by Finite Differences
mD = spdiags([-ones(numCols, 1), ones(numCols, 1)], [0, 1], numCols - 1, numCols);

vAB = mA.' * vB;
mAAInv = inv(mA.' * mA); %<! Useful Constant
vC = paramLambda * mD * mAAInv * vAB; %<! Useful Constant
mC = (paramLambda * paramLambda) * mD * mAAInv * mD.'; %<! Useful Constant

vXInit = pinv(mA) * vB;

hCalcXFun = @(vP) mAAInv * (vAB - paramLambda * mD.' * vP); %<! Calculates vX out of vP

% hObjFun = @(vX) (0.5 * sum(((mA * vX) - vB) .^ 2)) + (paramLambda * sum(abs(mD * vX)));
% hObjFun = @(vP) (0.5 * sum(((mA * hCalcXFun(vP)) - vB) .^ 2)) + (paramLambda * sum(abs(mD * hCalcXFun(vP))));

hObjFun = @(vP) (0.5 * sum(((mA * hCalcXFun(vP)) - vB) .^ 2)) + (paramLambda * vP.' * mD * hCalcXFun(vP));

mT0 = inv(mA.' * mA);
vT1 = mD.' * vP;
vT2 = mT0 * (mA.' * vB - paramLambda * vT1);


%% Numerical Solution

vGN = CalcFunGrad(vP, hObjFun, diffMode, epsVal);


%% Analytical Solution

vGA = (mC * vP) - vC; %<! Gradient of vP
vGA = -vGA;

% vGA = (paramLambda * mD * vT2) - (paramLambda * mD * mT0 * mA.' * (mA * vT2 - vB)) - (paramLambda * paramLambda * mD * mT0 * vT1);


%% Results Analysis

[vGN, vGA]

max(abs(vGN - vGA))


%% Restore Defaults

% set(0, 'DefaultFigureWindowStyle', 'normal');
% set(0, 'DefaultAxesLooseInset', defaultLoosInset);

