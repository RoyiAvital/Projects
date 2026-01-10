% Curve Smoothing
% Analysis of curve smoothing algorithms.
% References:
%   1.  A
% Remarks:
%   1.  B
% TODO:
% 	1.  C
% Release Notes Royi Avital RoyiAvital@yahoo.com
% - 1.0.000     09/01/2026
%   *   First release.


%% General Parameters

subStreamNumberDefault = 79;

run('InitScript.m');

figureIdx           = 0;
figureCounterSpec   = '%04d';

generateFigures = ON;


%% Constants


%% Parameters

% Data
% From https://huggingface.co/datasets/thuml/Time-Series-Library (https://huggingface.co/datasets/thuml/Time-Series-Library/tree/main/exchange_rate).
csvFileName = 'exchange_rate.csv';
varName     = 'Var2';
decFactor   = 50;

% Model
diffPow       = 2;
paramLambda   = 1.95;
numRefPtsFctr = 0.05;


%% Generate / Load Data

tT = readtable(csvFileName);
vY = tT.(varName)(2:decFactor:end);

numSamples = length(vY);
vT = 1:numSamples;

mD = GenDiffOp(diffPow, numSamples);

numRefPts = round(numRefPtsFctr * numSamples);
vI = sort(randperm(numSamples, numRefPts));

mA = GenMonoOp(vI, vY);

mE = sparse(1:numRefPts, vI, 1, numRefPts, numSamples);


%% Analysis

mH = speye(numSamples) + paramLambda * (mD.' * mD);
vF = -vY;
sOpt = optimoptions('quadprog', 'Display', 'off');
vXRef = quadprog(mH, vF, mA, zeros(size(mA, 1), 1), mE, vY(vI), [], [], [], sOpt);

[vX, isConv] = SplineQPSmooth(vY, mD, paramLambda, vI, mA);

fprintf('Converged: %d\n', isConv);

hF = @() SplineQPSmoothQuadProg(vY, mD, paramLambda, vI, mA);
runTime = TimeItMin(hF);
fprintf('The runtime of `quadprog(): %0.3f [Mili Sec]\n', runTime * 1000);

hF = @() SplineQPSmooth(vY, mD, paramLambda, vI, mA);
runTime = TimeItMin(hF);
fprintf('The runtime of SplineQPSmooth(): %0.3f [Mili Sec]\n', runTime * 1000);


%% Display Results

hF = figure();
hA = axes();
set(hA, 'NextPlot', 'add');
hP = plot(vT, vY, 'DisplayName', 'Data Samples');
set(hP, 'LineWidth', lineWidthNormal);
hP = plot(vT(vI), vY(vI), 'DisplayName', 'Refernce Samples');
set(hP, 'LineStyle', 'none', 'LineWidth', lineWidthNormal, 'Marker', 'o', 'MarkerSize', 8);
hP = plot(vT, vXRef, 'DisplayName', 'Reference Solution');
set(hP, 'LineWidth', lineWidthThin, 'LineStyle', '-.');
hP = plot(vT, vX, 'DisplayName', 'Manual Solution');
set(hP, 'LineWidth', lineWidthThin, 'LineStyle', '-.');
ClickableLegend();


%% Auxiliary Functions

function [ vX ] = SplineQPSmoothQuadProg( vY, mD, paramLambda, vI, mA )

numSamples = length(vY);
numRefPts  = length(vI);

mE = sparse(1:numRefPts, vI, 1, numRefPts, numSamples);

mH = speye(numSamples) + paramLambda * (mD.' * mD);
vF = -vY;
sOpt = optimoptions('quadprog', 'Display', 'off');
vX = quadprog(mH, vF, mA, zeros(size(mA, 1), 1), mE, vY(vI), [], [], [], sOpt);

end



%% Restore Defaults

% set(0, 'DefaultFigureWindowStyle', 'normal');
% set(0, 'DefaultAxesLooseInset', defaultLoosInset);

