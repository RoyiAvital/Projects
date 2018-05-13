% Least Squares Regularized by L1 Norm - Solvers Analysis
% Solving: \arg \min_{x} 0.5 * || A x - b ||_{2}^{2} + \lambda || x ||_{1}
% References:
%   1.  aa
% Remarks:
%   1.  
% Known Issues:
%   1.  IRLS Solver doesn't work.
% TODO:
% 	1.  Add Coordinate Descent based algorithm.
% Release Notes
% - 1.0.000     23/08/2017
%   *   First release.


%% General Parameters

run('InitScript.m');

figureIdx           = 0;
figureCounterSpec   = '%04d';

generateFigures = ON;


%% Simulation Parameters

numRows = 200;
numCols = 50; %<! Number of Vectors - i (K in the question)

paramLambda = 1.1;

numIterations   = 500;
stopThr         = 0;

solverIdx       = 0;
cLegendString   = {};

mObjFunVal  = zeros([numIterations, 1]);
mSolErrNorm = zeros([numIterations, 1]);


%% Generate Data

mA = randn([numRows, numCols]);
vB = randn([numRows, 1]);

hObjFun = @(vX) (0.5 * sum(((mA * vX) - vB) .^ 2)) + (paramLambda * sum(abs(vX)));


%% Solution by CVX

cvx_begin('quiet')
    cvx_precision('best');
    variable vXCvx(numCols)
    minimize( (0.5 * sum_square(mA * vXCvx - vB)) + (paramLambda * norm(vXCvx, 1)) )
cvx_end

disp([' ']);
disp(['CVX Solution Summary']);
disp(['The CVX Solver Status - ', cvx_status]);
disp(['The Optimal Value Is Given By - ', num2str(cvx_optval)]);
disp(['The Optimal Argument Is Given By - [ ', num2str(vXCvx.'), ' ]']);
disp([' ']);

sCvxSol.vXCvx     = vXCvx;
sCvxSol.cvxOptVal = cvx_optval;


%% Solution by Sub Gradient Method

[vX, mX] = SolveLsL1SubGradient(mA, vB, paramLambda, numIterations);

objVal = hObjFun(vX);

disp([' ']);
disp(['Sub Gradient Method Solution Summary']);
disp(['The Optimal Value Is Given By - ', num2str(objVal)]);
disp(['The Optimal Argument Is Given By - [ ', num2str(vX.'), ' ]']);
disp([' ']);

solverIdx                   = solverIdx + 1;
cLegendString{solverIdx}    = ['Sub Gradient'];

[mObjFunVal, mSolErrNorm] = UpdateAnalysisData(mObjFunVal, mSolErrNorm, mX, hObjFun, sCvxSol, solverIdx);


%% Solution by Smoothing (Huber Loss) Method

[vX, mX] = SolveLsL1Huber(mA, vB, paramLambda, numIterations);

objVal = hObjFun(vX);

disp([' ']);
disp(['Smoothing (Huber Loss) Method Solution Summary']);
disp(['The Optimal Value Is Given By - ', num2str(objVal)]);
disp(['The Optimal Argument Is Given By - [ ', num2str(vX.'), ' ]']);
disp([' ']);

solverIdx                   = solverIdx + 1;
cLegendString{solverIdx}    = ['Smoothing (Huber Loss)'];

[mObjFunVal, mSolErrNorm] = UpdateAnalysisData(mObjFunVal, mSolErrNorm, mX, hObjFun, sCvxSol, solverIdx);


%% Solution by Proximal Gradient Method

[vX, mX] = SolveLsL1Prox(mA, vB, paramLambda, numIterations);

objVal = hObjFun(vX);

disp([' ']);
disp(['Proximal Gradient Method Solution Summary']);
disp(['The Optimal Value Is Given By - ', num2str(objVal)]);
disp(['The Optimal Argument Is Given By - [ ', num2str(vX.'), ' ]']);
disp([' ']);

solverIdx                   = solverIdx + 1;
cLegendString{solverIdx}    = ['Proximal Gradient'];

[mObjFunVal, mSolErrNorm] = UpdateAnalysisData(mObjFunVal, mSolErrNorm, mX, hObjFun, sCvxSol, solverIdx);


%% Solution by Proximal Gradient Method with Line Search

[vX, mX] = SolveLsL1ProxLs(mA, vB, paramLambda, numIterations);

objVal = hObjFun(vX);

disp([' ']);
disp(['Proximal Gradient Method (LS) Solution Summary']);
disp(['The Optimal Value Is Given By - ', num2str(objVal)]);
disp(['The Optimal Argument Is Given By - [ ', num2str(vX.'), ' ]']);
disp([' ']);

solverIdx                   = solverIdx + 1;
cLegendString{solverIdx}    = ['Proximal Gradient - Line Search'];

[mObjFunVal, mSolErrNorm] = UpdateAnalysisData(mObjFunVal, mSolErrNorm, mX, hObjFun, sCvxSol, solverIdx);


%% Solution by Proximal Gradient Method Accelerated

[vX, mX] = SolveLsL1ProxAccel(mA, vB, paramLambda, numIterations);

objVal = hObjFun(vX);

disp([' ']);
disp(['Proximal Gradient Method (Accelerated) Solution Summary']);
disp(['The Optimal Value Is Given By - ', num2str(objVal)]);
disp(['The Optimal Argument Is Given By - [ ', num2str(vX.'), ' ]']);
disp([' ']);

solverIdx                   = solverIdx + 1;
cLegendString{solverIdx}    = ['Proximal Gradient - Accelerated'];

[mObjFunVal, mSolErrNorm] = UpdateAnalysisData(mObjFunVal, mSolErrNorm, mX, hObjFun, sCvxSol, solverIdx);


%% Solution by Proximal Gradient Method Accelerated with Line Search

[vX, mX] = SolveLsL1ProxAccelLs(mA, vB, paramLambda, numIterations);

objVal = hObjFun(vX);

disp([' ']);
disp(['Proximal Gradient Method (Accelerated + LS) Solution Summary']);
disp(['The Optimal Value Is Given By - ', num2str(objVal)]);
disp(['The Optimal Argument Is Given By - [ ', num2str(vX.'), ' ]']);
disp([' ']);

solverIdx                   = solverIdx + 1;
cLegendString{solverIdx}    = ['Proximal Gradient - Accelerated + Line Search'];

[mObjFunVal, mSolErrNorm] = UpdateAnalysisData(mObjFunVal, mSolErrNorm, mX, hObjFun, sCvxSol, solverIdx);


%% Solution by ADMM

[vX, mX] = SolveLsL1Admm(mA, vB, paramLambda, numIterations);

objVal = hObjFun(vX);

disp([' ']);
disp(['ADMM Solution Summary']);
disp(['The Optimal Value Is Given By - ', num2str(objVal)]);
disp(['The Optimal Argument Is Given By - [ ', num2str(vX.'), ' ]']);
disp([' ']);

solverIdx                   = solverIdx + 1;
cLegendString{solverIdx}    = ['ADMM'];

[mObjFunVal, mSolErrNorm] = UpdateAnalysisData(mObjFunVal, mSolErrNorm, mX, hObjFun, sCvxSol, solverIdx);


%% Solution by IRLS

[vX, mX] = SolveLsL1Irls(mA, vB, paramLambda, numIterations);

objVal = hObjFun(vX);

disp([' ']);
disp(['IRLS Solution Summary']);
disp(['The Optimal Value Is Given By - ', num2str(objVal)]);
disp(['The Optimal Argument Is Given By - [ ', num2str(vX.'), ' ]']);
disp([' ']);

solverIdx                   = solverIdx + 1;
cLegendString{solverIdx}    = ['IRLS'];

[mObjFunVal, mSolErrNorm] = UpdateAnalysisData(mObjFunVal, mSolErrNorm, mX, hObjFun, sCvxSol, solverIdx);


%% Solution by CD

[vX, mX] = SolveLsL1Cd(mA, vB, paramLambda, numIterations);

objVal = hObjFun(vX);

disp([' ']);
disp(['CD Solution Summary']);
disp(['The Optimal Value Is Given By - ', num2str(objVal)]);
disp(['The Optimal Argument Is Given By - [ ', num2str(vX.'), ' ]']);
disp([' ']);

solverIdx                   = solverIdx + 1;
cLegendString{solverIdx}    = ['CD'];

[mObjFunVal, mSolErrNorm] = UpdateAnalysisData(mObjFunVal, mSolErrNorm, mX, hObjFun, sCvxSol, solverIdx);


%% Display Results

figureIdx = figureIdx + 1;

hFigure     = figure('Position', figPosLarge);

hAxes       = subplot_tight(2, 1, 1, [0.09, 0.09]);
hLineSeries = plot(1:numIterations, 10 * log10(mObjFunVal));
set(hLineSeries, 'LineWidth', lineWidthNormal);
set(get(hAxes, 'Title'), 'String', ['Objective Function Value vs. Optimal Value (CVX)'], ...
    'FontSize', fontSizeTitle);
set(get(hAxes, 'XLabel'), 'String', 'Iteration Number', ...
    'FontSize', fontSizeAxis);
set(get(hAxes, 'YLabel'), 'String', 'Objective Function Value', ...
    'FontSize', fontSizeAxis);
set(hAxes, 'XLim', [1, numIterations]);
hLegend = ClickableLegend(cLegendString);

hAxes       = subplot_tight(2, 1, 2, [0.09, 0.09]);
hLineSeries = plot(1:numIterations, 10 * log10(mSolErrNorm));
set(hLineSeries, 'LineWidth', lineWidthNormal);
set(get(hAxes, 'Title'), 'String', ['Solution Error Norm'], ...
    'FontSize', fontSizeTitle);
set(get(hAxes, 'XLabel'), 'String', 'Iteration Number', ...
    'FontSize', fontSizeAxis);
set(get(hAxes, 'YLabel'), 'String', 'Objective Function Value', ...
    'FontSize', fontSizeAxis);
set(hAxes, 'XLim', [1, numIterations]);
hLegend = ClickableLegend(cLegendString);

if(generateFigures == ON)
    saveas(hFigure,['Figure', num2str(figureIdx, figureCounterSpec), '.png']);
end


%% Restore Defaults

% set(0, 'DefaultFigureWindowStyle', 'normal');
% set(0, 'DefaultAxesLooseInset', defaultLoosInset);

