% Least Squares Regularized by L1 Norm - Solvers Analysis
% Solving: \arg \min_{x} 0.5 * || A x - b ||_{2}^{2} + \lambda || x ||_{TV}
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

subStreamNumberDefault = 79;

run('InitScript.m');

figureIdx           = 0;
figureCounterSpec   = '%04d';

generateFigures = ON;


%% Simulation Parameters

numRows = 20;
numCols = 5;

paramLambda = 0.35;

numIterations   = 500;
stopThr         = 0;

solverIdx       = 0;
cLegendString   = {};

mObjFunVal  = zeros([numIterations, 1]);
mSolErrNorm = zeros([numIterations, 1]);


%% Generate Data

mA = randn([numRows, numCols]);
% mA = eye(numRows);
vB = randn([numRows, 1]);

% Generate the Diff Operator (1D Gradient) by Finite Differences
mD = spdiags([-ones(numCols, 1), ones(numCols, 1)], [0, 1], numCols - 1, numCols);

vXInit = pinv(mA) * vB;

hObjFun = @(vX) (0.5 * sum(((mA * vX) - vB) .^ 2)) + (paramLambda * sum(abs(mD * vX)));


%% Solution by CVX

cvx_begin('quiet')
    cvx_precision('best');
    variable vXCvx(numCols)
    minimize( (0.5 * sum_square(mA * vXCvx - vB)) + (paramLambda * norm(mD * vXCvx, 1)) )
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

[vX, mX] = SolveLsTvSubGradient(vXInit, mA, vB, mD, paramLambda, numIterations);

objVal = hObjFun(vX);

disp([' ']);
disp(['Sub Gradient Method Solution Summary']);
disp(['The Optimal Value Is Given By - ', num2str(objVal)]);
disp(['The Optimal Argument Is Given By - [ ', num2str(vX.'), ' ]']);
disp([' ']);

solverIdx                   = solverIdx + 1;
cLegendString{solverIdx}    = ['Sub Gradient'];

[mObjFunVal, mSolErrNorm] = UpdateAnalysisData(mObjFunVal, mSolErrNorm, mX, hObjFun, sCvxSol, solverIdx);


%% Solution by Chambolle Method

[vX, mX] = SolveLsTvChambolle(vXInit, mA, vB, mD, paramLambda, numIterations);

objVal = hObjFun(vX);

disp([' ']);
disp(['Chambolle Method Solution Summary']);
disp(['The Optimal Value Is Given By - ', num2str(objVal)]);
disp(['The Optimal Argument Is Given By - [ ', num2str(vX.'), ' ]']);
disp([' ']);

solverIdx                   = solverIdx + 1;
cLegendString{solverIdx}    = ['Chambolle'];

[mObjFunVal, mSolErrNorm] = UpdateAnalysisData(mObjFunVal, mSolErrNorm, mX, hObjFun, sCvxSol, solverIdx);


%% Solution by ADMM

[vX, mX] = SolveLsTvAdmm(vXInit, mA, vB, mD, paramLambda, numIterations);

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

% [vX, mX] = SolveLsL1Irls(mA, vB, paramLambda, numIterations);
% 
% objVal = hObjFun(vX);
% 
% disp([' ']);
% disp(['IRLS Solution Summary']);
% disp(['The Optimal Value Is Given By - ', num2str(objVal)]);
% disp(['The Optimal Argument Is Given By - [ ', num2str(vX.'), ' ]']);
% disp([' ']);
% 
% solverIdx                   = solverIdx + 1;
% cLegendString{solverIdx}    = ['IRLS'];
% 
% [mObjFunVal, mSolErrNorm] = UpdateAnalysisData(mObjFunVal, mSolErrNorm, mX, hObjFun, sCvxSol, solverIdx);


%% Solution by CD

% [vX, mX] = SolveLsL1Cd(mA, vB, paramLambda, numIterations);
% 
% objVal = hObjFun(vX);
% 
% disp([' ']);
% disp(['CD Solution Summary']);
% disp(['The Optimal Value Is Given By - ', num2str(objVal)]);
% disp(['The Optimal Argument Is Given By - [ ', num2str(vX.'), ' ]']);
% disp([' ']);
% 
% solverIdx                   = solverIdx + 1;
% cLegendString{solverIdx}    = ['CD'];
% 
% [mObjFunVal, mSolErrNorm] = UpdateAnalysisData(mObjFunVal, mSolErrNorm, mX, hObjFun, sCvxSol, solverIdx);


%% Display Results

figureIdx = figureIdx + 1;

hFigure     = figure('Position', figPosXLarge);

hAxes       = subplot_tight(2, 1, 1, [0.09, 0.09]);
hLineSeries = plot(1:numIterations, 10 * log10(mObjFunVal));
set(hLineSeries, 'LineWidth', lineWidthNormal);
set(get(hAxes, 'Title'), 'String', {['Objective Function Value vs. Optimal Value (CVX)'], ['SubStream Number - ', num2str(subStreamNumber)]}, ...
    'FontSize', fontSizeTitle);
set(get(hAxes, 'XLabel'), 'String', 'Iteration Number', ...
    'FontSize', fontSizeAxis);
set(get(hAxes, 'YLabel'), 'String', '$ 10 \log_{10} \left( \left| f \left( x \right) - f \left( {x}_{CVX} \right) \right| \right) $', ...
    'FontSize', fontSizeAxis, 'Interpreter', 'latex');
set(hAxes, 'XLim', [1, numIterations]);
hLegend = ClickableLegend(cLegendString);

hAxes       = subplot_tight(2, 1, 2, [0.09, 0.09]);
hLineSeries = plot(1:numIterations, 10 * log10(mSolErrNorm));
set(hLineSeries, 'LineWidth', lineWidthNormal);
set(get(hAxes, 'Title'), 'String', ['Solution Error Norm'], ...
    'FontSize', fontSizeTitle);
set(get(hAxes, 'XLabel'), 'String', 'Iteration Number', ...
    'FontSize', fontSizeAxis);
set(get(hAxes, 'YLabel'), 'String', '$ 10 \log_{10} \left( {\left\| x - {x}_{CVX} \right\|}_{1} \right) $', ...
    'FontSize', fontSizeAxis, 'Interpreter', 'latex');
set(hAxes, 'XLim', [1, numIterations]);
hLegend = ClickableLegend(cLegendString);

if(generateFigures == ON)
    % saveas(hFigure, ['Figure', num2str(figureIdx, figureCounterSpec), '.png']);
    print(hFigure, ['Figure', num2str(figureIdx, figureCounterSpec), '.png'], '-dpng', '-r0'); %<! Saves as Screen Resolution
end


%% Restore Defaults

% set(0, 'DefaultFigureWindowStyle', 'normal');
% set(0, 'DefaultAxesLooseInset', defaultLoosInset);

