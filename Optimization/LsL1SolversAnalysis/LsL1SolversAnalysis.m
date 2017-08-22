% Mathematics Q561696
% https://math.stackexchange.com/questions/561696
% Solving Non Negative Least Squares by Analogy with Least Squares (MATLAB)
% References:
%   1.  aa
% Remarks:
%   1.  sa
% TODO:
% 	1.  ds
% Release Notes
% - 1.0.000     05/08/2017
%   *   First release.


%% General Parameters

run('InitScript.m');

figureIdx           = 0; %<! Continue from Question 1
figureCounterSpec   = '%04d';

generateFigures = OFF;


%% Simulation Parameters

numRows = 100;
numCols = 40; %<! Number of Vectors - i (K in the question)

paramLambda = 0.1;

numIterations   = 5000;
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


%% Solution by Proximal Gradient Method

[vX, mX] = SolveLsL1Prox(mA, vB, paramLambda, numIterations);
% [vX, mX] = SolveLsL1ProxAccel(mA, vB, paramLambda, numIterations);

objVal = hObjFun(vX);

disp([' ']);
disp(['Projected Proximal Gradient Method Solution Summary']);
disp(['The Optimal Value Is Given By - ', num2str(objVal)]);
disp(['The Optimal Argument Is Given By - [ ', num2str(vX.'), ' ]']);
disp([' ']);

solverIdx                   = solverIdx + 1;
cLegendString{solverIdx}    = ['Proximal Gradient'];

for ii = 1:numIterations

    mObjFunVal(ii, solverIdx)   = abs(hObjFun(mX(:, ii)) - cvx_optval);
    mSolErrNorm(ii, solverIdx)  = norm(mX(:, ii) - vXCvx);

end


%% Solution by Proximal Gradient Method with Line Search

[vX, mX] = SolveLsL1ProxLs(mA, vB, paramLambda, numIterations);

objVal = hObjFun(vX);

disp([' ']);
disp(['Projected Proximal Gradient Method (LS) Solution Summary']);
disp(['The Optimal Value Is Given By - ', num2str(objVal)]);
disp(['The Optimal Argument Is Given By - [ ', num2str(vX.'), ' ]']);
disp([' ']);

solverIdx                   = solverIdx + 1;
cLegendString{solverIdx}    = ['Proximal Gradient - Line Search'];

for ii = 1:numIterations

    mObjFunVal(ii, solverIdx)   = abs(hObjFun(mX(:, ii)) - cvx_optval);
    mSolErrNorm(ii, solverIdx)  = norm(mX(:, ii) - vXCvx);

end


%% Solution by Proximal Gradient Method Accelerated

[vX, mX] = SolveLsL1ProxAccel(mA, vB, paramLambda, numIterations);

objVal = hObjFun(vX);

disp([' ']);
disp(['Projected Proximal Gradient Method (Accelerated) Solution Summary']);
disp(['The Optimal Value Is Given By - ', num2str(objVal)]);
disp(['The Optimal Argument Is Given By - [ ', num2str(vX.'), ' ]']);
disp([' ']);

solverIdx                   = solverIdx + 1;
cLegendString{solverIdx}    = ['Proximal Gradient - Accelerated'];

for ii = 1:numIterations

    mObjFunVal(ii, solverIdx)   = abs(hObjFun(mX(:, ii)) - cvx_optval);
    mSolErrNorm(ii, solverIdx)  = norm(mX(:, ii) - vXCvx);

end


%% Solution by Proximal Gradient Method Accelerated with Line Search

[vX, mX] = SolveLsL1ProxAccelLs(mA, vB, paramLambda, numIterations);

objVal = hObjFun(vX);

disp([' ']);
disp(['Projected Proximal Gradient Method (Accelerated + LS) Solution Summary']);
disp(['The Optimal Value Is Given By - ', num2str(objVal)]);
disp(['The Optimal Argument Is Given By - [ ', num2str(vX.'), ' ]']);
disp([' ']);

solverIdx                   = solverIdx + 1;
cLegendString{solverIdx}    = ['Proximal Gradient - Accelerated + Line Search'];

for ii = 1:numIterations

    mObjFunVal(ii, solverIdx)   = abs(hObjFun(mX(:, ii)) - cvx_optval);
    mSolErrNorm(ii, solverIdx)  = norm(mX(:, ii) - vXCvx);

end


%% Display Results

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

