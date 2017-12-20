% Least Squares Regularized by L0 Norm - Solvers Analysis
% Solving: \arg \min_{x} 0.5 * || A x - b ||_{2}^{2} + \lambda || x ||_{0}
% References:
%   1.  aa
% Remarks:
%   1.  
% Known Issues:
%   1.  IRLS Solver doesn't work.
% TODO:
% 	1.  Add Coordinate Descent based algorithm.
% Release Notes
% - 1.0.000     08/12/2017
%   *   First release.


%% General Parameters

run('InitScript.m');

figureIdx           = 0; %<! Continue from Question 1
figureCounterSpec   = '%04d';

generateFigures = OFF;


%% Simulation Parameters

numRows = 6;
numCols = 25; %<! Number of Vectors - i (K in the question)

paramLambda = 0.3;
tolVal      = 1e-6;

numIterations   = 500;
stopThr         = 0;

solverIdx       = 0;
cLegendString   = {};

mObjFunVal  = zeros([numIterations, 1]);
mSolErrNorm = zeros([numIterations, 1]);


%% Generate Data

mA = randn([numRows, numCols]);
vB = randn([numRows, 1]);

hL0Norm = @(vX) sum(abs(vX) > tolVal);
hObjFun = @(vX) (0.5 * sum(((mA * vX) - vB) .^ 2)) + (paramLambda * hL0Norm(vX));


%% Solution by Matching Pursuit

[vX, mX] = SolveLsL0Mp(mA, vB, paramLambda, numIterations, tolVal);

objVal = hObjFun(vX);

disp([' ']);
disp(['MP Solution Summary']);
disp(['The Optimal Value Is Given By - ', num2str(objVal)]);
disp(['The Optimal Argument L0 Norm - ', num2str(hL0Norm(vX))]);
disp(['The Optimal Argument Is Given By - [ ', num2str(vX.'), ' ]']);
disp([' ']);

vIdx = find(vX);


%% Solution by Orthogonal Matching Pursuit

[vX, mX] = SolveLsL0Omp(mA, vB, paramLambda, numIterations, tolVal);

objVal = hObjFun(vX);

disp([' ']);
disp(['OMP Solution Summary']);
disp(['The Optimal Value Is Given By - ', num2str(objVal)]);
disp(['The Optimal Argument L0 Norm - ', num2str(hL0Norm(vX))]);
disp(['The Optimal Argument Is Given By - [ ', num2str(vX.'), ' ]']);
disp([' ']);


%% Display Results


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

