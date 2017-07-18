% Levinson Recursion Run Time Analysis
% Solving mT * vX = vY where mT is a Toeplitz Matrix.
% Remarks:
%   1.  Matrix is assumed to be n x n in size.
% TODO:
% 	1.  Add Eigen based implementation.
% Release Notes
% - 1.0.000     18/07/2017
%   *   First release.


%% General Parameters

run('InitScript.m');

figureIdx           = 0;
figureCounterSpec   = '%04d';

generateFigures = ON;

LIB_NAME = 'LevinsonRecursionToeplitzMatrix';
LIB_NAME = 'LevinsonRecursionToeplitzMatrixGcc';

if(~libisloaded(LIB_NAME))
    [a, b] = loadlibrary([LIB_NAME, '.dll'], 'LevinsonRecursionToeplitzMatrix.h');
    libfunctions(LIB_NAME, '-full');
end

SEC_TO_MILI_SEC_FCTR = 1000;


%% Simulation Parameters

numIterations   = 20;
vMatrixSize     = [1, 10, 25, 50, 75, 100, 200, 300, 400, 500, 1000, 2000, 5000];


%% Analysis

numMtxSize  = length(vMatrixSize);
% Run Time Matrix
% Column I - MATLAB Solver, Column II - Levinson Recursion MATLAB, Column III - Levinson Recursion C
mRunTime    = zeros([numMtxSize, 3]);

for ii = 1:numMtxSize
    
    mtxSize = vMatrixSize(ii);
    
    vR      = randn([mtxSize, 1]);
    vC      = randn([mtxSize, 1]);
    vC(1)   = vR(1);
    mT      = toeplitz(vC, vR);
    
    vY = rand(mtxSize, 1);
    
    hRunTime = tic();
    for jj = 1:numIterations
        vX = mT \ vY;
    end
    runTime = toc(hRunTime) / numIterations;
    
    mRunTime(ii, 1) = runTime;
    
    hRunTime = tic();
    for jj = 1:numIterations
        vX = LevinsonRecursion(mT, vY);
    end
    runTime = toc(hRunTime) / numIterations;
    
    mRunTime(ii, 2) = runTime;
    
    hRunTime = tic();
    for jj = 1:numIterations
        [~, ~, vX] = calllib(LIB_NAME, 'LevinsonRecursionToeplitzMatrix', mT, vY, vX, mtxSize);
    end
    runTime = toc(hRunTime) / numIterations;
    
    mRunTime(ii, 3) = runTime;
    
end


%% Dsipaly Results

figureIdx   = figureIdx + 1;
hFigure     = figure('Position', figPosMedium);
hAxes       = axes();
% set(hAxes, 'NextPlot', 'add');
hLineSeries = plot(vMatrixSize, SEC_TO_MILI_SEC_FCTR * mRunTime);
set(hLineSeries, 'LineWidth', lineWidthNormal);
set(get(hAxes, 'Title'), 'String', {['Linear Equation Solver Run Time - Toeplitz Matrix']}, ...
    'FontSize', fontSizeTitle);
set(get(hAxes, 'XLabel'), 'String', {['Matrix Dimension']}, ...
    'FontSize', fontSizeAxis);
set(get(hAxes, 'YLabel'), 'String', {['Run Time [Mili Sec]']}, ...
    'FontSize', fontSizeAxis);
hLegend = ClickableLegend({['MATLAB \\ Solver'], ['Levinson Recursion (MATLAB)'], ['Levinson Recursion (C)']});
set(hAxes, 'LooseInset', [0.05, 0.05, 0.05, 0.05]);

if(generateFigures == ON)
    saveas(hFigure,['Figure', num2str(figureIdx, figureCounterSpec), '.png']);
end


%% Restore Defaults

unloadlibrary(LIB_NAME);

% set(0, 'DefaultFigureWindowStyle', 'normal');
% set(0, 'DefaultAxesLooseInset', defaultLoosInset);

