% Prox of Total Variation (TV) Norm Analysis
% Analyzing the solution of different methods for calculation of the Prox
% of the Total Variation (TV) Norm. The problem given by:
% $$ \arg \min_{ x \in \mathbb{R}^{n} } \frac{1}{2} {\left\| x - y \right|}_{2}^{2} + \lambda {\left\| D x \right\|}_{1} $$
% References:
%   1.  aa
% Remarks:
%   1.  sa
% TODO:
% 	1.  ds
% Release Notes
% - 1.0.000     27/11/2019
%   *   First release.


%% General Parameters

subStreamNumberDefault = 79;

run('InitScript.m');

figureIdx           = 0; %<! Continue from Question 1
figureCounterSpec   = '%04d';

generateFigures = OFF;


%% Simulation Parameters

numElements = 40;
paramLambda = 0.5;

numIterations   = 5000;

cMethodString   = {['CVX'], ['Sub Gradient Method'], ['Chambolle Method'], ['ADMM']};
methodIdx       = 0;


%% Generate Data

vY = 10 * randn(numElements, 1);

% Generate the Diff Operator (1D Gradient) by Finite Differences
mD = spdiags([-ones(numElements, 1), ones(numElements, 1)], [0, 1], numElements - 1, numElements);

% Objective Function
hObjFun = @(vX) (0.5 * sum( (vX - vY) .^ 2)) + (paramLambda * sum(abs(mD * vX)));

numMethods  = size(cMethodString, 2);
vObjVal     = zeros(numMethods, 1);
mX          = zeros(numElements, numMethods);


%% Solution by CVX

methodIdx = methodIdx + 1;

cvx_begin('quiet')
    cvx_precision('best');
    variable vX(numElements);
    minimize( (0.5 * pow_pos(norm(vX - vY, 2), 2)) + (paramLambda * norm(mD * vX, 1)));
cvx_end

disp([' ']);
disp([cMethodString{methodIdx}, ' Solution Summary']);
disp(['The ', cMethodString{methodIdx}, ' Solver Status - ', cvx_status]);
disp(['The Optimal Value Is Given By - ', num2str(cvx_optval)]);
disp(['The Optimal Argument Is Given By - [ ', num2str(vX(:).'), ' ]']);
disp([' ']);

vObjVal(methodIdx)  = hObjFun(vX);
mX(:, methodIdx)    = vX;


%% Solution by Sub Gradient Descent
%{
Solving $ \arg \min_x \frac{1}{2} {\left\| x - y \right\|}_{2}^{2} +
{\lambda}_{1} {\left\| x \right\|}_{1} + {\lambda}_{2} {\left\| D x \right\|}_{1} $
%}

methodIdx = methodIdx + 1;

vX = SolveProxTvChambolle(vY, mD, paramLambda, numIterations);

disp([' ']);
disp([cMethodString{methodIdx}, ' Solution Summary']);
disp(['The Optimal Value Is Given By - ', num2str(hObjFun(vX))]);
disp(['The Optimal Argument Is Given By - [ ', num2str(vX(:).'), ' ]']);
disp([' ']);

vObjVal(methodIdx)  = hObjFun(vX);
mX(:, methodIdx)    = vX;


%% Solution by Chambolle's Method
%{
Solving
%}

methodIdx = methodIdx + 1;

vX = SolveProxTvChambolle(vY, mD, paramLambda, numIterations);

disp([' ']);
disp([cMethodString{methodIdx}, ' Solution Summary']);
disp(['The Optimal Value Is Given By - ', num2str(hObjFun(vX))]);
disp(['The Optimal Argument Is Given By - [ ', num2str(vX(:).'), ' ]']);
disp([' ']);

vObjVal(methodIdx)  = hObjFun(vX);
mX(:, methodIdx)    = vX;


%% Solution by ADMM
%{
Solving
%}

methodIdx = methodIdx + 1;

vX = SolveProxTvAdmm(vY, mD, paramLambda, numIterations);

disp([' ']);
disp([cMethodString{methodIdx}, ' Solution Summary']);
disp(['The Optimal Value Is Given By - ', num2str(hObjFun(vX))]);
disp(['The Optimal Argument Is Given By - [ ', num2str(vX(:).'), ' ]']);
disp([' ']);

vObjVal(methodIdx)  = hObjFun(vX);
mX(:, methodIdx)    = vX;


%% Analysis

vObjFunValErr   = abs(vObjVal(1) - vObjVal(2:end));
vObjFunValErr   = vObjFunValErr(:);
vSolErrNorm2    = vecnorm(mX(:, 1) - mX(:, 2:end), 2);
vSolErrNorm2    = vSolErrNorm2(:);


%% Display Results

figureIdx = figureIdx + 1;

hFigure     = figure('Position', figPosLarge);
hAxes       = axes();
hLineSeries = plot([1:(numMethods - 1)], 10 * log10([vObjFunValErr, vSolErrNorm2]));
set(hLineSeries, 'LineWidth', lineWidthNormal);
set(hLineSeries, 'LineStyle', 'none');
set(hLineSeries, 'Marker', '*');
set(get(hAxes, 'Title'), 'String', {['Prox Solvers for the Total Variation Norm'], ['SubStream - ', num2str(subStreamNumber)]}, ...
    'FontSize', fontSizeTitle);
set(get(hAxes, 'XLabel'), 'String', 'Method', ...
    'FontSize', fontSizeAxis);
set(hAxes, 'XTick', [1:(numMethods - 1)]);
set(hAxes, 'XTickLabel', cMethodString(2:end));
set(hAxes, 'XTickLabelRotation', 45);
set(get(hAxes, 'YLabel'), 'String', 'Value vs. CVX [dB]', ...
    'FontSize', fontSizeAxis);
hLegend = ClickableLegend({['$ 10 \log_{10} \left( \left| f \left( x \right) - f \left( {x}_{CVX} \right) \right| \right) $'], ...
    ['$ 10 \log_{10} \left( {\left\| x - {x}_{CVX} \right\|}_{2} \right) $']}, 'Interpreter', 'latex');
set(hAxes, 'LooseInset', [0.07, 0.07, 0.07, 0.07]);

if(generateFigures == ON)
    % saveas(hFigure,['Figure', num2str(figureIdx, figureCounterSpec), '.png']);
    print(hFigure, ['Figure', num2str(figureIdx, figureCounterSpec), '.png'], '-dpng', '-r0'); %<! Saves as Screen Resolution
end


%% Restore Defaults

% set(0, 'DefaultFigureWindowStyle', 'normal');
% set(0, 'DefaultAxesLooseInset', defaultLoosInset);

