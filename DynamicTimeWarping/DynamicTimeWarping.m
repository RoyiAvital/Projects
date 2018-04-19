% Dynamic Time Warping
% Display of the Dynamic Time Warping concept.
% References:
%   1.  aa
% Remarks:
%   1.  sa
% TODO:
% 	1.  ds
% Release Notes
% - 1.0.000     09/05/2016
%   *   First release.


%% General Parameters and Initialization

run('InitScript.m');

figureIdx           = 0;
figureCounterSpec   = '%04d';

generateFigures = OFF;


%%

vTimeSupport1 = [0:0.1:1].';
vTimeSupport2 = [0:0.2:1.6].';


vY1 = sin(2 * pi * vTimeSupport1);
vY2 = sin(2 * pi * vTimeSupport2);

figure();
plot(vTimeSupport1, vY1, vTimeSupport2, vY2);


mDistMtx = bsxfun(@minus, vY1, vY2.') .^ 2;

numRows = size(mDistMtx, 1);
numCols = size(mDistMtx, 2);


% hFigure     = figure();
% hAxes       = axes();
% DisplayMatrixValues(mDistMtx, hAxes);
% 
% set(hAxes, 'DataAspectRatio', [1, 1, 1])


figure('Position', [100, 100, 900, 900]);

hAxes(1) = subplot(3, 3, [2, 3, 5, 6]);
hImageObj   = DisplayMatrixValues(mDistMtx, hAxes(1));
% set(hAxes(1), 'DataAspectRatio', [1, 1, 1]);
set(get(hAxes(1), 'Title'), 'String', ['Distance Matrix'], ...
    'FontSize', fontSizeTitle);

hAxes(2) = subplot(3, 3, [1, 4]);
hLineSeries   = plot(vY1);
set(hLineSeries, 'Marker', '*', 'MarkerSize', markerSizeMedium, ...
    'MarkerEdgeColor', 'r', 'MarkerFaceColor', get(hLineSeries, 'Color'));
set(hAxes(2), 'View', [90, 90]);
set(hAxes(2), 'XLim', [(0 + 0.5), (numRows + 0.5)]);
set(get(hAxes(2), 'Title'), 'String', ['Signal #001'], ...
    'FontSize', fontSizeTitle);

hAxes(3) = subplot(3, 3, [8, 9]);
hLineSeries   = plot(vY2);
set(hLineSeries, 'Marker', '*', 'MarkerSize', markerSizeMedium, ...
    'MarkerEdgeColor', 'r', 'MarkerFaceColor', get(hLineSeries, 'Color'));
set(hAxes(3), 'XLim', [(0 + 0.5), (numCols + 0.5)]);
set(get(hAxes(3), 'Title'), 'String', ['Signal #002'], ...
    'FontSize', fontSizeTitle);

hRunTimer   = tic();
[shortDist, mShortDist, mShortPath] = CalcShortestPathDyn(mDistMtx, numRows, numCols);
runTime     = toc(hRunTimer);

disp([' ']);
disp(['Dynamic Programming Method']);
disp(['Run Time - ', num2str(runTime), ' [Sec]']);
disp(['Shortest Path Value - ', num2str(shortDist), ' [Sec]']);

mShortPathImg = ExtractShortestPath(mShortDist, mShortPath);

figure('Position', [100, 100, 900, 900])

hAxes(1) = subplot(3, 3, [2, 3, 5, 6]);
hImageObj   = DisplayMatrixValues(mShortDist, hAxes(1));
% set(hAxes(1), 'DataAspectRatio', [1, 1, 1]);
set(get(hAxes(1), 'Title'), 'String', ['Paths Matrix'], ...
    'FontSize', fontSizeTitle);

hAxes(2) = subplot(3, 3, [1, 4]);
hLineSeries   = plot(vY1);
set(hLineSeries, 'Marker', '*', 'MarkerSize', markerSizeMedium, ...
    'MarkerEdgeColor', 'r', 'MarkerFaceColor', get(hLineSeries, 'Color'));
set(hAxes(2), 'View', [90, 90]);
set(hAxes(2), 'XLim', [(0 + 0.5), (numRows + 0.5)]);
set(get(hAxes(2), 'Title'), 'String', ['Signal #001'], ...
    'FontSize', fontSizeTitle);

hAxes(3) = subplot(3, 3, [8, 9]);
hLineSeries   = plot(vY2);
set(hLineSeries, 'Marker', '*', 'MarkerSize', markerSizeMedium, ...
    'MarkerEdgeColor', 'r', 'MarkerFaceColor', get(hLineSeries, 'Color'));
set(hAxes(3), 'XLim', [(0 + 0.5), (numCols + 0.5)]);
set(get(hAxes(3), 'Title'), 'String', ['Signal #002'], ...
    'FontSize', fontSizeTitle);


vShortPathIdx = find(mShortPathImg);
[vShortPathRowIdx, vShortPathColIdx] = ind2sub([numRows, numCols], vShortPathIdx);

figure('Position', [100, 100, 900, 900]);

hAxes(1) = subplot(3, 3, [2, 3, 5, 6]);
hImageObj   = DisplayMatrixValues(mDistMtx, hAxes(1));
% set(hAxes(1), 'DataAspectRatio', [1, 1, 1]);
set(get(hAxes(1), 'Title'), 'String', ['Distance Matrix'], ...
    'FontSize', fontSizeTitle);
set(hAxes(1), 'NextPlot', 'add');
hScatterSeris = scatter(vShortPathColIdx, vShortPathRowIdx, 80, ...
    'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'r', 'MarkerFaceAlpha', 0.5);

hAxes(2) = subplot(3, 3, [1, 4]);
hLineSeries   = plot(vY1);
set(hLineSeries, 'Marker', '*', 'MarkerSize', markerSizeMedium, ...
    'MarkerEdgeColor', 'r', 'MarkerFaceColor', get(hLineSeries, 'Color'));
set(hAxes(2), 'View', [90, 90]);
set(hAxes(2), 'XLim', [(0 + 0.5), (numRows + 0.5)]);
set(get(hAxes(2), 'Title'), 'String', ['Signal #001'], ...
    'FontSize', fontSizeTitle);

hAxes(3) = subplot(3, 3, [8, 9]);
hLineSeries   = plot(vY2);
set(hLineSeries, 'Marker', '*', 'MarkerSize', markerSizeMedium, ...
    'MarkerEdgeColor', 'r', 'MarkerFaceColor', get(hLineSeries, 'Color'));
set(hAxes(3), 'XLim', [(0 + 0.5), (numCols + 0.5)]);
set(get(hAxes(3), 'Title'), 'String', ['Signal #002'], ...
    'FontSize', fontSizeTitle);


%%

numSamples  = 60;
noiseStd    = 0.05;
shiftVal    = 5;

vX                  = [-1:0.05:1];
vGaussianDerivative = CalcGaussianGradient(3, 4);
vSigmoidSignal      = 1 ./ (1 + exp(-5 * vX));
vBaseSignal         = [zeros([1, 15]), vSigmoidSignal, ...
    ones([1, 20]), (1 + vGaussianDerivative), ones([1, 20])];
vBaseSignal         = 10 * vBaseSignal(:);
numSamplesBase      = size(vBaseSignal, 1);
vXBase              = [1:numSamplesBase].';

vY1Idx = sort(randperm(numSamplesBase, numSamples));
vY2Idx = sort(randperm(numSamplesBase, numSamples));

vY1 = vBaseSignal(vY1Idx);
vY2 = vBaseSignal(vY2Idx);

vY1 = vY1 + (noiseStd * randn([numSamples, 1]));
vY2 = vY2 + (noiseStd * randn([numSamples, 1]));

hFigure = figure('Position', [100, 100, 900, 600]);
hAxes   = axes();
hLineSeries = plot(vXBase, vBaseSignal, ...
    vXBase(vY1Idx), vBaseSignal(vY1Idx), ...
    vXBase(vY2Idx), vBaseSignal(vY2Idx));
set(hLineSeries(1), 'LineWidth', lineWidthNormal);
set(hLineSeries(2), 'LineStyle', 'none', 'Marker', 'o', ...
    'MarkerSize', markerSizeLarge);
set(hLineSeries(3), 'LineStyle', 'none', 'Marker', '*', ...
    'MarkerSize', markerSizeLarge);
set(get(hAxes, 'Title'), 'String', ['Signals'], ...
    'FontSize', fontSizeTitle);
set(get(hAxes, 'XLabel'), 'String', ['Index Number'], ...
    'FontSize', fontSizeAxis);
set(get(hAxes, 'YLabel'), 'String', ['Value'], ...
    'FontSize', fontSizeAxis);
hLegend = ClickableLegend({['Base Signal'], ...
    ['Signal 001'], ['Signal 002']});
set(hLegend, 'FontSize', fontSizeAxis);

mDistMtx = bsxfun(@minus, vY1, vY2.') .^ 2;

numRows = size(mDistMtx, 1);
numCols = size(mDistMtx, 2);


% hFigure     = figure();
% hAxes       = axes();
% DisplayMatrixValues(mDistMtx, hAxes);
% 
% set(hAxes, 'DataAspectRatio', [1, 1, 1])


figure('Position', [100, 100, 900, 900]);

hAxes(1) = subplot(3, 3, [2, 3, 5, 6]);
hImageObj   = imagesc(hAxes(1), mDistMtx);
% set(hAxes(1), 'DataAspectRatio', [1, 1, 1]);
set(get(hAxes(1), 'Title'), 'String', ['Distance Matrix'], ...
    'FontSize', fontSizeTitle);

hAxes(2) = subplot(3, 3, [1, 4]);
hLineSeries   = plot(vY1);
set(hLineSeries, 'Marker', '*', 'MarkerSize', markerSizeMedium, ...
    'MarkerEdgeColor', 'r', 'MarkerFaceColor', get(hLineSeries, 'Color'));
set(hAxes(2), 'View', [90, 90]);
set(hAxes(2), 'XLim', [(0 + 0.5), (numRows + 0.5)]);
set(get(hAxes(2), 'Title'), 'String', ['Signal #001'], ...
    'FontSize', fontSizeTitle);

hAxes(3) = subplot(3, 3, [8, 9]);
hLineSeries   = plot(vY2);
set(hLineSeries, 'Marker', '*', 'MarkerSize', markerSizeMedium, ...
    'MarkerEdgeColor', 'r', 'MarkerFaceColor', get(hLineSeries, 'Color'));
set(hAxes(3), 'XLim', [(0 + 0.5), (numCols + 0.5)]);
set(get(hAxes(3), 'Title'), 'String', ['Signal #002'], ...
    'FontSize', fontSizeTitle);

hRunTimer   = tic();
[shortDist, mShortDist, mShortPath] = CalcShortestPathDyn(mDistMtx, numRows, numCols);
runTime     = toc(hRunTimer);

disp([' ']);
disp(['Dynamic Programming Method']);
disp(['Run Time - ', num2str(runTime), ' [Sec]']);
disp(['Shortest Path Value - ', num2str(shortDist), ' [Sec]']);

mShortPathImg = ExtractShortestPath(mShortDist, mShortPath);

figure('Position', [100, 100, 900, 900])

hAxes(1) = subplot(3, 3, [2, 3, 5, 6]);
hImageObj   = imagesc(hAxes(1), mShortDist);
% set(hAxes(1), 'DataAspectRatio', [1, 1, 1]);
set(get(hAxes(1), 'Title'), 'String', ['Paths Matrix'], ...
    'FontSize', fontSizeTitle);

hAxes(2) = subplot(3, 3, [1, 4]);
hLineSeries   = plot(vY1);
set(hLineSeries, 'Marker', '*', 'MarkerSize', markerSizeMedium, ...
    'MarkerEdgeColor', 'r', 'MarkerFaceColor', get(hLineSeries, 'Color'));
set(hAxes(2), 'View', [90, 90]);
set(hAxes(2), 'XLim', [(0 + 0.5), (numRows + 0.5)]);
set(get(hAxes(2), 'Title'), 'String', ['Signal #001'], ...
    'FontSize', fontSizeTitle);

hAxes(3) = subplot(3, 3, [8, 9]);
hLineSeries   = plot(vY2);
set(hLineSeries, 'Marker', '*', 'MarkerSize', markerSizeMedium, ...
    'MarkerEdgeColor', 'r', 'MarkerFaceColor', get(hLineSeries, 'Color'));
set(hAxes(3), 'XLim', [(0 + 0.5), (numCols + 0.5)]);
set(get(hAxes(3), 'Title'), 'String', ['Signal #002'], ...
    'FontSize', fontSizeTitle);


vShortPathIdx = find(mShortPathImg);
[vShortPathRowIdx, vShortPathColIdx] = ind2sub([numRows, numCols], vShortPathIdx);

figure('Position', [100, 100, 900, 900]);

hAxes(1) = subplot(3, 3, [2, 3, 5, 6]);
hImageObj   = imagesc(hAxes(1), mDistMtx);
% set(hAxes(1), 'DataAspectRatio', [1, 1, 1]);
set(get(hAxes(1), 'Title'), 'String', ['Distance Matrix'], ...
    'FontSize', fontSizeTitle);
set(hAxes(1), 'NextPlot', 'add');
hScatterSeris = scatter(vShortPathColIdx, vShortPathRowIdx, 80, ...
    'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'r', 'MarkerFaceAlpha', 0.5);

hAxes(2) = subplot(3, 3, [1, 4]);
hLineSeries   = plot(vY1);
set(hLineSeries, 'Marker', '*', 'MarkerSize', markerSizeMedium, ...
    'MarkerEdgeColor', 'r', 'MarkerFaceColor', get(hLineSeries, 'Color'));
set(hAxes(2), 'View', [90, 90]);
set(hAxes(2), 'XLim', [(0 + 0.5), (numRows + 0.5)]);
set(get(hAxes(2), 'Title'), 'String', ['Signal #001'], ...
    'FontSize', fontSizeTitle);

hAxes(3) = subplot(3, 3, [8, 9]);
hLineSeries   = plot(vY2);
set(hLineSeries, 'Marker', '*', 'MarkerSize', markerSizeMedium, ...
    'MarkerEdgeColor', 'r', 'MarkerFaceColor', get(hLineSeries, 'Color'));
set(hAxes(3), 'XLim', [(0 + 0.5), (numCols + 0.5)]);
set(get(hAxes(3), 'Title'), 'String', ['Signal #002'], ...
    'FontSize', fontSizeTitle);


hFigure = figure('Position', [100, 100, 900, 600]);
hAxes   = axes();
hLineSeries = plot([vY1, (vY2 + shiftVal)]);
set(hLineSeries(1:2), 'LineWidth', lineWidthNormal);
DrawLinesBetweenCurves(vShortPathRowIdx, vY1, vShortPathColIdx, (vY2 + shiftVal), hAxes);
set(get(hAxes, 'Title'), 'String', ['Signals'], ...
    'FontSize', fontSizeTitle);
set(get(hAxes, 'XLabel'), 'String', ['Index Number'], ...
    'FontSize', fontSizeAxis);
set(get(hAxes, 'YLabel'), 'String', ['Value'], ...
    'FontSize', fontSizeAxis);
hLegend = ClickableLegend({['Signal 001'], ['Signal 002']});
set(hLegend, 'FontSize', fontSizeAxis);


%% Restore Defaults
% set(0, 'DefaultFigureWindowStyle', 'normal');
% set(0, 'DefaultAxesLooseInset', defaultLooseInset);

