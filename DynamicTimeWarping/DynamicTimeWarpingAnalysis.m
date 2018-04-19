% -------------------------------- Dynamic Time Warping Analysis -------------------------------- %

%% General Parameters and Initialization

clear();
close('all');

% set(0, 'DefaultFigureWindowStyle', 'docked');
defaultLooseInset = get(0, 'DefaultAxesLooseInset');
% set(0, 'DefaultAxesLooseInset', [0.05, 0.05, 0.05, 0.05]);

titleFontSize   = 14;
axisFotnSize    = 12;
stringFontSize  = 12;

thinLineWidth   = 2;
normalLineWidth = 3;
thickLineWidth  = 4;

smallMarkerSize     = 6;
mediumMarkerSize    = 8;
largeMarkerSize     = 10;

smallSizeData   = 36;
mediumSizeData  = 48;
bigSizeData     = 60;

randomNumberStream = RandStream('mlfg6331_64', 'NormalTransform', 'Ziggurat');
subStreamNumber = 57162;
subStreamNumber = 2108;
subStreamNumber = round(sum(clock()));
set(randomNumberStream, 'Substream', subStreamNumber);
RandStream.setGlobalStream(randomNumberStream);


%% Setting Constants

FALSE   = 0;
TRUE    = 1;

OFF = 0;
ON  = 1;


%% Setting Parameters

numSamples1 = 150;
numSamples2 = 125;

vRefSignal1 = 1 * sin(2 * pi * (40 / 1e3) * [9:109].');
vRefSignal2 = 1 * sin(2 * pi * (20 / 4e2) * [26:189].');


%% Create Data

% numRows = numSamples1;
% numCols = numSamples2;

numRows = size(vRefSignal1, 1);
numCols = size(vRefSignal2, 1);

mDistMtx = [bsxfun(@minus, vRefSignal1, vRefSignal2.') .^ 2];

hRunTimer   = tic();
[shortDist, mShortDist, mShortPath] = CalcShortestPathDyn(mDistMtx, numRows, numCols);
runTime     = toc(hRunTimer);

disp([' ']);
disp(['Dynamic Programming Method']);
disp(['Run Time - ', num2str(runTime), ' [Sec]']);
disp(['Shortest Path Value - ', num2str(shortDist), ' [Sec]']);

mShortPathImg = ExtractShortestPath(mShortDist, mShortPath);


hFigure         = figure();
set(hFigure, 'Units', 'pixels');
set(hFigure, 'Position', [75, 75, 750, 500]);

hAxes       = axes;
hImageObj   = imagesc(mShortPathImg);
set(hAxes, 'DataAspectRatio', [1, 1, 1]);
set(get(hAxes, 'Title'), 'String', ['Shortest Path'], ...
    'FontSize', titleFontSize);


%% Display Signals

figure('Position', [100, 100, 900, 900])

hAxes(1) = subplot(3, 3, [2, 3, 5, 6]);
hImageObj   = imagesc(mDistMtx);
% set(hAxes(1), 'DataAspectRatio', [1, 1, 1]);
set(get(hAxes(1), 'Title'), 'String', ['Distance Matrix'], ...
    'FontSize', titleFontSize);

hAxes(2) = subplot(3, 3, [1, 4]);
hLineSeries   = plot(vRefSignal1);
set(hAxes(2), 'View', [90, 90]);
set(hAxes(2), 'XLim', [1, numRows]);
set(get(hAxes(2), 'Title'), 'String', ['Signal #001'], ...
    'FontSize', titleFontSize);

hAxes(3) = subplot(3, 3, [8, 9]);
hLineSeries   = plot(vRefSignal2);
set(hAxes(3), 'XLim', [1, numCols]);
set(get(hAxes(3), 'Title'), 'String', ['Signal #002'], ...
    'FontSize', titleFontSize);





figure('Position', [100, 100, 900, 900])

hAxes(1) = subplot(3, 3, [2, 3, 5, 6]);
hImageObj   = imagesc(mShortPathImg);
% set(hAxes(1), 'DataAspectRatio', [1, 1, 1]);
% set(get(hAxes(1), 'Title'), 'String', ['Shortest Path'], ...
%     'FontSize', titleFontSize);

hAxes(2) = subplot(3, 3, [1, 4]);
hImageObj   = plot(vRefSignal1);
set(hAxes(2), 'View', [90, 90]);
set(hAxes(2), 'XLim', [1, numRows]);

hAxes(3) = subplot(3, 3, [8, 9]);
hImageObj   = plot(vRefSignal2);
set(hAxes(3), 'XLim', [1, numCols]);
% set(hAxes(1), 'DataAspectRatio', [1, 1, 1]);
% set(get(hAxes(1), 'Title'), 'String', ['Shortest Path'], ...
%     'FontSize', titleFontSize);


%% Restore Defaults
% set(0, 'DefaultFigureWindowStyle', 'normal');
% set(0, 'DefaultAxesLooseInset', defaultLooseInset);

