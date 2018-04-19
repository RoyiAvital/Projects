% ----------------------------------- Shortest Path Analysis ------------------------------------ %

%% General Parameters and Initialization

clear();
close('all');

% set(0, 'DefaultFigureWindowStyle', 'docked');
% defaultLooseInset = get(0, 'DefaultAxesLooseInset');
% set(0, 'DefaultAxesLooseInset', [0.05, 0.05, 0.05, 0.05]);

fontSizeTitle   = 14;
fontSizeAxis    = 12;
fontSizeString  = 12;

lineWidthThin   = 2;
lineWidthNormal = 3;
lineWidthThick  = 4;

markerSizeSmall     = 6;
markerSizeMedium    = 8;
markerSizeLarge     = 10;

dataSizeSmall   = 36;
dataSizeMedium  = 48;
dataSizeLarge   = 60;

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

numRows = 1000;
numCols = 1000;


%% Create Data

mDistMtx        = randi([2, 10], [numRows, numCols]);
% mDistMtx(:, numCols) = 1;

% for ii = 1:(numRows - 1)
%     mDistMtx(ii, ii) = 1;
% end
% 
% mDistMtx(numRows, (numCols - 1)) = 1;


%% Display Results

% hRunTimer   = tic();
% shortDist   = CalcShortestPathRec(mDistMtx, numRows, numCols);
% runTime     = toc(hRunTimer);
% 
% disp([' ']);
% disp(['Recursive Method']);
% disp(['Run Time - ', num2str(runTime), ' [Sec]']);
% disp(['Shortest Path Value - ', num2str(shortDist), ' [Sec]']);

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
set(get(hAxes(1), 'Title'), 'String', ['Shortest Path'], ...
    'FontSize', fontSizeTitle);



% set(hLineSeries(1), 'LineWidth', normalLineWidth);
% set(hLineSeries(2), 'LineStyle', 'none');
% set(hLineSeries(2), 'Marker', 'o');
% set(hLineSeries(2), 'LineWidth', normalLineWidth);
% set(get(hAxes(1), 'Title'), 'String', 'Unit Test - Local Peak Detector', ...
%     'FontSize', titleFontSize);
% set(get(hAxes(1), 'XLabel'), 'String', 'Sample Index', ...
%     'FontSize', axisFotnSize);
% set(get(hAxes(1), 'YLabel'), 'String', 'Sample Value', ...
%     'FontSize', axisFotnSize);
% hLegend = legend(hAxes, {['Input Signal'], ['Local Peaks']});
% set(hLegend, 'FontSize', axisFotnSize);


%% Restore Defaults
% set(0, 'DefaultFigureWindowStyle', 'normal');
% set(0, 'DefaultAxesLooseInset', defaultLooseInset);

