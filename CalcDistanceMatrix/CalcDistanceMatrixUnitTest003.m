% ----------------------------------------------------------------------------------------------- %
% Calc Distanc Matrix Unit Test 001
% Reference:
%   1. Checks only the Iterative K-Means.
% Remarks:
%   1.  Working within the [0, 1] range in RGB.
% TODO:
%   1.  A
%   Release Notes:
%   -   1.0.000     06/02/2017  Royi Avital
%       *   First release version.
% ----------------------------------------------------------------------------------------------- %

%% Setting Enviorment Parameters

run('InitScript.m');

% if(libisloaded('ClusterImageColorsAuxSseDll') == FALSE)
%     % loadlibrary('x64\Release\ClusterImageColorsAuxSseDll.dll', 'ClusterImageColorsAuxSse\ClusterImageColorsAuxSseDll.h');
%     loadlibrary('x64\Debug\ClusterImageColorsAuxSseDll.dll', 'ClusterImageColorsAuxSse\ClusterImageColorsAuxSseDll.h');
%     libfunctions('ClusterImageColorsAuxSseDll', '-full')
% end


%% Setting Parameters

% Data Parameters
numFeatures     = 10;
numSamples001   = 10000;
numSamples002   = 25;

numIterations = 50;


%% Generating Data

mDataSamples001 = randn([numFeatures, numSamples001], 'single');
mDataSamples002 = randn([numFeatures, numSamples002], 'single');


%% Running Algorithms

% Method 001
hRunTimer = tic();
for ii = 1:numIterations
    mDistMat001 = CalcDistanceMatrix001(mDataSamples001, mDataSamples002);
end
runTime = toc(hRunTimer);

disp(['Method 001 Run Time - ', num2str(runTime), ' [Sec]']);
disp([' ']);

% Method 002
hRunTimer = tic();
for ii = 1:numIterations
    mDistMat002 = CalcDistanceMatrix002(mDataSamples001, mDataSamples002);
end
runTime = toc(hRunTimer);

disp(['Method 002 Run Time - ', num2str(runTime), ' [Sec]']);
disp([' ']);


%% Displaying Results


%% Restoring Defaults

% if(libisloaded('ClusterImageColorsAuxSseDll') == TRUE)
%     unloadlibrary('ClusterImageColorsAuxSseDll');
% end
% 
% if(libisloaded('ClusterImageColorsAuxSseDll') == TRUE)
%     disp('Error UnLoading ClusterImageColorsAuxSseDll Library');
% end

% set(0, 'DefaultFigureWindowStyle', 'normal');
% set(0, 'DefaultAxesLooseInset', defaultLooseInset);

