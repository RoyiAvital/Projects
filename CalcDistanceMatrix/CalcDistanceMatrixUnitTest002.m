% ----------------------------------------------------------------------------------------------- %
% Calculate Dsitance Matrix Unit Test 009 - 'CalcDistanceMatrix()'
% Reference:
%   1. fd
% Remarks:
%   1.  Working on Float (Single).
%   2.  M
% TODO:
%   1.  A
%   Release Notes:
%   -   1.0.000     19/04/2018  Royi Avital
%       *   First release version.
% ----------------------------------------------------------------------------------------------- %

%% Setting Enviorment Parameters

run('InitScript.m');

COMPILING_MODE_DEBUG    = 1;
COMPILING_MODE_RELEASE  = 2;
COMPILING_MODE_GCC      = 3;
COMPILING_MODE_ICC      = 4;

LIB_NAME            = 'CalcDistanceMatrixDll';
H_FILE_NAME         = 'CalcDistanceMatrixDll';
LIB_PATH_DEBUG      = 'x64\Debug\';
LIB_PATH_RELEASE    = 'x64\Release\';
LIB_PATH_GCC        = 'x64\GCC\';
LIB_PATH_ICC        = 'x64\ICC\'; %<! Intel Compiler
H_FILE_PATH         = 'CalcDistanceMatrix\';


%% Settings

compilingMode   = COMPILING_MODE_GCC;


%% Loading Library

switch(compilingMode)
    case(COMPILING_MODE_DEBUG)
        libFullPath = [LIB_PATH_DEBUG, LIB_NAME, '.dll'];
    case(COMPILING_MODE_RELEASE)
        libFullPath = [LIB_PATH_RELEASE, LIB_NAME, '.dll'];
    case(COMPILING_MODE_GCC)
        libFullPath = [LIB_PATH_GCC, LIB_NAME, '.dll'];
    case(COMPILING_MODE_ICC)
        libFullPath = [LIB_PATH_ICC, LIB_NAME, '.dll'];
end

headerFullPath = [H_FILE_PATH, H_FILE_NAME, '.h'];

if(libisloaded(LIB_NAME) == FALSE)
    loadlibrary(libFullPath, headerFullPath);
    libfunctions(LIB_NAME, '-full');
end


%% Analysis 'ConvertFromUint8Fa()'

% funName         = 'CalcDistanceMatrixVanilla';
% funName         = 'CalcDistanceMatrixSse';
% funName         = 'CalcDistanceMatrixAvx';
funName         = 'CalcDistanceMatrixEigen';

vecDim      = 80;
numColsA    = 8000;
numColsB    = 6000;

mA = randn([vecDim, numColsA], 'single');
mB = randn([vecDim, numColsB], 'single');
mD = zeros([numColsA, numColsB], 'single');

hRunTime = tic();
% mDRef = CalcDistanceMatrix001(mA, mB); %<! Slower
mDRef = CalcDistanceMatrix002(mA, mB); %<! Faster
matlabRunTime = toc(hRunTime);

vecDim = int32(vecDim);
numColsA = int32(numColsA);
numColsB = int32(numColsB);

hRunTime = tic();
mD = calllib(LIB_NAME, funName, mD, mA, mB, vecDim, numColsA, numColsB);
dllRunTime = toc(hRunTime);

vE = mD(:) - mDRef(:);
maxAbsErr = max(abs(vE));
rmseErr = sqrt(mean(vE .* vE));

disp([' ']);
disp([funName, ' Unit Test']);
disp(['Max Abs Error - ', num2str(maxAbsErr)]);
disp(['RMSE Error - ', num2str(rmseErr)]);
disp(['MATLAB Run Time - ', num2str(matlabRunTime)]);
disp(['DLL Run Time - ', num2str(dllRunTime)]);
disp([' ']);


%% Restore Defaults

if(libisloaded(LIB_NAME) == TRUE)
    unloadlibrary(LIB_NAME);
end

if(libisloaded(LIB_NAME) == TRUE)
    disp('Error UnLoading ', LIB_NAME, ' Library');
end

% set(0, 'DefaultFigureWindowStyle', 'normal');
% set(0, 'DefaultAxesLooseInset', defaultLoosInset);

