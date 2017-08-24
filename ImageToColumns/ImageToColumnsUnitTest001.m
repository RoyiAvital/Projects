% ----------------------------------------------------------------------------------------------- %
% Image Convolution Unit Test 001
% Reference:
%   1. fd
% Remarks:
%   1.  Working on Float (Single).
% TODO:
%   1.  A
%   Release Notes:
%   -   1.0.000     02/06/2017  Royi Avital
%       *   First release version.
% ----------------------------------------------------------------------------------------------- %

%% Setting Enviorment Parameters

run('InitScript.m');

COMPILING_MODE_DEBUG    = 1;
COMPILING_MODE_RELEASE  = 2;
COMPILING_MODE_GCC      = 3;

DYNAMIC_LIB_POSTFIX_WIN     = '.dll';
DYNAMIC_LIB_POSTFIX_LINUX   = '.so';
DYNAMIC_LIB_POSTFIX_MACOS   = '.dylib';

LIB_NAME            = 'ImageToColumnsDll';
H_FILE_NAME         = 'ImageToColumnsDll';
LIB_PATH_DEBUG      = 'x64/Debug/';
LIB_PATH_RELEASE    = 'x64/Release/';
LIB_PATH_GCC        = 'GCC/';
H_FILE_PATH         = 'ImageToColumnsDll/';


%% Settings

funName         = 'ImageToColumns';
compilingMode   = COMPILING_MODE_RELEASE;

numRows     = 1000;
numCols     = 1000;
blockRadius = 5;

numIter = 5;


%% Loading Library

if(ispc)
    dyLibPostfix = DYNAMIC_LIB_POSTFIX_WIN;
elseif(isunix)
    dyLibPostfix = DYNAMIC_LIB_POSTFIX_LINUX;
elseif(ismac)
    dyLibPostfix = DYNAMIC_LIB_POSTFIX_MACOS;
else
    disp('Platform Is not Supported')
end

switch(compilingMode)
    case(COMPILING_MODE_DEBUG)
        libFullPath = [LIB_PATH_DEBUG, LIB_NAME, dyLibPostfix];
    case(COMPILING_MODE_RELEASE)
        libFullPath = [LIB_PATH_RELEASE, LIB_NAME, dyLibPostfix];
    case(COMPILING_MODE_GCC)
        libFullPath = [LIB_PATH_GCC, LIB_NAME, dyLibPostfix];
end

headerFullPath = [H_FILE_PATH, H_FILE_NAME, '.h'];

if(libisloaded(LIB_NAME) == FALSE)
    loadlibrary(libFullPath, headerFullPath);
    libfunctions(LIB_NAME, '-full');
end


%% Generating Data

mI                  = rand([numRows, numCols], 'single');
% mI                  = rand([numRows, numCols], 'double');
blockSize           = (2 * blockRadius) + 1;
numPixelsToProcess  = (numRows - blockSize + 1) * (numCols - blockSize + 1);


%% Analysis

hRunTime = tic();

for ii = 1:numIter
    mORef = ImageToColumnsSliding(mI, [blockSize, blockSize]);
    % mORef = im2col(mI, [blockSize, blockSize], 'sliding');
end

matlabRunTime = toc(hRunTime) / numIter;

hRunTime = tic();
for ii = 1:numIter
    mO = zeros([(blockSize * blockSize), numPixelsToProcess], 'single');
    % mO = zeros([(blockSize * blockSize), numPixelsToProcess], 'double');
    mO = calllib(LIB_NAME, funName, mO, mI, numCols, numRows, blockRadius);
end
dllRunTime = toc(hRunTime) / numIter;

mE      = mO - mORef;
maxErr  = max(abs(mE(:)));
disp(['Max Error        - ', num2str(maxErr)]);
disp(['MATLAB Run Time  - ', num2str(matlabRunTime), ' [Sec]']);
disp(['DLL Run Time     - ', num2str(dllRunTime), ' [Sec]']);


%% Restore Defaults

if(libisloaded(LIB_NAME) == TRUE)
    unloadlibrary(LIB_NAME);
end

if(libisloaded(LIB_NAME) == TRUE)
    disp('Error UnLoading ', LIB_NAME, ' Library');
end

% set(0, 'DefaultFigureWindowStyle', 'normal');
% set(0, 'DefaultAxesLooseInset', defaultLoosInset);

