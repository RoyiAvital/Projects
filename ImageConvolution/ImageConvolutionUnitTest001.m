% ----------------------------------------------------------------------------------------------- %
% Solve Linear Equations Eigen
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

LIB_NAME            = 'FaSseLibraryDll';
H_FILE_NAME         = 'FaSseLibDll';
LIB_PATH_DEBUG      = 'x64\Debug\';
LIB_PATH_RELEASE    = 'x64\Release\';
H_FILE_PATH         = 'FaSseLibrary\';


%% Settings

compilingMode = COMPILING_MODE_RELEASE;


%% Loading Library


switch(compilingMode)
    case(COMPILING_MODE_DEBUG)
        libFullPath = [LIB_PATH_DEBUG, LIB_NAME, '.dll'];
    case(COMPILING_MODE_RELEASE)
        libFullPath = [LIB_PATH_RELEASE, LIB_NAME, '.dll'];
end

headerFullPath = [H_FILE_PATH, H_FILE_NAME, '.h'];

if(libisloaded(LIB_NAME) == FALSE)
    loadlibrary(libFullPath, headerFullPath);
    libfunctions(LIB_NAME, '-full');
end


%% Analysis

numRows     = 800;
numCols     = 800;
gammaFctr   = single(2.1);

mA = rand([numRows, numCols], 'single');
mB = rand([numRows, numCols], 'single');

tic();
mB = calllib(LIB_NAME, 'ApplyGammaCurveFa' ,mB, mA, numCols, numRows, numCols, gammaFctr);
toc();

mBRef = mA .^ gammaFctr;

maxErr = max(abs(mB(:) - mBRef(:)));

disp(['Max Error - ', num2str(maxErr)]);


%% Restore Defaults

if(libisloaded(LIB_NAME) == TRUE)
    unloadlibrary(LIB_NAME);
end

if(libisloaded(LIB_NAME) == TRUE)
    disp('Error UnLoading ', LIB_NAME, ' Library');
end

% set(0, 'DefaultFigureWindowStyle', 'normal');
% set(0, 'DefaultAxesLooseInset', defaultLoosInset);

