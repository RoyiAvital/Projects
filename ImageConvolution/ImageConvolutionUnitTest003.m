% ----------------------------------------------------------------------------------------------- %
% Image Convolution Unit Test 003
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

LIB_NAME            = 'ImageConvolutionDll';
H_FILE_NAME         = 'ImageConvolutionDll';
LIB_PATH_DEBUG      = 'x64\Debug\';
LIB_PATH_RELEASE    = 'x64\Release\';
H_FILE_PATH         = 'ImageConvolutionDll\';


%% Settings

funName         = 'ImageConvolutionGaussianKernel';
compilingMode   = COMPILING_MODE_DEBUG;

numRows     = 800; %<! Must be a factor of 4
numCols     = 800;

gaussianStd         = 5;
stdToRadiusFactor   = 5;


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

kernelRadius = ceil(stdToRadiusFactor * gaussianStd);
kernelLength = (2 * kernelRadius) + 1;

vGaussianGrid   = [-kernelRadius:kernelRadius];
vGaussianKernel = exp(-(vGaussianGrid .^ 2) / (2 * gaussianStd * gaussianStd));
vGaussianKernel = vGaussianKernel / sum(vGaussianKernel);

mI = rand([numRows, numCols], 'single');

mO = zeros([numRows, numCols], 'single');
mTmp = zeros([numRows, numCols], 'single');

% Remember MATLAB is Column Wise and C is Row Wise
% void ImageConvolutionGaussianKernel(float* mO, float* mI, float* mTmp, int numRows, int numCols, float gaussianStd, int stdToRadiusFactor)
mORef = single(imfilter(mI, double(vGaussianKernel.' * vGaussianKernel), 'replicate', 'same', 'corr'));
mO = calllib(LIB_NAME, funName, mO, mI, mTmp, numCols, numRows, gaussianStd, stdToRadiusFactor);

maxErr = max(abs(mO(:) - mORef(:)));
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

