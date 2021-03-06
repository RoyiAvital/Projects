% ----------------------------------------------------------------------------------------------- %
% Image Convolution Unit Test 002
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

LIB_NAME            = 'ImageBilateralFilterDll';
H_FILE_NAME         = 'ImageBilateralFilterDll';
LIB_PATH_DEBUG      = 'x64/Debug/';
LIB_PATH_RELEASE    = 'x64/Release/';
LIB_PATH_GCC        = 'x64/GCC/';
H_FILE_PATH         = 'ImageBilateralFilter/';


%% Settings

funName         = 'ImageConvolutionSeparableKernel';
compilingMode   = COMPILING_MODE_GCC;

if(ispc())
    dyLibPostfix = DYNAMIC_LIB_POSTFIX_WIN;
elseif(isunix())
    dyLibPostfix = DYNAMIC_LIB_POSTFIX_LINUX;
elseif(ismac())
    dyLibPostfix = DYNAMIC_LIB_POSTFIX_MACOS;
else
    disp('Platform Is not Supported')
end


%% Loading Library


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


%% Analysis

numRows     = 8000; %<! Must be a factor of 4
numCols     = 8000; 

rowKernelLength = 29; %<! Must Be Odd
colKernelLength = 51; %<! Must Be Odd

mI = rand([numRows, numCols], 'single');

mO = zeros([numRows, numCols], 'single');
mTmp = zeros([numRows, numCols], 'single');

vRowKernel = randn([rowKernelLength, 1]);
% vRowKernel = (vRowKernel + rot90(vRowKernel, 2)) / 2;


vColKernel = randn([colKernelLength, 1]);
% vColKernel = (vColKernel + rot90(vColKernel, 2)) / 2;

vRowKernel = vRowKernel / sum(vRowKernel(:));
vColKernel = vColKernel / sum(vColKernel(:));


% mORef = single(imfilter(mI, double(vRowKernel * vColKernel.'), 'replicate', 'same', 'corr'));
% Remember MATLAB is Column Wise and C is Row Wise
% void ImageConvolutionSeparableKernel(float* mO, float* mI, float* mTmp, int numRows, int numCols, float* vRowKernel, int rowKernelLength, float* vColKernel, int colKernelLength)
% mO = calllib(LIB_NAME, funName, mO, mI, mTmp, numCols, numRows, vRowKernel, rowKernelLength, vColKernel, colKernelLength);

% Remember MATLAB is Column Wise and C is Row Wise
% void ImageConvolutionSeparableKernel(float* mO, float* mI, float* mTmp, int numRows, int numCols, float* vRowKernel, int rowKernelLength, float* vColKernel, int colKernelLength)
mORef = single(imfilter(mI, double(vColKernel * vRowKernel.'), 'replicate', 'same', 'corr'));
mO = calllib(LIB_NAME, funName, mO, mI, mTmp, numCols, numRows, vColKernel, colKernelLength, vRowKernel, rowKernelLength);


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

