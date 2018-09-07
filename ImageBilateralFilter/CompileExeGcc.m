% Compile DLL Using GCC
%
% References:
%   1.  https://linux.die.net/man/1/gcc.
%   2.  http://www.transmissionzero.co.uk/computing/building-dlls-with-mingw/.
% Remarks:
%   1.  Matrix is assumed to be n x n in size.
% TODO:
% 	1.  Add Eigen based implementation.
% Release Notes
% - 1.0.000     18/07/2017
%   *   First release.


%% General Parameters

run('InitScript.m');

figureIdx           = 0;
figureCounterSpec   = '%04d';

generateFigures = OFF;

GCC_FOLDER_PATH                 = 'D:\Applications\Programming\MinGW\';
GCC_BIN_FOLDER_PATH             = 'D:\Applications\Programming\MinGW\bin\';


%% Set System Enviorment (PATH Variable)

orgSystemPath = getenv('PATH');

if(contains(orgSystemPath, GCC_BIN_FOLDER_PATH) == FALSE)
    setenv('PATH', [orgSystemPath, GCC_BIN_FOLDER_PATH, ';']);
end
setenv('C_INCLUDE_PATH', [GCC_FOLDER_PATH, 'include\']);
setenv('CPLUS_INCLUDE_PATH', [GCC_FOLDER_PATH, 'include\']);


%% Set Files Names


%% Set Compiling Flags


%% Compile

% gccCommand01 = ['g++ -Ofast -fopenmp IpaLib/IpaLibMain.cpp IpaLib/IpaLibUnitTest.cpp IpaLib/BilateralFilterFa.cpp -o IpaLibExeGcc.exe'];
gccCommand01 = ['gcc -Ofast -fopenmp -mavx2 ImageBilateralFilter/ImageBilateralFilterMain.c ImageBilateralFilter/ImageBilateralFilterUnitTest.c ImageBilateralFilter/ImageBilateralFilter.c ImageBilateralFilter/ImageConvolution.c ImageBilateralFilter/ImageConvolutionGaussianKernel.c ImageBilateralFilter/ImageConvolutionSeparableKernel.c -o ImageBilateralFilter.exe'];

system(gccCommand01);

movefile('ImageBilateralFilter.exe', 'x64/GCC/');


%% Restore Defaults

% set(0, 'DefaultFigureWindowStyle', 'normal');
% set(0, 'DefaultAxesLooseInset', defaultLoosInset);

