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

fileName = 'BilateralFilterFa.cpp';


%% Set Compiling Flags


%% Compile

gccCommand01 = ['gcc -Ofast -mavx2 -fopenmp -c -o ImageBilateralFilterDllGcc.o ImageBilateralFilter/ImageBilateralFilter.c -D _USRDLL'];
gccCommand02 = ['gcc -Ofast -mavx2 -fopenmp -c -o ImageConvolution.o ImageBilateralFilter/ImageConvolution.c -D _USRDLL'];
gccCommand03 = ['gcc -Ofast -mavx2 -fopenmp -c -o ImageConvolutionGaussianKernel.o ImageBilateralFilter/ImageConvolutionGaussianKernel.c -D _USRDLL'];
gccCommand04 = ['gcc -Ofast -mavx2 -fopenmp -c -o ImageConvolutionSeparableKernel.o ImageBilateralFilter/ImageConvolutionSeparableKernel.c -D _USRDLL'];
gccCommand05 = ['gcc -fopenmp  -o ImageBilateralFilterDll.dll -shared -s ImageBilateralFilterDllGcc.o ImageConvolution.o ImageConvolutionGaussianKernel.o ImageConvolutionSeparableKernel.o'];

% gccCommand01 = ['g++ -Ofast -fopenmp -c -o IpaLibDllGcc.o IpaLib/BilateralFilterFa.cpp -D _USRDLL'];
% gccCommand02 = ['g++ -Ofast -fopenmp -o IpaLibDllGcc.dll -shared -s IpaLibDllGcc.o'];

system(gccCommand01);
system(gccCommand02);
system(gccCommand03);
system(gccCommand04);
system(gccCommand05);

movefile('ImageBilateralFilterDll.dll', 'x64/GCC/');
delete('Image*.o');


%% Restore Defaults

% set(0, 'DefaultFigureWindowStyle', 'normal');
% set(0, 'DefaultAxesLooseInset', defaultLoosInset);

