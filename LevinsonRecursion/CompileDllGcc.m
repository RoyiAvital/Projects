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

GCC_FOLDER_PATH                 = 'D:\Applications\Programming\MinGW64\';
GCC_BIN_FOLDER_PATH             = 'D:\Applications\Programming\MinGW64\bin\';


%% Set System Enviorment (PATH Variable)

orgSystemPath = getenv('PATH');

setenv('PATH', [orgSystemPath, GCC_BIN_FOLDER_PATH, ';']);
setenv('C_INCLUDE_PATH', [GCC_FOLDER_PATH, 'include\']);
setenv('CPLUS_INCLUDE_PATH', [GCC_FOLDER_PATH, 'include\']);


%% Set Compiling Flags


%% Compile

gccCommand01 = 'gcc -O3 -c -o LevinsonRecursionToeplitzMatrixGcc.o LevinsonRecursionToeplitzMatrix.cpp -D _USRDLL';
gccCommand02 = 'gcc -O3 -o LevinsonRecursionToeplitzMatrixGcc.dll -shared -s LevinsonRecursionToeplitzMatrixGcc.o';

% gccCommand01 = 'g++ -O3 -c -o LevinsonRecursionToeplitzMatrixGcc.o LevinsonRecursionToeplitzMatrix.cpp -D _USRDLL';
% gccCommand02 = 'g++ -O3 -o LevinsonRecursionToeplitzMatrixGcc.dll -shared -s LevinsonRecursionToeplitzMatrixGcc.o';

system(gccCommand01);
system(gccCommand02);


%% Restore Defaults

% set(0, 'DefaultFigureWindowStyle', 'normal');
% set(0, 'DefaultAxesLooseInset', defaultLoosInset);

