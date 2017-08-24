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
% - 1.0.000     04/08/2017
%   *   First release.


%% General Parameters

run('../InitScript.m');

figureIdx           = 0;
figureCounterSpec   = '%04d';

generateFigures = OFF;


%% Set System Enviorment (PATH Variable)




%% Set Compiling Flags


%% Compile

gccCommand01 = 'gcc -O3 -fopenmp -fpic ../ImageConvolution/ImageConvolution.c ../ImageConvolution/ImageConvolutionGaussianKernel.c ../ImageConvolution/ImageConvolutionSeparableKernel.c -c -D _USRDLL';
gccCommand02 = 'gcc -O3 -fopenmp -o ImageConvolutionDll.so -shared -s ImageConvolution.o ImageConvolutionSeparableKernel.o ImageConvolutionGaussianKernel.o';

system(gccCommand01);
system(gccCommand02);


%% Restore Defaults

% set(0, 'DefaultFigureWindowStyle', 'normal');
% set(0, 'DefaultAxesLooseInset', defaultLoosInset);

