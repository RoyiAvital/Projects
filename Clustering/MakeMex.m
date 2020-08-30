% ----------------------------------------------------------------------------------------------- %
% MakeMex - Generating x64 MEX File for ClusterLpStabilityMex
% Generates the MEX file of IncompleteCholeskyDecompositionMex by compiling
% the MEX wrapper.
% Reference:
%   1. See https://github.com/RoyiAvital/IncompleteCholeskyDecompositionThreshold.
% Remarks:
%   1.  Was verified on computer with MSVC 2019 Professional.
% TODO:
%   1.  A
%   Release Notes:
%   -   1.0.000     10/07/2020  Royi Avital
%       *   First release version.
% ----------------------------------------------------------------------------------------------- %


%% General Parameters

subStreamNumberDefault = 0; %<! Set to 0 for Random

run('InitScript.m');

SIMD_LEVEL_X64  = 1;
SIMD_LEVEL_SSE2 = 1;
SIMD_LEVEL_SSE4 = 2;
SIMD_LEVEL_AVX2 = 3;

OPTIMIZATION_LEVEL_O_FAST   = 1;
OPTIMIZATION_LEVEL_O_3      = 2;
OPTIMIZATION_LEVEL_O_2      = 3;

WIN_CRT_LIB_DYNAMIC = 1;
WIN_CRT_LIB_STATIC  = 2;

MEX_API_R2018A              = 1; %<! Large Arrays (64 Bit Indices), Interleaved Complex Arrays
MEX_API_R2017B              = 2; %<! Large Arrays (64 Bit Indices), Separate Complex Arrays
MEX_API_LARGE_ARRAY         = 3; %<! Large Arrays (64 Bit Indices), Separate Complex Arrays (Basically like R2017b)
MEX_API_COMPATIBLE_ARRAY    = 4; %<! 32 Bit Indices, Separate Complex Arrays


%% Parameters

enableFastMath      = ON;
simdLevel           = SIMD_LEVEL_AVX2;
optimizationLevel   = OPTIMIZATION_LEVEL_O_FAST;
enableOpenMP        = OFF;
enableAutoParallel  = OFF; %<! Windows only
winCrtLib           = WIN_CRT_LIB_STATIC; %<! Windows only

% MATLAB Configuration
mexApi = MEX_API_LARGE_ARRAY;

% Source File List
cFileList = {['ClusterLpStabilityMex.c']};


%% Setting Compilation Flags

switch(mexApi)
    case(MEX_API_R2018A)
        mexApiString = '-R2018a';
    case(MEX_API_R2017B)
        mexApiString = '-R2017b';
    case(MEX_API_LARGE_ARRAY)
        mexApiString = '-largeArrayDims';
    case(MEX_API_COMPATIBLE_ARRAY)
        mexApiString = '-compatibleArrayDims';
end

compStr = computer();
cCompFlags = [];

switch(compStr)
    case('PCWIN64')
        if(enableFastMath == ON)
            cCompFlags = [cCompFlags, ' /fp:fast'];
        end
        switch(simdLevel)
            case({SIMD_LEVEL_X64, SIMD_LEVEL_SSE2, SIMD_LEVEL_SSE4})
                % MSVC Defaults to SSE2, Doesn't support SSE4
                cCompFlags = cCompFlags;
            case(SIMD_LEVEL_AVX2)
                cCompFlags = [cCompFlags, ' /arch:AVX2'];
        end
        switch(optimizationLevel)
            case({OPTIMIZATION_LEVEL_O_2})
                cCompFlags = [cCompFlags, ' /O2'];
            case({OPTIMIZATION_LEVEL_O_3, OPTIMIZATION_LEVEL_O_FAST})
                cCompFlags = [cCompFlags, ' /O2 /Ob3 /Oi'];
        end
        if(winCrtLib == WIN_CRT_LIB_STATIC)
            cCompFlags = [cCompFlags, ' /MT'];
        end
        if(enableOpenMP == ON)
            cCompFlags = [cCompFlags, ' /openmp'];
        end
    case({'GLNXA64', 'MACI64'})
        if(enableFastMath == ON)
            cCompFlags = [cCompFlags, ' -ffast-math'];
        end
        switch(simdLevel)
            case({SIMD_LEVEL_X64, SIMD_LEVEL_SSE2})
                % MSVC Defaults to SSE2, Doesn't support SSE4
                cCompFlags = cCompFlags;
            case(SIMD_LEVEL_SSE4)
                cCompFlags = [cCompFlags, ' -msse4.2'];
            case(SIMD_LEVEL_AVX2)
                cCompFlags = [cCompFlags, ' -mavx2'];
        end
        switch(optimizationLevel)
            case({OPTIMIZATION_LEVEL_O_2})
                cCompFlags = [cCompFlags, ' -O2'];
            case({OPTIMIZATION_LEVEL_O_3, OPTIMIZATION_LEVEL_O_FAST})
                cCompFlags = [cCompFlags, ' -Ofast'];
        end
        if(enableOpenMP == ON)
            cCompFlags = [cCompFlags, ' -fopenmp'];
        end
end

cCompFlags = cCompFlags(2:end);

switch(compStr)
    case('PCWIN64')
        cCompFlags = sprintf('COMPFLAGS="$COMPFLAGS %s"', cCompFlags);
    case({'GLNXA64', 'MACI64'})
        cCompFlags = sprintf('CFLAGS=$CFLAGS %s', cCompFlags);
end


%% Compilation


for ii = 1:length(cFileList)
    mex('-v', mexApiString, cCompFlags, cFileList{ii});
end



%% Restore Defaults

% set(0, 'DefaultFigureWindowStyle', 'normal');
% set(0, 'DefaultAxesLooseInset', defaultLoosInset);

