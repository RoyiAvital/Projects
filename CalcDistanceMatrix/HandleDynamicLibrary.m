function [  ] = HandleDynamicLibrary( compilingMode, libHandleMode )
%LOADDYNAMICLIBRARY Summary of this function goes here
%   Detailed explanation goes here

FALSE   = 0;
TRUE    = 1;

OFF     = 0;
ON      = 1;

COMPILING_MODE_DEBUG    = 1;
COMPILING_MODE_RELEASE  = 2;
COMPILING_MODE_MSVC     = COMPILING_MODE_RELEASE;
COMPILING_MODE_GCC      = 3;
COMPILING_MODE_ICC      = 4;

LIB_HANDLE_MODE_LOAD    = 1;
LIB_HANDLE_MODE_UNLOAD  = 2;

LIB_NAME            = 'CalcDistanceMatrixDll';
H_FILE_NAME         = 'CalcDistanceMatrixDll';
LIB_PATH_DEBUG      = 'x64\Debug\';
LIB_PATH_RELEASE    = 'x64\Release\';
LIB_PATH_GCC        = 'x64\GCC\';
LIB_PATH_ICC        = 'x64\ICC\'; %<! Intel Compiler
H_FILE_PATH         = 'CalcDistanceMatrix\';

switch(compilingMode)
    case(COMPILING_MODE_DEBUG)
        libFullPath = [LIB_PATH_DEBUG, LIB_NAME, '.dll'];
    case(COMPILING_MODE_MSVC)
        libFullPath = [LIB_PATH_RELEASE, LIB_NAME, '.dll'];
    case(COMPILING_MODE_GCC)
        libFullPath = [LIB_PATH_GCC, LIB_NAME, '.dll'];
    case(COMPILING_MODE_ICC)
        libFullPath = [LIB_PATH_ICC, LIB_NAME, '.dll'];
end

headerFullPath = [H_FILE_PATH, H_FILE_NAME, '.h'];

if(libHandleMode == LIB_HANDLE_MODE_LOAD)
    if(libisloaded(LIB_NAME) == FALSE)
        loadlibrary(libFullPath, headerFullPath);
        libfunctions(LIB_NAME, '-full');
    end
end

if(libHandleMode == LIB_HANDLE_MODE_UNLOAD)
    if(libisloaded(LIB_NAME) == TRUE)
        unloadlibrary(LIB_NAME);
    end
    
    if(libisloaded(LIB_NAME) == TRUE)
        disp('Error UnLoading ', LIB_NAME, ' Library');
    end
end


end

