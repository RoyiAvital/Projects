% ----------------------------------------------------------------------------------------------- %
% Calculate Dsitance Matrix Unit Test 002 - 'CalcDistanceMatrix()'
% Run Time Analysis
% Reference:
%   1. fd
% Remarks:
%   1.  Working on Float (Single).
%   2.  M
% TODO:
%   1.  A
%   Release Notes:
%   -   1.0.000     19/04/2018  Royi Avital
%       *   First release version.
% ----------------------------------------------------------------------------------------------- %

%% Setting Enviorment Parameters

run('InitScript.m');

SEC_TO_MILI_SEC_FCTR = 1e3;

COMPILING_MODE_DEBUG    = 1;
COMPILING_MODE_RELEASE  = 2;
COMPILING_MODE_MSVC     = COMPILING_MODE_RELEASE;
COMPILING_MODE_GCC      = 3;
COMPILING_MODE_ICC      = 4;

LIB_HANDLE_MODE_LOAD    = 1;
LIB_HANDLE_MODE_UNLOAD  = 2;

LIB_NAME            = 'CalcDistanceMatrixDll';

generateFigures = OFF;


%% Settings

vCompilingMode  = [COMPILING_MODE_GCC; COMPILING_MODE_ICC];
cCompilerString = {['GCC Compiler'], ['ICC Compiler']};
cFunName        = {['CalcDistanceMatrixVanilla'], ['CalcDistanceMatrixSse'], ['CalcDistanceMatrixAvx'], ['CalcDistanceMatrixEigen']};
mDataDim        = [5, 10, 10; 10, 100, 100; 20, 1000, 1000; 40, 2000, 2000; 80, 6000, 4000]; %<! Each Row - vecDim, numColsA, numColsB
mDataDim        = [20, 1000, 1000; 40, 2000, 2000; 80, 4000, 4000]; %<! Each Row - vecDim, numColsA, numColsB

numIterations = 50;


%% Setting Data

numCompilingModes   = length(vCompilingMode);
numFunctions        = length(cFunName);
numDim              = size(mDataDim, 1);

mRunTime = zeros([numIterations, numDim, numFunctions, numCompilingModes]);


%% Runtime Analysis

for compModIdx = 1:numCompilingModes
    
    compilingMode   = vCompilingMode(compModIdx);
    HandleDynamicLibrary(compilingMode, LIB_HANDLE_MODE_LOAD);
    
    for funIdx = 1:numFunctions
        
        funName         = cFunName{funIdx};
        
        for dimComb = 1:numDim
            vecDim      = mDataDim(dimComb, 1);
            numColsA    = mDataDim(dimComb, 2);
            numColsB    = mDataDim(dimComb, 3);
            
            mA = randn([vecDim, numColsA], 'single');
            mB = randn([vecDim, numColsB], 'single');
            mD = zeros([numColsA, numColsB], 'single');
            for ii = 1:numIterations
                
                vecDim = int32(vecDim);
                numColsA = int32(numColsA);
                numColsB = int32(numColsB);
                
                hRunTime = tic();
                mD = calllib(LIB_NAME, funName, mD, mA, mB, vecDim, numColsA, numColsB);
                runTime = toc(hRunTime);
                
                mRunTime(ii, dimComb, funIdx, compModIdx) = runTime;
                
            end
        end
    end
    
    HandleDynamicLibrary(compilingMode, LIB_HANDLE_MODE_UNLOAD);
    
end

mRunTime = mRunTime * SEC_TO_MILI_SEC_FCTR;


%% Data Analysis

cMeasureString = {['Mean Run Time'], ['Median Run Time'], ['Max Run Time'], ['Min Run Time']};

mMeanRunTime    = zeros([1, numDim, numFunctions, numCompilingModes]);
mMeaidanRunTime = zeros([1, numDim, numFunctions, numCompilingModes]);
mMaxRunTime     = zeros([1, numDim, numFunctions, numCompilingModes]);
mMinRunTime     = zeros([1, numDim, numFunctions, numCompilingModes]);

for ii = 1:numCompilingModes
    for jj = 1:numFunctions
        mMeanRunTime(1, :, jj, ii)     = mean(mRunTime(:, :, jj, ii), 1);
        mMeaidanRunTime(1, :, jj, ii)  = median(mRunTime(:, :, jj, ii), 1);
        mMaxRunTime(1, :, jj, ii)      = max(mRunTime(:, :, jj, ii), [], 1);
        mMinRunTime(1, :, jj, ii)      = min(mRunTime(:, :, jj, ii), [], 1);
    end
end

cRunTime = {mMeanRunTime, mMeaidanRunTime, mMaxRunTime, mMinRunTime};


%% Display Results

numMeasures = length(cMeasureString);

for ii = 1:numCompilingModes
    hFigure = figure('Position', figPosX2Large);
    for jj = 1:numMeasures
        hAxes = subplot(2, 2, jj);
        
        for kk = 1:numFunctions
            hLineObj = line(1:numDim, cRunTime{jj}(1, :, kk, ii));
            set(hLineObj, 'LineWidth', lineWidthNormal);
            set(hLineObj, 'Color', mColorOrder(kk, :));
        end
        set(get(hAxes, 'Title'), 'String', {[cMeasureString{jj}, ' - ', cCompilerString{ii}]}, ...
            'FontSize', fontSizeTitle);
        set(get(hAxes, 'XLabel'), 'String', {['Data Dimension Combination']}, ...
            'FontSize', fontSizeAxis);
        set(get(hAxes, 'YLabel'), 'String', {['Run Time [Mili Sec]']}, ...
            'FontSize', fontSizeAxis);
        set(hAxes, 'YLim', [0, 50]);
        hLegend = ClickableLegend(cFunName);

        if(generateFigures == ON)
            saveas(hFigure,['Figure', num2str(figureIdx, figureCounterSpec), '.png']);
        end
        
        
    end
end


%% Restore Defaults

% set(0, 'DefaultFigureWindowStyle', 'normal');
% set(0, 'DefaultAxesLooseInset', defaultLoosInset);

