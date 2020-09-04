% ----------------------------------------------------------------------------------------------- %
% 'ClusterLpStability()' Unit Test
% Reference:
%   1. fd
% Remarks:
%   1.  A
% TODO:
%   1.  A
%   Release Notes:
%   -   1.0.000     30/08/2020  Royi Avital
%       *   First release version.
% ----------------------------------------------------------------------------------------------- %

%% Setting Environment Parameters

subStreamNumberDefault = 192;
% subStreamNumberDefault = 542428563;
run('InitScript.m');

figureIdx           = 0;
figureCounterSpec   = '%04d';


%% Settings 

numTests        = 4;
numSamples      = 100;
maxNumMedoids   = 8;


%% Validation

maxDiff = 0;

for ii = 1:numTests
    
    paramMu = 5 * rand(1);
    
    if((ii == 1) || (ii == numTests))
        load('mD.mat');
    elseif(ii == 2)
        mA = [0, 0; 0.1, 0.1; -0.1, -0.1; 2, 2; 2.1, 2.1; 1.9, 1.9];
        mD = squareform(pdist(mA));
    else
        % mD = 5 * rand(numSamples, numSamples);
        % mD = mD + mD.';
        % SetDiag(mD, 0);
        
        mA = 5 * rand(numSamples, 2);
        mD = squareform(pdist(mA));
    end
    
    [vClusterIdx, vMedoidIdx]   = ClusterLpStability(mD, paramMu, maxNumMedoids);
    vMedoidIdxM = ClusterLpStabilityMex(mD, paramMu, maxNumMedoids);
    % [vClusterIdxM, vMedoidIdxM] = ClusterLpStabilityMex(mD, paramMu, maxNumMedoids);
    
    % assert(isequal(vClusterIdxM, vClusterIdx));
    % assert(isequal(vMedoidIdxM(:), vMedoidIdx(:)));
    % assert(isequal(vMedoidIdxM(1), vMedoidIdx(1)));
    % assert(isequal(sort(vMedoidIdxM(:)), sort(vMedoidIdx(:))));
    
    maxDiff = max(maxDiff, (CalcObjVal(mD, vMedoidIdx) - CalcObjVal(mD, vMedoidIdxM)));
    
    disp(['Finished Test #', num2str(ii, figureCounterSpec), ' Out of ', num2str(numTests), ' Tests']);
    
end

maxDiff


%% Run Time Analysis
% load('mD.mat');

hF = @() ClusterLpStabilityMex(mD, paramMu, numSamples);
hG = @() ClusterLpStability(mD, paramMu, numSamples);

runTimeRoyi = TimeItMin(hF, 1);
runTimeOrg  = TimeItMin(hG, 2);

disp(['MEX version is x', num2str(runTimeOrg / runTimeRoyi), ' faster']);


%% Restore Defaults

% set(0, 'DefaultFigureWindowStyle', 'normal');
% set(0, 'DefaultAxesLooseInset', defaultLoosInset);

