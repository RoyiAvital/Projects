% Task Assignment Problem Analysis
% https://math.stackexchange.com/questions/2415617
% Analysis of the problem using MATLAB Simulation
% References:
%   1.  aa
% Remarks:
%   1.  sa
% TODO:
% 	1.  ds
% Release Notes
% - 1.0.000     03/09/2017
%   *   First release.


%% General Parameters

run('InitScript.m');
addpath('ImageToHTML/');

figureIdx           = 0; %<! Continue from Question 1
figureCounterSpec   = '%04d';

generateFigures = OFF;


%% Data Parameters

numTasks    = 5;
numWorkers  = 4;

hDecToBin = @(inputNum, numBits) (rem(floor(inputNum * pow2((1 - numBits):0)) ,2) == 1);


%% Generate Data

mG = randi([0, 100], [numTasks, numWorkers]); %<! Work Gain
mW = randi([0, 10], [numTasks, numWorkers]); %<! Work Resource

vC = randi([0, 100], [numWorkers, 1]); %<! Workers Capacity

NumWorkersMin = 3;
NumWorkersMax = 3;


%% Brute Force

numElements     = numTasks * numWorkers;
numCombinations = 2 ^ numElements;

sumGainOpt  = 0;
mX          = zeros([numTasks, numWorkers]);
mXOpt       = zeros([numTasks, numWorkers]);

hRunTime = tic();

for ii = 1:numCombinations
    validComb = TRUE;
    
    mX(:) = hDecToBin(ii - 1, numElements);
    
    % Capacity Validation
    validComb = validComb * all(sum(mX, 1) <= vC);
    if(validComb == FALSE)
        continue;
    end
    
    vT = any(mX, 2);
    
    % Minimum Number of Workers per Valid Task
    validComb = validComb * all(sum(mX, 2) >= NumWorkersMin);
    if(validComb == FALSE)
        continue;
    end
    
    % Maximum Number of Workers per Valid Task
    validComb = validComb * all(sum(mX, 2) <= NumWorkersMax);
    if(validComb == FALSE)
        continue;
    end
    
    sumGain = sum(sum(mG(vT, :) .* mX(vT, :)));
    
    if(sumGain > sumGainOpt)
        sumGainOpt  = sumGain;
        mXOpt       = mX;
    end
    
end

runTime = toc(hRunTime);

disp(['Total Run Time - ', num2str(runTime), ' [Sec]']);
disp(['Optimal Value Achieved - ', num2str(sumGainOpt)]);




figureIdx   = figureIdx + 1;
hFigure     = figure('Position', [100, 100, 1400, 900]);

for ii = 1:length(vEnergyThr)
    
    energyThr   = vEnergyThr(ii);
    mO          = CompressImageSvd(mI, energyThr, blockRadius);
    
    hAxes           = subplot_tight(2, 3, ii, [0.06, 0.06]);
    % hAxes           = subplot(1, 3, ii);
    hImageObject    = image(mO);
    set(get(hAxes, 'XLabel'), 'String', {['Energy Threhold = ', num2str(energyThr)]}, ...
        'FontSize', fontSizeAxis);
    set(hAxes, 'XTick', []);
    set(hAxes, 'XTickLabel', []);
    set(hAxes, 'YTick', []);
    set(hAxes, 'YTickLabel', []);
    set(hAxes, 'DataAspectRatio', [1, 1, 1]);
    % set(hAxes, 'LooseInset', [0.05, 0.05, 0.05, 0.05]);
    drawnow();
end

saveas(hFigure,['Figure', num2str(figureIdx, figureCounterSpec), '.png'], 'png');


%% Restore Defaults

% set(0, 'DefaultFigureWindowStyle', 'normal');
% set(0, 'DefaultAxesLooseInset', defaultLoosInset);

