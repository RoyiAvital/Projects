% Deep Learning Matched Filter - Generate Data
% Generating the data for Training and Validating procedures.
% References:
%   1.  aa
% Remarks:
%   1.  SNR is defined as Mean of Squred Signal Samples divided by Mean of
%       Squared Noise Samples.
% Known Issues:
%   1.  ds
% TODO:
% 	1.  ds
% Release Notes
% - 1.0.000     30/12/2017
%   *   First release.


%% General Parameters

run('InitScript.m');

figureIdx           = 0; %<! Continue from Question 1
figureCounterSpec   = '%04d';


%% Simulation Parameters

kernelRadius        = 45;
kernelStd           = 9;
kernelMeanEnergy    = 10;

numSamples      = 500;
numSignalsTrain = 1000; %<! Per SNR
numSignalsValid = 100; %<! Per SNR
vSnrLevel       = [10 ^ -2, 1, 10 ^ 2, 10 ^ 4];


%% Generate Data

vKernelGrid = [-kernelRadius:kernelRadius];
vKernelGrid = vKernelGrid(:);

vKernelSamples = exp(-(vKernelGrid .* vKernelGrid) ./ (2 * kernelStd * kernelStd));
vKernelSamples = sqrt(kernelMeanEnergy) * (vKernelSamples / sqrt(mean(vKernelSamples .^ 2)));

kernelLength = length(vKernelGrid);

numSnrLevels            = length(vSnrLevel);
numSignalsTrainTotal    = numSignalsTrain * numSnrLevels;
numSignalsValidTotal    = numSignalsValid * numSnrLevels;

mDataTrain = zeros([numSamples, numSignalsTrainTotal]);
mDataValid = zeros([numSamples, numSignalsValidTotal]);

signalIdxTrain = 0;
signalIdxValid = 0;

for ii = 1:numSnrLevels
    snrLevel = vSnrLevel(ii);
    noiseStd = sqrt(kernelMeanEnergy / snrLevel);
    for jj = 1:numSignalsTrain
        signalIdxTrain = signalIdxTrain + 1;
        vDataSignal = noiseStd * randn([numSamples, 1]);
        if(mod(signalIdxTrain, 2) == 1)
            kernelIdx = randi([1, (numSamples - kernelLength + 1)]);
            vDataSignal(kernelIdx:(kernelIdx + kernelLength - 1)) = vDataSignal(kernelIdx:(kernelIdx + kernelLength - 1)) + vKernelSamples;
        end
        mDataTrain(:, signalIdxTrain) = vDataSignal;
    end
    for jj = 1:numSignalsValid
        signalIdxValid = signalIdxValid + 1;
        vDataSignal = noiseStd * randn([numSamples, 1]);
        if(mod(signalIdxValid, 2) == 1)
            kernelIdx = randi([1, (numSamples - kernelLength + 1)]);
            vDataSignal(kernelIdx:(kernelIdx + kernelLength - 1)) = vDataSignal(kernelIdx:(kernelIdx + kernelLength - 1)) + vKernelSamples;
        end
        mDataValid(:, signalIdxValid) = vDataSignal;
    end
end

vDataLabelsTrain = mod(1:numSignalsTrainTotal, 2);
vDataLabelsValid = mod(1:numSignalsValidTotal, 2);


%% Save Data

save('Data\DeepLearningMatchedFilterData.m', 'mDataTrain', 'mDataValid', 'vDataLabelsTrain', 'vDataLabelsValid', 'kernelRadius', 'kernelStd', 'kernelMeanEnergy', 'numSamples', 'numSignalsTrain', 'numSignalsValid', 'vSnrLevel');


%% Display Results


% hAxes       = subplot_tight(2, 1, 2, [0.09, 0.09]);
% hLineSeries = plot(1:numIterations, 10 * log10(mSolErrNorm));
% set(hLineSeries, 'LineWidth', lineWidthNormal);
% set(get(hAxes, 'Title'), 'String', ['Solution Error Norm'], ...
%     'FontSize', fontSizeTitle);
% set(get(hAxes, 'XLabel'), 'String', 'Iteration Number', ...
%     'FontSize', fontSizeAxis);
% set(get(hAxes, 'YLabel'), 'String', 'Objective Function Value', ...
%     'FontSize', fontSizeAxis);
% set(hAxes, 'XLim', [1, numIterations]);
% hLegend = ClickableLegend(cLegendString);
% 
% if(generateFigures == ON)
%     saveas(hFigure,['Figure', num2str(figureIdx, figureCounterSpec), '.png']);
% end


%% Restore Defaults

% set(0, 'DefaultFigureWindowStyle', 'normal');
% set(0, 'DefaultAxesLooseInset', defaultLoosInset);

