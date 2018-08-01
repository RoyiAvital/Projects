% DTOA Localization CRLB Analysis
%
% References:
%   1.  aa
% Remarks:
%   1.  sa
% TODO:
% 	1.  ds
% Release Notes
% - 1.0.000     01/08/2018
%   *   First release.


%% General Parameters

run('InitScript.m');

figureIdx           = 0;
figureCounterSpec   = '%04d';

generateFigures = OFF;


%% Analysis Parameters


gridLength = 6000; %<! [Meter]
gridHeight = 2000;
gridNumPts = 25; %<! [Meter]

dtoaNoiseStd    = 10e-9; %<! [Sec]
waveSpeed       = 3e8; %<! [Meter / Sec]

mSensorLocation = [1500, 2500, 3500, 4500; 0, 0, 0, 0; 0, 0, 0, 0];

mSensorLocation = [1500, 2250, 3000, 3750, 4500; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0];


%% Load / Generate Data

numSensors = size(mSensorLocation, 2);

vGridXY     = linspace(0, gridLength, gridNumPts);
vGridZ      = linspace(10, gridHeight, gridNumPts);
mCrlbStd    = zeros(gridNumPts, gridNumPts, size(mSensorLocation, 1));

mG = zeros(size(mSensorLocation, 1), numSensors);

for jj = 1:gridNumPts

    for ii = 1:gridNumPts
        
        if(vGridXY(ii) < 50)
            continue;
        end
        
        for kk = 1:gridNumPts
            
            vP = [vGridXY(jj); vGridXY(ii); vGridZ(kk)];
            
            idxVal = 0;
            for mm = 1:numSensors - 1
                for nn = (mm + 1):numSensors
                    idxVal = idxVal + 1;
                    vGi = vP - mSensorLocation(:, mm);
                    vGi = vGi / norm(vGi);
                    vGj = vP - mSensorLocation(:, nn);
                    vGj = vGj / norm(vGj);
                    mG(:, idxVal) = vGi - vGj;
                end
            end
            
            mCrlbVar = ((dtoaNoiseStd * waveSpeed) ^ 2) * pinv(mG * mG.');
            
            %             if(any(mCrlbVar([1, 4, 9]) < 0))
            %                 keyboard;
            %             end
            
            mCrlbStd(ii, jj, 1) = sqrt(mCrlbVar(1));
            mCrlbStd(ii, jj, 2) = sqrt(mCrlbVar(4));
            mCrlbStd(ii, jj, 3) = sqrt(mCrlbVar(9));
            
        end
    end
end

mCrlbStd = abs(mCrlbStd);


%% Process Data

hFigure = figure('Position', figPosLarge);
hAxes = axes();
% set(hAxes, 'NextPlot', 'add');
hImageObj   = imagesc(vGridXY, vGridXY, mCrlbStd(:, :, 1));
set(hAxes, 'NextPlot', 'add');
hLineObj    = line(mSensorLocation(1, :), mSensorLocation(2, :));
set(hLineObj, 'Color', 'r', 'LineStyle', 'none', 'Marker', 'o', 'MarkerSize', markerSizeNormal, 'MarkerFaceColor', 'r');
set(hAxes, 'DataAspectRatio', [1, 1, 1]);
set(get(hAxes, 'Title'), 'String', {['CRLB of TDOA for X Coordinate Estimation [STD]']}, ...
    'FontSize', fontSizeTitle);
set(get(hAxes, 'XLabel'), 'String', 'X Coordinate [Meter]', ...
    'FontSize', fontSizeAxis);
set(get(hAxes, 'YLabel'), 'String', 'Y Coordinate [Meter]', ...
    'FontSize', fontSizeAxis);
hColorBar = colorbar(hAxes);
hLegend = ClickableLegend({['Sensors Location']});

hFigure = figure('Position', figPosLarge);
hAxes = axes();
% set(hAxes, 'NextPlot', 'add');
hImageObj   = imagesc(vGridXY, vGridXY, mCrlbStd(:, :, 2));
set(hAxes, 'NextPlot', 'add');
hLineObj    = line(mSensorLocation(1, :), mSensorLocation(2, :));
set(hLineObj, 'Color', 'r', 'LineStyle', 'none', 'Marker', 'o', 'MarkerSize', markerSizeNormal, 'MarkerFaceColor', 'r');
set(hAxes, 'DataAspectRatio', [1, 1, 1]);
set(get(hAxes, 'Title'), 'String', {['CRLB of TDOA for Y Coordinate Estimation [STD]']}, ...
    'FontSize', fontSizeTitle);
set(get(hAxes, 'XLabel'), 'String', 'X Coordinate [Meter]', ...
    'FontSize', fontSizeAxis);
set(get(hAxes, 'YLabel'), 'String', 'Y Coordinate [Meter]', ...
    'FontSize', fontSizeAxis);
hColorBar = colorbar(hAxes);
hLegend = ClickableLegend({['Sensors Location']});


%% Analyze Results


%% Display Results


%% Restore Defaults

% set(0, 'DefaultFigureWindowStyle', 'normal');
% set(0, 'DefaultAxesLooseInset', defaultLoosInset);
