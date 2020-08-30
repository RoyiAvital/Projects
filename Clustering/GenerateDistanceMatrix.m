
numSamples = 50;
mMu = [0, 0; 1, 1; -1, 1; -1, -1; 1, -1];
paramSigma = 0.15;

numClusters = size(mMu, 1);

mA = zeros(numSamples * numClusters, 2);

startIdx    = 1;
endIdx      = numSamples;
for ii = 1:numClusters
    mA(startIdx:endIdx, :) = (paramSigma * randn(numSamples, 2)) + mMu(ii, :);
    
    startIdx    = startIdx + numSamples;
    endIdx      = endIdx + numSamples;
end

figure(); scatter(mA(:, 1), mA(:, 2));

mD = squareform(pdist(mA));
save('mD', 'mD');


% From http://cs.joensuu.fi/sipu/datasets/
% hFile = fopen('s1.txt', 'r');
% formatSpec = '%d';
% mA = reshape(fscanf(hFile, formatSpec), 2, 5000).';
% fclose(hFile);
% mA = mA ./ max(mA(:));
% % figure(); scatter(mA(:, 1), mA(:, 2));
% mD = squareform(pdist(mA));




