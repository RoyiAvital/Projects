function [ vClusterIdx, vMedoidIdx ] = ClusterLpStability( mD, paramMu, maxNumMedoids, debugMode )
% ----------------------------------------------------------------------------------------------- %
%[ vClusterIdx, vMedoidIdx ] = ClusterLpStability( mD, paramMu, debugMode )
% Clusters data according to its Distance Matrix by solving the relaxation
% of the Integer Linear Programming of Clustering with regularization on
% the number of clusters.
% Input:
%   - mD            -   Distance Matrix.
%                       The distance matrix of the data. The 'mD(ii, jj)'
%                       elements represent the distance between the 'ii'
%                       and the 'jj' samples. The matrix is assumed to be
%                       Symmetric.
%                       Structure: Matrix (numSamples x numSamples).
%                       Type: 'Single' / 'Double'.
%                       Range: [0, inf).
%   - paramMu       -   Parameter Mu.
%                       Scales the Median of distances in order to set the
%                       penalty / regularization on the number of clusters
%                       / medoids. The method use the penalty to measure
%                       the stability of the clusters. Hence it shouldn't
%                       be to sensitive to this.
%                       Structure: Scalar.
%                       Type: 'Single' / 'Double'.
%                       Range: (0, inf).
%   - maxNumMedoids -   Maximum Number of Medoids.
%                       Scales the maximum number of medoids. In case
%                       'paramMu' is low enough this will be the number of
%                       medoids (The K most dominant medoids).
%                       Structure: Scalar.
%                       Type: 'Single' / 'Double'.
%                       Range: (0, inf).
%   - debugMode     -   Debug Mode.
%                       Sets Debug Mode which prints to screen the
%                       information about the Margin between the Primal and
%                       Dual problems and finding new cluster.
%                       Structure: Scalar.
%                       Type: 'Single' / 'Double'.
%                       Range: {0, 1}.
% Output:
%   - vClusterIdx   -   Cluster Index.
%                       Describe per sample to which cluster it was
%                       assigned. If 'vClusterIdx(ii) = k' the the 'ii'
%                       samples is assigned to the 'k' cluster represented
%                       by the 'k' sample which is a medoid. 
%                       Number of clusters is 'length(unique(vClusterIdx))'
%                       or 'length(vMedoidIdx)'.
%                       Structure: Vector (numSamples x 1).
%                       Type: 'Single' / 'Double'.
%                       Range: {1, 2, ..., numSamples}.
%   - vMedoidIdx    -   Medoid Index.
%                       Vector of samples indices which are assigned as
%                       medoids. If 'vMedoidIdx(ii) = jj' it means that the
%                       'jj' sample is the medoid of the cluster it was
%                       assigned to.
%                       Structure: Vector (numClusters x 1).
%                       Type: 'Single' / 'Double'.
%                       Range: {1, 2, ..., numSamples}.
% References:
%   1.  Nikos Komodakis - Clustering via LP Based Stabilities (NIPS 2009).
%   2.  Lucas Fidon's code - https://github.com/LucasFidon/Clustering-via-LP-based-Stabilities. 
% Remarks:
%   1.  The algorithm extracts medoids by their "Dominance". Namely how
%       stable they are as a medoid.
%   2.  In practice 'paramMu' sets the number of clusters. Lower value
%       means more clusters while higher value means less clusters. It acts
%       as the "Cost" of new medoids.
%   3.  In case of low enough value of 'paramMu', namely the cost of adding
%       a centroid is low, all the samples will become medoids. In case
%       'maxNumMedoids' is lower than the number of sample and 'paramMu' is
%       low enough the algorithm will practically extract the most dominant
%       medoids.
% TODO:
%   1.  Optimize the data structure and allocations for better performance.
% Release Notes:
%   -   1.1.006     24/08/2020  Royi Avital
%       *   Preallocating the field 'sSolState.Q'.
%   -   1.1.005     21/08/2020  Royi Avital
%       *   Took a constant term calculation out of the loop.
%       *   Optimized performance for Set Operations.
%   -   1.1.004     04/06/2020  Royi Avital
%       *   Optimized some loops operations to make them O(n).
%   -   1.1.003     20/05/2020  Royi Avital
%       *   Optimized and vectroized some loop operations to make them
%           O(n) instead of O(n^2).
%   -   1.1.002     18/05/2020  Royi Avital
%       *   Optimized and vectroized some loop operations to make them
%           O(n) instead of O(n^2).
%   -   1.1.001     07/05/2020  Royi Avital
%       *   Optimized code to prevent some allocations.
%       *   Replacing 'inf' by large numbers.
%   -   1.1.000     25/04/2020  Royi Avital
%       *   Added option to set the maximum number of medoids by
%           'maxNumMedoids'.
%   -   1.0.001     25/04/2020  Royi Avital
%       *   Fixed crash when the value of paramMu was low enough to make
%           all samples centroids.
%   -   1.0.000     10/04/2020  Royi Avital
%       *   First release version.
% ----------------------------------------------------------------------------------------------- %

arguments
    mD {mustBeNumeric, mustBeReal, mustBeSquareMatrix}
    paramMu (1, 1) {mustBeNumeric, mustBeReal, mustBePositive} = 1
    maxNumMedoids (1, 1) {mustBeNumeric, mustBeReal, mustBeInteger, mustBePositive} = size(mD, 1);
    debugMode (1, 1) {mustBeMember(debugMode, [0, 1])} = 0
end

% Initialization
numSamples      = size(mD, 1);
medoidPenalty   = paramMu * median(mD(:));

maxNumMedoids = min(maxNumMedoids, numSamples);

sSolState               = struct();
sSolState.numSamples    = numSamples;
sSolState.mD            = mD + diag(medoidPenalty * ones(numSamples, 1));
sSolState.h             = sSolState.mD;
sSolState.numMedoids    = 0;
sSolState.Q             = zeros(1, maxNumMedoids);

% Running the algorithm
% The q parameter is the index of a new medoid
[sSolState, q] = SearchStablePoint(sSolState, debugMode);
while((CalcMargin(sSolState, q) >= 0)) %>! While adding new medoid makes sense
    sSolState.numMedoids = sSolState.numMedoids + 1;
    sSolState.Q(sSolState.numMedoids) = q;
    if(sSolState.numMedoids == maxNumMedoids)
        break;
    end
    if(debugMode)
        fprintf('\nadd point %d to the centroids set (total: %d centroids)\n', q, length(sSolState.Q));
    end
    sSolState       = ProjectMedoid(sSolState, q);
    [sSolState, q]  = SearchStablePoint(sSolState, debugMode);
end

% Extracting the medoids and cluster assignments
numMedoids  = sSolState.numMedoids;
vMedoidIdx  = sSolState.Q(1:numMedoids);
mD          = mD(vMedoidIdx, :);
[~, vIdx]   = min(mD, [], 1);

vClusterIdx = vMedoidIdx(vIdx); %<! Cluster labels are the medoid sample index


end

function [ sSolState ] = ProjectMedoid( sSolState, q )

numSamples  = sSolState.numSamples;
numMedoids  = sSolState.numMedoids;
mD          = sSolState.mD;
h           = sSolState.h;
Q           = sSolState.Q(1:numMedoids);

vP = SetDiffQ(Q, numSamples);
for pp = vP
    h(pp, pp) = h(pp, pp) + h(q, pp) - mD(q, pp);
    h(q, pp)  = mD(q, pp);
    h(pp, q)  = mD(pp, q);
end

h(q, q) = mD(q, q);
sSolState.h  = h;


end

function [ sSolState, q ] = SearchStablePoint( sSolState, debugMode )

valEps = 1e-5;

numSamples  = sSolState.numSamples;
numMedoids  = sSolState.numMedoids;
Q           = sSolState.Q(1:numMedoids);
Q_c         = SetDiffQ(Q, numSamples);
Nqc         = length(Q_c);
vMargin     = zeros(1, Nqc);
for ii = 1:Nqc
    vMargin(ii) = CalcMargin(sSolState, Q_c(ii));
end

dualObjVal  = CalcDualObjective(sSolState);
maxMargin   = max(vMargin);

if(debugMode)
    fprintf('\ndual sSolStateective = %.2f, max margin = %.2f\n', dualObjVal, maxMargin);
end

dualObjValPrev = 1e20; %<! Prevent using 'inf'
while((maxMargin < 0) && ((abs(dualObjVal - dualObjValPrev) / numSamples) > valEps))
    dualObjValPrev  = dualObjVal;
    sSolState       = Distribute(sSolState);
    dualObjVal      = CalcDualObjective(sSolState);
    vMargin(:)      = 0;
    for ii = 1:Nqc
        vMargin(ii) = CalcMargin(sSolState, Q_c(ii));
    end
    maxMargin  = max(vMargin);
    if(debugMode)
        fprintf('\ndual sSolStateective = %.2f, max margin = %.2f\n', dualObjVal, maxMargin);
    end
end

[~, maxIdx] = max(vMargin);
q           = Q_c(maxIdx);


end

function [ sSolState ] = Distribute( sSolState )

numSamples  = sSolState.numSamples;
numMedoids  = sSolState.numMedoids;
mD          = sSolState.mD;
h           = sSolState.h;
Q           = sSolState.Q(1:numMedoids);

Q_c         = SetDiffQ(Q, numSamples);
[h_min, nn] = min(h, [], 2);
h_thres     = h;
for pp = Q_c
    h_thres(pp, nn(pp)) = 1e30; %<! Prevent using Inf
end
h_hat = min(h_thres, [], 2);

Nqc     = length(Q_c);
vMargin = zeros(Nqc, 1);
for ii = 1:Nqc
    vMargin(ii) = CalcMargin(sSolState, Q_c(ii));
end

L_Q     = Q_c(IsMemberInt(nn(Q_c), Q));
vV      = ~(IsMemberInt(Q_c, L_Q))';
vIdx    = false(Nqc, 1);

for ii = 1:Nqc
    qq          = Q_c(ii);
    margin_q    = vMargin(ii);
    vIdx(:)     = vV & (h_min(Q_c) >= mD(Q_c, qq));
    
    card_V_q = sum(vIdx);
    if(~IsMemberScalar(qq, Q_c(vIdx)))
        card_V_q = card_V_q + 1;
    end        
    
    for pp = Q_c
        if((pp ~= qq) && (IsMemberScalar(pp, L_Q) || (h_min(pp) < mD(pp, qq))))
            h(pp, qq) = max(h_min(pp), mD(pp, qq));
        elseif(h(pp, qq) > h_min(pp))
            h(pp, qq) = h_min(pp) - margin_q / card_V_q;
        elseif(h(pp, qq) == h_min(pp))
            h(pp, qq) = h_hat(pp) - margin_q / card_V_q;
        end
    end
end

sSolState.h = h;


end

function [ dualObjVal ] = CalcDualObjective( sSolState )

dualObjVal = sum(min(sSolState.h, [], 2));


end

function [ valDelta ] = CalcMargin( sSolState, q )

valDelta    = 0;
numSamples  = sSolState.numSamples;
numMedoids  = sSolState.numMedoids;
mD          = sSolState.mD;
h           = sSolState.h;
Q           = sSolState.Q(1:numMedoids);

[~, nn] = min(h, [], 2);

[~, nn_notq]  = min(h(:, [1:(q - 1), (q + 1):numSamples]), [], 2);
vIdx          = nn_notq >= q;
nn_notq(vIdx) = nn_notq(vIdx) + 1;

vPP = SetDiffQ(Q, numSamples); %<! Not in Q

for pp = vPP
    if(h(pp, q) == h(pp, nn(pp)))
        valDelta = valDelta + h(pp, nn_notq(pp)) - h(pp, q);
    end
    if(pp ~= q)
        valDelta = valDelta - (h(pp, q) - max(h(pp, nn(pp)), mD(pp, q)));
    end
end

valDelta = valDelta - (h(q, q) - h(q, nn(q)));


end

function [ vD ] = SetDiffQ( vQ, numSamples )

numElements = length(vQ);

if(numElements == 0)
    vD = 1:numSamples;
    return;
end

vD  = zeros(1, numSamples - numElements);
vQs = sort(vQ, 'ascend');

idxQ = 1;
idxD = 1;
for ii = 1:numSamples
    if(vQs(idxQ) == ii)
        idxQ = idxQ + 1;
        idxQ = min(idxQ, numElements);
    else
        vD(idxD) = ii;
        idxD = idxD + 1;
    end
end


end

function [ isMember ] = IsMemberScalar( valIn, vA )

isMember = false();

if(isempty(vA))
    return;
end

isMember = any(valIn == vA);


end

function [ vF ] = IsMemberInt( vA, vB )

numElements = length(vA);

vBs = sort(vB, 'ascend');
vF  = false(size(vA));

for ii = 1:numElements
    vF(ii) = IsMemberValSorted(vBs, vA(ii));
end


end

function [ isMember ] = IsMemberValSorted( vA, valIn )
% Simple Binary Search to find if 'valIn' exists in 'vA'.


idxLeft     = 1;
idxRight    = length(vA);
isMember    = false();

while idxLeft <= idxRight
    idxMid = ceil((idxLeft + idxRight) / 2);
    
    if(vA(idxMid) == valIn)
        % If the index is needed
        % idxMem = idxMid;
        isMember = true();
        break;
    else
        if(vA(idxMid) > valIn)
            idxRight = idxMid - 1;
        else
            idxLeft = idxMid + 1;
        end
    end
end

% If the index is needed
% if(~isMember)
%     idxMem = -1;
% end


end


function mustBeSquareMatrix( mA )

% Test for Square Matrix

if(any(size(mA) < 1) || ~ismatrix(mA) || (size(mA, 1) ~= size(mA, 2)))
    error('Input must be a sqaure matrix')
end


end

