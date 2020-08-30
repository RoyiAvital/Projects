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
%   2.  Pre sort 'vQ' and pre generate 'vP' so all other functions can
%       assume they are pre defined in the struct.
% Release Notes:
%   -   1.1.007     30/08/2020  Royi Avital
%       *   Optimized the use of 'CalcMargin()' in 'Distribute()' and 
%           'SearchStablePoint()' by removing 'vMargin'.
%       *   Pre defined 'vP' and sort 'vQ' as 'vQs'.
%   -   1.1.006     24/08/2020  Royi Avital
%       *   Preallocating the field 'sSolState.vQ'.
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
%   -   1.0.000     10/04/2020  Or Yair
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
sSolState.mH            = sSolState.mD;
sSolState.numMedoids    = 0;
sSolState.vQ            = zeros(1, maxNumMedoids);
sSolState.vQs           = zeros(1, maxNumMedoids); % Sorted
sSolState.numNotMedoids = numSamples;
sSolState.vP            = 1:numSamples;

% Running the algorithm
% The q parameter is the index of a new medoid
[sSolState, qq] = SearchStablePoint(sSolState, debugMode);
while(CalcMargin(sSolState, qq) >= 0) %>! While adding new medoid makes sense
    sSolState.numMedoids = sSolState.numMedoids + 1;
    sSolState.vQ(sSolState.numMedoids)  = qq;
    sSolState.vQs(sSolState.numMedoids) = qq;
    if(sSolState.numMedoids == maxNumMedoids)
        break;
    end
    if(debugMode)
        fprintf('\nadd point %d to the centroids set (total: %d centroids)\n', qq, length(sSolState.vQ));
    end
    sSolState.vQs(1:sSolState.numMedoids) = sort(sSolState.vQs(1:sSolState.numMedoids), 'ascend');
    sSolState.numNotMedoids = sSolState.numNotMedoids - 1;
    sSolState.vP(:) = UpdateP(sSolState.vP, sSolState.numNotMedoids, qq);
    sSolState       = ProjectMedoid(sSolState, qq);
    [sSolState, qq] = SearchStablePoint(sSolState, debugMode);
end

% Extracting the medoids and cluster assignments
numMedoids  = sSolState.numMedoids;
vMedoidIdx  = sSolState.vQ(1:numMedoids);
mD          = mD(vMedoidIdx, :);
[~, vIdx]   = min(mD, [], 1);

vClusterIdx = vMedoidIdx(vIdx); %<! Cluster labels are the medoid sample index


end

function [ sSolState ] = ProjectMedoid( sSolState, q )

mD          = sSolState.mD;
mH          = sSolState.mH;
numNotMedoids   = sSolState.numNotMedoids;
vP              = sSolState.vP(1:numNotMedoids);

for pp = vP
    mH(pp, pp) = mH(pp, pp) + mH(q, pp) - mD(q, pp);
    mH(q, pp)  = mD(q, pp);
    mH(pp, q)  = mD(pp, q);
end

mH(q, q) = mD(q, q);
sSolState.mH  = mH;


end

function [ sSolState, qq ] = SearchStablePoint( sSolState, debugMode )

valEps = 1e-5;

numSamples  = sSolState.numSamples;
numNotMedoids = sSolState.numNotMedoids;
vQc         = sSolState.vP(1:numNotMedoids);
Nqc         = numNotMedoids;

maxMargin   = CalcMargin(sSolState, vQc(1));
maxIdx      = 1;
for ii = 2:Nqc
    currMargin = CalcMargin(sSolState, vQc(ii));
    if(currMargin > maxMargin)
        maxMargin = currMargin;
        maxIdx = ii;
    end
end

dualObjVal          = CalcDualObjective(sSolState);

if(debugMode)
    fprintf('\ndual sSolStateective = %.2f, max margin = %.2f\n', dualObjVal, maxMargin);
end

dualObjValPrev = 1e20; %<! Prevent using 'inf'
while((maxMargin < 0) && ((abs(dualObjVal - dualObjValPrev) / numSamples) > valEps))
    dualObjValPrev  = dualObjVal;
    sSolState       = Distribute(sSolState);
    dualObjVal      = CalcDualObjective(sSolState);
    maxMargin   = CalcMargin(sSolState, vQc(1));
    maxIdx      = 1;
    for ii = 2:Nqc
        currMargin = CalcMargin(sSolState, vQc(ii));
        if(currMargin > maxMargin)
            maxMargin = currMargin;
            maxIdx = ii;
        end
    end
    if(debugMode)
        fprintf('\ndual sSolStateective = %.2f, max margin = %.2f\n', dualObjVal, maxMargin);
    end
end

qq          = vQc(maxIdx);


end

function [ sSolState ] = Distribute( sSolState )

numMedoids  = sSolState.numMedoids;
mD          = sSolState.mD;
mH          = sSolState.mH;
vQ          = sSolState.vQs(1:numMedoids);

numNotMedoids = sSolState.numNotMedoids;
vQc         = sSolState.vP(1:numNotMedoids);
[vNVal, vN] = min(mH, [], 2);
mHH         = mH;
for pp = vQc
    mHH(pp, vN(pp)) = 1e30; %<! Prevent using Inf
end
vH = min(mHH, [], 2);

Nqc = numNotMedoids;

vLq     = vQc(IsMemberInt(vN(vQc), vQ, true));
vV      = ~(IsMemberInt(vQc, vLq, true))';
vIdx    = false(Nqc, 1);

for ii = 1:Nqc
    qq          = vQc(ii);
    valMarginQ  = CalcMargin(sSolState, qq);
    vIdx(:)     = vV & (vNVal(vQc) >= mD(vQc, qq));
    
    valCardVq = sum(vIdx);
    if(~IsMemberScalar(qq, vQc(vIdx)))
        valCardVq = valCardVq + 1;
    end
    
    for pp = vQc
        if((pp ~= qq) && (IsMemberScalar(pp, vLq) || (vNVal(pp) < mD(pp, qq))))
            mH(pp, qq) = max(vNVal(pp), mD(pp, qq));
        elseif(mH(pp, qq) > vNVal(pp))
            mH(pp, qq) = vNVal(pp) - (valMarginQ / valCardVq);
        elseif(mH(pp, qq) == vNVal(pp))
            mH(pp, qq) = vH(pp) - (valMarginQ / valCardVq);
        end
    end
end

sSolState.mH = mH;


end

function [ dualObjVal ] = CalcDualObjective( sSolState )

dualObjVal = sum(min(sSolState.mH, [], 2));


end

function [ valMargin ] = CalcMargin( sSolState, qq )

valMargin   = 0;
numSamples  = sSolState.numSamples;
mD          = sSolState.mD;
mH          = sSolState.mH;

[~, vN] = min(mH, [], 2);

[~, vNq]    = min(mH(:, [1:(qq - 1), (qq + 1):numSamples]), [], 2);
vIdx        = vNq >= qq; % Compensation for removing the qq column
vNq(vIdx)   = vNq(vIdx) + 1; % Compensation for removing the qq column

numNotMedoids = sSolState.numNotMedoids;
vP = sSolState.vP(1:numNotMedoids);

for pp = vP
    if(mH(pp, qq) == mH(pp, vN(pp)))
        valMargin = valMargin + mH(pp, vNq(pp)) - mH(pp, qq);
    end
    if(pp ~= qq)
        valMargin = valMargin - (mH(pp, qq) - max(mH(pp, vN(pp)), mD(pp, qq)));
    end
end

valMargin = valMargin - (mH(qq, qq) - mH(qq, vN(qq)));


end

function [ vP ] = UpdateP( vP, numNotMedoids, qq )

shiftFlag = false;

for ii = 1:numNotMedoids
    shiftFlag = shiftFlag || (vP(ii) == qq);
    if(shiftFlag)
        vP(ii) = vP(ii + 1);
    end
end


end

function [ vD ] = SetDiffQ( vQ, numSamples )

numElements = length(vQ);

if(numElements == 0)
    vD = 1:numSamples;
    return;
end

vD = zeros(1, numSamples - numElements);
vQ = sort(vQ, 'ascend');

idxQ = 1;
idxD = 1;
for ii = 1:numSamples
    if(vQ(idxQ) == ii)
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

function [ vF ] = IsMemberInt( vA, vB, isBSorted )

numElements = length(vA);

if(~isBSorted)
    vB(:) = sort(vB, 'ascend');
end
vF  = false(size(vA));

for ii = 1:numElements
    vF(ii) = IsMemberValSorted(vA(ii), vB);
end


end

function [ isMember ] = IsMemberValSorted( valIn, vA )
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

