function [ objVal ] = CalcObjVal( mD, vMedoidIdx )
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
%   1.  V
% TODO:
%   1.  C
% Release Notes:
%   -   1.0.000     04/09/2020  Or Yair
%       *   First release version.
% ----------------------------------------------------------------------------------------------- %

arguments
    mD {mustBeNumeric, mustBeReal, mustBeSquareMatrix}
    vMedoidIdx {mustBeNumeric, mustBeReal, mustBePositive}
end

[~, vIdx]   = min(mD(vMedoidIdx, :), [], 1);
vClusterIdx = vMedoidIdx(vIdx); %<! Cluster labels are the medoid sample index


objVal = 0;

for ii = 1:size(mD, 1)
    objVal = objVal + mD(ii, vClusterIdx(ii));
end


end


function mustBeSquareMatrix( mA )

% Test for Square Matrix

if(any(size(mA) < 1) || ~ismatrix(mA) || (size(mA, 1) ~= size(mA, 2)))
    error('Input must be a sqaure matrix')
end


end

