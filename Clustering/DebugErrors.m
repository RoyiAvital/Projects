
clear('all'); clc();
load('Error005.mat');
% load('Error004.mat'); % Different first item

[vClusterIdx, vMedoidIdx]   = ClusterLpStability(mD, paramMu, maxNumMedoids);
vMedoidIdxM = ClusterLpStabilityMex(mD, paramMu, maxNumMedoids);

assert(isequal(sort(vMedoidIdxM(1)), sort(vMedoidIdx(1))));
assert(isequal(sort(vMedoidIdxM(:)), sort(vMedoidIdx(:))));
assert(isequal((vMedoidIdxM(:)), (vMedoidIdx(:))));
