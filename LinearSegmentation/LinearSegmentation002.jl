# [![Visitors](https://hits.seeyoufarm.com/api/count/incr/badge.svg?url=https%3A%2F%2Fgithub.com%2FRoyiAvital%2FStackExchangeCodes&count_bg=%2379C83D&title_bg=%23555555&icon=&icon_color=%23E7E7E7&title=Visitors+%28Daily+%2F+Total%29&edge_flat=false)](https://github.com/RoyiAvital/Julia100Exercises)
# [![Royi Avital](https://i.imgur.com/ghq7NUE.png)](https://github.com/RoyiAvital/StackExchangeCodes)
# 
# # Linear Segmentation
# This document describe a [_Dynamic Programming_]() based method to solve the _linear segmentation_ problem.
# 
# > Notebook by:
# > - Royi Avital RoyiAvital@yahoo.com
#
# References:
#  1.   A
#
# Remarks:
#  1.   B
#
# To Do:
#  1.   C
# 
# ## Revision History
# 
# | Version | Date       | User        |Content / Changes                                                                         |
# |---------|------------|-------------|------------------------------------------------------------------------------------------|
# | 0.1.000 | 12/07/2023 | Royi Avital | First version                                                                            |

## Packages

# Internal
using DelimitedFiles;
using Printf;
using Statistics;           

# External
using CairoMakie;
using StableRNGs;
# using UnicodePlots;

## Constants & Configuration

oRng = StableRNG(123);
figureIdx = 0;

## Auxiliary Functions

include("LinearSegmentationFun.jl");


## Functions

# Loss Fun -> Minimize (Like Distance)
hLossFunMse(vY, vYY) = mean(abs2, vY - vYY); #<! vY Ground Truth, vYY - Estimation
# AffinityFun -> Maximize (Like Affinity)
hAffFunR2(vY, vYY) = 1 - (sum(abs2, vY .- vYY) / sum(abs2, mean(vY) .- vYY)); #<! vY Ground Truth, vYY - Estimation
hLossFunR2(vY, vYY) = -hAffFunR2(vY, vYY);

# Parameters

# Data
fileUrl   = raw"https://raw.githubusercontent.com/FixelAlgorithmsTeam/FixelCourses/refs/heads/master/DataSets/PieceWiseLinearData.csv";
decFactor = 3; #<! Decimation Factor

# Model
minSegLen = 10.0;
maxSegLen = 1000.0;
maxRmse   = 1.75;
maxDist   = 1e6;

# ## Load / Generate Data

mData, tuHeader = readdlm(download(fileUrl), ','; header = true);
vX = mData[1:decFactor:end, 1];
vY = mData[1:decFactor:end, 2];

# Display Data
figureIdx += 1;
figureFileName = @sprintf("%04d.png", figureIdx);

hF = Figure(size = (700, 700));
hA = Axis(hF, bbox = Rect2i((60, 60), (600, 600)), title = "Piece Wise Linear Model", xlabel = "x", ylabel = "y");
scatter!(hA, vX, vY;);
display(hF);
# save(figureFileName, hF);


# Build Cost Matrix / Distance Matrix
# mD = CalcDistMatReg(vX, vY, hLossFunMse; minLen = minSegLen, maxLen = maxSegLen, maxLoss = maxRmse * maxRmse, maxDist = maxDist);
mD = CalcDistMatReg(vX, vY, hLossFunMse);
mS, mP = SolveMinCostPartitionIntervals(mD, 5);
vP = ExtractPath(mS, mP); #<! Doesn't support NaN

vS = zeros(length(vX));

for ii in 1:length(vP)
    vS[vP[ii][1]:vP[ii][2]] .= ii;
end

# Display Data
figureIdx += 1;
figureFileName = @sprintf("%04d.png", figureIdx);

hF = Figure(size = (700, 700));
hA = Axis(hF, bbox = Rect2i((60, 60), (600, 600)), title = "Estimated Segments", xlabel = "x", ylabel = "y");
scatter!(hA, vX, vY; markersize = 20, color = vS);
display(hF);
# save(figureFileName, hF);
