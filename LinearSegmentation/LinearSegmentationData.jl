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
using Statistics;           

# External
using StableRNGs;

## Constants & Configuration

oRng = StableRNG(123);

## Auxiliary Files
include("LinearSegmentationFun.jl");

## Functions

function GenLinearSineData( ampFiltSize, phaseFiltSize, vSeg )
    
    numSamples = vSeg[end] - 1;

    vAmp    = rand(oRng, numSamples);
    vAmp    = 0.2 * Conv1D(vAmp, ones(ampFiltSize) / ampFiltSize; convMode = "same");
    vPhase  = 0.2 * rand(oRng, numSamples);
    vPhase  = Conv1D(vPhase, ones(phaseFiltSize) / phaseFiltSize; convMode = "same");
    vPhase  = cumsum(vPhase);

    vX = LinRange(0, numSamples - 1, numSamples);

    vC = vAmp .* cos.(2 * pi * vPhase); #<! Harmonic Signal

    vL = zeros(numSamples); #<! Piece Wise Linear
    vL[vSeg[1]:(vSeg[2] - 1)] .= 0;
    vL[vSeg[2]:(vSeg[3] - 1)] .= 1;
    vL[vSeg[3]:(vSeg[4] - 1)] .= collect(LinRange(0.5, 1.0, vSeg[4] - vSeg[3]));
    vL[vSeg[4]:(vSeg[5] - 1)] .= collect(LinRange(1.0, 0.4, vSeg[5] - vSeg[4]));

    vY = vY = vC .+ vL;

    return vX, vY;

end


# Data for Introduction
vSeg = [001, 051, 101, 151, 201];
noiseStd = 0.01;

numSamples = vSeg[end] - 1;

vL = zeros(numSamples);
vL[vSeg[1]:(vSeg[2] - 1)] .= collect(LinRange(0.05, -0.05, vSeg[4] - vSeg[3]));
vL[vSeg[2]:(vSeg[3] - 1)] .= collect(LinRange(0.95, 1.05, vSeg[4] - vSeg[3]));
vL[vSeg[3]:(vSeg[4] - 1)] .= collect(LinRange(0.5, 1.0, vSeg[4] - vSeg[3]));
vL[vSeg[4]:(vSeg[5] - 1)] .= collect(LinRange(1.0, 0.4, vSeg[5] - vSeg[4]));

vY = vL + (noiseStd .* randn(numSamples));

vS = [ii for ii in 1:(length(vSeg) - 1) for jj in vSeg[ii]:(vSeg[ii + 1] - 1)];

vX = LinRange(0, numSamples - 1, numSamples);

writedlm("IntroData.csv", hcat(vX, vY, vS), ',');



# Data foe Example

vSeg = [1, 151, 301, 401, 501];

ampFiltSize   = 20;
phaseFiltSize = 50;


vX, vY = GenLinearSineData(ampFiltSize, phaseFiltSize, vSeg);
vS = [ii for ii in 1:(length(vSeg) - 1) for jj in vSeg[ii]:(vSeg[ii + 1] - 1)];

writedlm("ExampleData.csv", hcat(vX, vY, vS), ',');