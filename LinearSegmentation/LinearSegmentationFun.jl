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
using Statistics;           

# External

## Constants & Configuration

## Functions

function Conv1D( vA :: Vector{T}, vB :: Vector{T}; convMode :: String = "full" ) :: Vector{T} where {T <: Real}

    lenA = length(vA);
    lenB = length(vB);

    if (convMode == "full")
        startIdx    = 1;
        endIdx      = lenA + lenB - 1;
    elseif (convMode == "same")
        startIdx    = 1 + floor(Int, lenB / 2);
        endIdx      = startIdx + lenA - 1;
    elseif (convMode == "valid")
        startIdx    = lenB;
        endIdx      = lenA;
    end

    vO = zeros(T, lenA + lenB - 1);

    for idxB in 1:lenB
        @simd for idxA in 1:lenA
            @inbounds vO[idxA + idxB - 1] += vA[idxA] * vB[idxB];
        end
    end

    return vO[startIdx:endIdx];
end

function SolveMinCostPartitionIntervals( mD :: Matrix{T}, maxPartitions :: S  ) where {T, S <: Integer}

    numSamples = size(mD, 1);
    maxPartitions = min(maxPartitions, numSamples);
    mS = maximum(mD) .* ones(T, maxPartitions, numSamples); #<! Cost per segment
    mP = zeros(Int, maxPartitions, numSamples); #<! Path matrix

    mS[1, :] .= mD[1, :];
    mP[1, :] .= 1;
    for ii ∈ 2:maxPartitions
        for jj ∈ (ii + 1):numSamples
            minCost = 1e6;
            kkMin = 1;
            for kk ∈ 1:(ii - 1)
                currCost = mS[kk, ii - 1] + mD[ii, jj];
                if (currCost < minCost)
                    kkMin   = kk;
                    minCost = currCost;
                end
            mS[ii, jj] = minCost;
            mP[ii, jj] = kkMin;
            end
        end
    end

    return mS, mP;

end

function ExtractPath( mS :: Matrix{T}, mP :: Matrix{S} ) where {T, S <: Integer}

    numRows, numCols    = size(mS);
    vS                  = Vector{Vector{Int}}(undef, 0);

    startIdx = argmin(mS[:, numCols]);
    endIdx   = numCols;

    while ((startIdx > 0) && (endIdx > 0))
        prepend!(vS, [[startIdx, endIdx]]);
        colIdx      = mP[startIdx, endIdx];
        endIdx      = startIdx - 1;
        startIdx    = colIdx;
    end

    return vS;

end

function CalcDistMat( vX :: Vector{T} ) where {T}
    
    # numSamples  = length(vX);
    mD = ((vX .- vX') .^ 2);

    return mD;

end

function PolyFit( vX :: AbstractVector{T}, vY :: AbstractVector{T}, polyDeg :: Int ) where {T}

    numSamples = length(vX);
    mM = zeros(T, numSamples, polyDeg + 1);

    for ii = 1:(polyDeg + 1)
        mM[:, ii] = vX .^ (ii - 1);
    end

    return (mM \ vY);

end

function PolyVal( vX :: AbstractVector{T}, vP :: Vector{S}, polyDeg :: Int ) where {T, S}

    numSamples = length(vX);
    mM = zeros(T, numSamples, polyDeg + 1);

    for ii = 1:(polyDeg + 1)
        mM[:, ii] = vX .^ (ii - 1);
    end

    return mM * vP;

end


function CalcDistMatReg(vX :: Vector{T}, vY :: Vector{T}, hLossFun :: Function; minLen :: S = 0.0, maxLen :: S = 1.0, maxLoss :: S = 0.9, maxDist :: S = inf) where {T, S}
    # TODO: Use symmetric matrix

    numSamples = length(vX);
    mD = zeros(numSamples, numSamples);

    for ii in 1:numSamples, jj in ii:numSamples
        if (abs(vX[ii] - vX[jj]) < minLen)
            mD[ii, jj] = maxDist;
            mD[jj, ii] = maxDist;
            continue;
        end

        if (abs(vX[ii] - vX[jj]) > maxLen)
            mD[ii, jj] = maxDist;
            mD[jj, ii] = maxDist;
            continue;
        end

        @views vP = PolyFit(vX[ii:jj], vY[ii:jj], 1);
        @views vE = PolyVal(vX[ii:jj], vP, 1);
        @views estLoss = hLossFun(vY[ii:jj], vE);

        if (estLoss > maxLoss)
            mD[ii, jj] = maxDist;
            mD[jj, ii] = maxDist;
        else
            mD[ii, jj] = estLoss;
            mD[jj, ii] = estLoss;
        end
    end

    return mD;

end

