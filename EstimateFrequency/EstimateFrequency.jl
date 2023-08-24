# Estimate Sine Signal Frequency
# https://github.com/RoyiAvital/Projects/tree/master/EstimateFrequency
# Estimating the Frequency of a Sinusoidal Signal in White Noise
# References:
#   1.  
# Remarks:
#   1.  Use in Julia as following:
#       -   Move to folder using `cd(raw"<PathToFolder>");`.
#       -   Activate the environment using `] activate .`.
#       -   Instantiate the environment using `] instantiate`.
#   2.  fd
# TODO:
# 	1.  C
# Release Notes Royi Avital RoyiAvital@yahoo.com
# - 1.0.000     21/07/2023  Royi Avital
#   *   First release.

## Packages

# Internal
using LinearAlgebra;
using Printf;
using Statistics;
# External
using Optim;
using PlotlyJS;
using FFTW;


## Constants & Configuration

## External
# juliaInitPath = joinpath(".", "..", "..", "JuliaCode", "JuliaInit.jl")
# include(juliaInitPath)

## General Parameters

figureIdx       = 0;
exportFigures   = true;

## Functions

function EstFrequencyMLGrid( vY :: Vector{T}, samplingFreq :: T, vN :: S, mB :: Matrix{T}; numGridPts :: R = 1024 ) where {T <: AbstractFloat, S <: UnitRange, R <: Integer}
    # Estimates the frequency using the ML Estimator.
    # Does not assume knowledge of the amplitude and phase.
    # Mutate vP.

    numSamples = length(vY);

    vF = LinRange(0, 0.5, numGridPts);
    vF = LinRange(0.23, 0.27, 1000);

    lsValMin = T(Inf);
    estFreq  = 0;

    for (ii, ff) in enumerate(vF[2:(end - 1)]) #<! 0 < ff < 0.5

        lsVal = CalcObjFun(ff, vY, vN, mB);

        if (lsVal < lsValMin)
            lsValMin = lsVal;
            estFreq = ff;
        end
    end

    return estFreq * samplingFreq;

end

function EstFrequencyMLOpt( vY :: Vector{T}, samplingFreq :: T, vN :: S, mB :: Matrix{T} ) where {T <: AbstractFloat, S <: UnitRange}
    # Estimates the frequency using the ML Estimator.
    # Does not assume knowledge of the amplitude and phase.
    # Mutate vP.

    numSamples = length(vY);

    vYY = rfft(vY);
    ~, idxK = findmax(abs2, @view vYY[2:(end - 1)]);
    idxK += 1;

    # Frequency boundary
    # Index should start from 0
    fLow = max((idxK - 4) * (1 / numSamples), 0.025);
    fUp  = min((idxK + 4) * (1 / numSamples), 0.475);

    # It seems it uses a variant of the Brent algorithm which interpolates the derivative.
    # The bisection part is interpolated by 3 points instead of using the derivative.
    oSol = optimize(f -> CalcObjFun(f, vY, vN, mB), fLow, fUp, Brent()); #<! Faster and better
    # oSol = optimize(f -> CalcObjFun(f, vY, vN, mB), fLow, fUp, GoldenSection());

    estFreq = oSol.minimizer;    

    return estFreq * samplingFreq;

end

function EstFrequencyMLGridOpt( vY :: Vector{T}, samplingFreq :: T, vN :: S, mB :: Matrix{T}; numGridPts :: R = 1024 ) where {T <: AbstractFloat, S <: UnitRange, R <: Integer}
    # Estimates the frequency using the ML Estimator.
    # Does not assume knowledge of the amplitude and phase.
    # Mutate vP.

    numSamples = length(vY);

    vF = LinRange(0.025, 0.475, numGridPts);
    vV = 1e6 * ones(T, numGridPts); #<! TODO: Use the buffer

    # lsValMin = T(Inf);
    estFreq  = 0;

    for (ii, ff) in enumerate(vF) #<! 0 <= ff <= 0.5
        if ((ii == 1) || (ii == numGridPts))
            continue;
        end

        vV[ii] = CalcObjFun(ff, vY, vN, mB);
    end

    idxK    = argmin(vV);
    oSol    = optimize(f -> CalcObjFun(f, vY, vN, mB), vF[idxK - 1], vF[idxK + 1], Brent()); #<! Faster and better
    estFreq = oSol.minimizer;  

    return estFreq * samplingFreq;

end

# The 1st formulation
# function CalcObjFun( valF :: T, vY :: Vector{T}, vN :: S, mB :: Matrix{T} ) where {T <: AbstractFloat, S <: UnitRange}

#     numSamples = length(vY);

#     mX      = @view mB[1:numSamples, :];
#     mXX     = @view mB[(numSamples + 1):(numSamples + 2), :];
#     mXXInv  = @view mB[(numSamples + 3):(numSamples + 4), :];
#     mT211   = @view mB[(numSamples + 5):(numSamples + 6), 1];
#     mT212   = @view mB[(numSamples + 5):(numSamples + 6), 2];
#     mTN11   = @view mB[(numSamples + 7):(numSamples + 7 + numSamples - 1), 1];

#     mX[:, 1] .= sin.(2π .* valF .* vN);
#     mX[:, 2] .= cos.(2π .* valF .* vN);
#     # The loop adds allocations!
#     # for nn in 0:(numSamples - 1)
#     #     ii = nn + 1;
#     #     mX[ii, 1] = sin(2π * valF * nn);
#     #     mX[ii, 2] = cos(2π * valF * nn);
#     # end
#     mXX[1] = dot(mX[:, 1], mX[:, 1]);
#     mXX[2] = dot(mX[:, 1], mX[:, 2]);
#     mXX[3] = mXX[2];
#     mXX[4] = dot(mX[:, 2], mX[:, 2]);
#     detXXInv = 1 / ((mXX[1] * mXX[4]) - (mXX[2] * mXX[2]));
#     mXXInv[1] = detXXInv * mXX[4];
#     mXXInv[2] = -detXXInv * mXX[2];
#     mXXInv[3] = mXXInv[2];
#     mXXInv[4] = detXXInv * mXX[1];

#     mul!(mT211, mX', vY);
#     mul!(mT212, mXXInv, mT211);
#     mul!(mTN11, mX, mT212);
#     mTN11 .-= vY;

#     lsVal = sum(abs2, mTN11);

#     return lsVal;

# end

# The 2nd formulation (Shorter)
function CalcObjFun( valF :: T, vY :: Vector{T}, vN :: S, mB :: Matrix{T} ) where {T <: AbstractFloat, S <: UnitRange}

    numSamples = length(vY);

    mX      = @view mB[1:numSamples, :];
    mXX     = @view mB[(numSamples + 1):(numSamples + 2), :];
    mL      = @view mB[(numSamples + 3):(numSamples + 4), :];
    mT121   = @view mB[(numSamples + 5):(numSamples + 5), 1:2]; #<! 2:2 -> Makes is a 1x2 matrix
    mTN21   = @view mB[(numSamples + 6):(numSamples + 6 + numSamples - 1), :];

    mX[:, 1] .= sin.(2π .* valF .* vN);
    mX[:, 2] .= cos.(2π .* valF .* vN);

    mXX[1] = dot(mX[:, 1], mX[:, 1]);
    mXX[2] = dot(mX[:, 1], mX[:, 2]);
    mXX[3] = mXX[2];
    mXX[4] = dot(mX[:, 2], mX[:, 2]);
    
    CholeskyDecompositionInv2x2!(mL, mXX);

    mul!(mTN21, mX, mL); #<! numSamples x 2
    mul!(mT121, vY', mTN21); #<! 1 x 2

    lsVal = dot(mT121, mT121);

    return -lsVal;

end

function CholeskyDecomposition2x2!( mLR :: AbstractMatrix{T}, mA :: AbstractMatrix{T}; trType :: String = "TrLeft" ) where {T <: AbstractFloat}

    mLR[1] = sqrt(mA[1]);
    offDiagVal = mA[2] / mLR[1];
    if (trType == "TrLeft")
        mLR[2] = offDiagVal;
        mLR[3] = zero(T);
    elseif (trType == "TrRight")
        mLR[2] = zero(T);
        mLR[3] = offDiagVal;
    end
    mLR[4] = sqrt(mA[4] - (offDiagVal * offDiagVal));

    return mLR;

end

function CholeskyDecompositionInv2x2!( mLR :: AbstractMatrix{T}, mA :: AbstractMatrix{T}; trType :: String = "TrLeft" ) where {T <: AbstractFloat}

    detInv = 1 / ((mA[1] * mA[4]) - (mA[2] * mA[2]));

    mLR[1] = sqrt(detInv * mA[4]);
    offDiagVal = -detInv * mA[2] / mLR[1];
    if (trType == "TrLeft")
      mLR[2] = offDiagVal;
      mLR[3] = zero(T);
    elseif (trType == "TrRight")
      mLR[2] = zero(T);
      mLR[3] = offDiagVal;
    end
    mLR[4] = sqrt(detInv * mA[1] - (offDiagVal * offDiagVal));
    
    return mLR;
    
end

function EstDataFreqCedron( vX :: Vector{T}, samplingFreq :: T ) where {T <: AbstractFloat}
    # Exact Frequency Formula for a Pure Real Tone
    # Cedron Dawg
    # https://www.dsprelated.com/showarticle/773.php

    N = length(vX);
    vXX = rfft(vX);
    ~, idxK = findmax(abs2, @view vXX[2:(end - 1)]);
    idxK += 1;
    vXK = @view vXX[(idxK - 1):(idxK + 1)];
    r = exp(-1im * 2π / N);
    vCosB = cos.(2π / N .* [idxK - 2, idxK - 1, idxK]); #<! Zero based like DFT
    num = (-vXK[1] * vCosB[1]) + (vXK[2] * (1 + r) * vCosB[2]) - (vXK[3] * r * vCosB[3]);
    den = -vXK[1] + (vXK[2] * (1 + r)) - vXK[3] * r;
    f = real(acos(num / den)) / (2π);

    return f * samplingFreq;

end

function EstimateSineFreqCedron3Bin( vX :: Vector{T}, samplingFreq :: T ) where {T <: AbstractFloat}
    # Improved Three Bin Exact Frequency Formula for a Pure Real Tone in a DFT
    # Cedron Dawg
    # https://www.dsprelated.com/showarticle/1108.php
    
    numSamples = length(vX);
    
    vXK = rfft(vX);
    
    ~, idxK = findmax(abs2, @view vXK[2:(end - 1)]);
    idxK += 1;
    
    vZ = @view vXK[(idxK - 1):(idxK + 1)];
    
    vR = real.(vZ);
    vI = imag.(vZ);
    iRoot2 = 1 / sqrt(2);
    vBetas = [idxK - 2, idxK - 1, idxK] * (2π / numSamples); #<! Zero based like DFT
    vCosB = cos.(vBetas);
    vSinB = sin.(vBetas);
    
    vA = [iRoot2 * (vR[2] - vR[1]); iRoot2 * (vR[2] - vR[3]); vI[1]; vI[2]; vI[3]];
    vB = [iRoot2 * (vCosB[2] * vR[2] - vCosB[1] * vR[1]); iRoot2 * (vCosB[2] * vR[2] - vCosB[3] * vR[3]); vCosB[1] * vI[1]; vCosB[2] * vI[2]; vCosB[3] * vI[3]];
    vC = [iRoot2 * (vCosB[2] - vCosB[1]); iRoot2 * (vCosB[2] - vCosB[3]); vSinB[1]; vSinB[2]; vSinB[3]];
    
    normC = sqrt(sum(abs2, vC)); #<! Using dot() will be faster
    vP = vC / normC;
    vD = vA .+ vB;
    vK = vD .- (vD' * vP) * vP;
    num = vK' * vB;
    den = vK' * vA;
    
    ratio = max(min(num / den, 1), -1);
    estFreq = acos(ratio) / (2π) * samplingFreq;

    return estFreq;
    
    
end


## Parameters

# Signal Parameters
numSamples      = 100;
samplingFreq    = 1.0; #<! The CRLB is for Normalized Frequency

# Sine Signal Parameters (Non integers divisors of N requires much more realizations).
# sineFreq    = 0.24; #<! Do for [0.05, 0.10, 0.25] For no integer use 0.37.
sineAmp     = 10; #<! High value to allow high SNR
vT          = 0:(numSamples - 1);

# Buffers
mB = zeros(2numSamples + 7, 2); #<! 3 samples (Around the peak)
vS = zeros(numSamples);
vW = zeros(numSamples);
vX = zeros(numSamples);

# Analysis Parameters
numRealizations = 15000;
# SNR of the Analysis (dB)
vSnrdB = LinRange(-15, 50, 150);
# vSnrdB = LinRange(30, 50, 150);

# Methods

tuMethods = (("ML DFT + Opt", (vX, Fₛ) -> EstFrequencyMLOpt(vX, Fₛ, vT, mB)), 
             ("ML Grid + Opt", (vX, Fₛ) -> EstFrequencyMLGridOpt(vX, Fₛ, vT, mB; numGridPts = 512)),
             ("Cedron 3 Bins", (vX, Fₛ) -> EstimateSineFreqCedron3Bin(vX, Fₛ)));

## Load / Generate Data

numNoiseStd = length(vSnrdB);
vNoiseStd   = zeros(numNoiseStd);

for ii = 1:numNoiseStd
    vNoiseStd[ii] = sqrt((sineAmp * sineAmp) / (2 * 10 ^ (vSnrdB[ii] / 10))); 
end

numMethods = length(tuMethods);

tFreqErr = zeros(numRealizations, numNoiseStd, numMethods);

## Analysis

for jj in 1:numNoiseStd
    noiseStd = vNoiseStd[jj];
    for ii in 1:numRealizations
        sineFreq = 0.05 + (0.35 * rand());
        # sineFreq = 0.25;
        angFreq = 2π * (sineFreq / samplingFreq); #<! Make sure (sineFreq / samplingFreq < 0.5)
        sinePhase = 2π * rand();
        vS .= sineAmp .* sin.((angFreq .* vT) .+ sinePhase);
        vW .= vNoiseStd[jj] .* randn(numSamples);
        vX .= vS .+ vW;
        for kk in 1:numMethods
            tFreqErr[ii, jj, kk] = sineFreq - tuMethods[kk][2](vX, samplingFreq);
        end
    end
end

mFreqErr = dropdims(mean(tFreqErr .^ 2; dims = 1), dims = 1);

sineMse     = (sineAmp * sineAmp) / 2;
vNoiseVar   = vNoiseStd .^ 2;
vSnr        = sineMse ./ vNoiseVar;

# CRLB
# 12 -> Sine / Cosine
# 6 -> Exp
vFreqMseCrlb = (12 * samplingFreq * samplingFreq) ./ (((2π) ^ 2) .* vSnr .* ((numSamples ^ 3) - numSamples));

## Display Results

figureIdx += 1;

vTr = [scatter(x = vSnrdB, y = 10 .* log10.(vFreqMseCrlb), mode = "lines", name = "CRLB")];

for ii = 1:numMethods
    push!(vTr, scatter(x = vSnrdB, y = 10 .* log10.(mFreqErr[:, ii]), mode = "lines", name = tuMethods[ii][1]));
end

oLyt = Layout(width = 800, height = 800, title = "Frequency Estimation with $(numRealizations) Realizations", xaxis_title = "SNR [dB]", yaxis_title = "MSE [dB]");

hP = plot(vTr, oLyt);

display(hP);

# hF = Figure(resolution = (800, 800));
# hA = Axis(hF[1, 1], title = "Frequency Estimation", xlabel = "SNR [dB]", ylabel = "MSE [dB]");
# lines!(vSnrdB, 10 .* log10.(vFreqMseCrlb), linewidth = 3, label = "CRLB");
# lines!(vSnrdB, 10 .* log10.(vFreqErr), linewidth = 3, label = "Quadratic Interpolation");
# axislegend();
# display(hF);
# display(hP);

if (exportFigures)
    figFileNme = @sprintf("Figure%04d.png", figureIdx);
    savefig(hP, figFileNme);
end
