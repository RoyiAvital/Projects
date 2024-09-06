# Image Tracking - Optical Flow
# Optical Flow in Lucas Kanade style.
# References:
#   1.  A
# Remarks:
#   1.  Use in Julia as following:
#       -   Move to folder using `cd(raw"<PathToFolder>");`.
#       -   Activate the environment using `] activate .`.
#       -   Instantiate the environment using `] instantiate`.
#   3. 
# TODO:
# 	1.  C
# Release Notes
# - 1.0.000     04/07/2024  Royi Avital RoyiAvital@yahoo.com
#   *   First release.

## Packages

# Internal
using FileIO; #<! Read ImageIO.jl to read images
using LinearAlgebra;
using Printf;
using Random;
# External
# using Convex;
# using Dierckx;
using Interpolations;
# using PlotlyJS;
# using SCS;
using StableRNGs;



## Constants & Configuration
RNG_SEED = 1234;

PROJECT_BASE_FOLDER = "Projects";
AUX_FUN_FOLDER      = "AuxiliaryFunctions";
AUX_FUN             = ["JuliaInit.jl", "JuliaImageProcessing.jl", "JuliaVisualization.jl"];

REG_PATTERN = raw"\\(" * "$(PROJECT_BASE_FOLDER)" * raw")\\\\";

currPwd = pwd();

regMatch = match(Regex(REG_PATTERN), currPwd);
projectFolderPath = currPwd[1:(regMatch.offset + length(PROJECT_BASE_FOLDER))];


for auxFun in AUX_FUN
    include(joinpath(projectFolderPath, AUX_FUN_FOLDER, auxFun));
end


## General Parameters

figureIdx = 0;

exportFigures = false;

oRng = StableRNG(1234);

## Functions

function CalcImgGrad( mI :: Matrix{T} ) where {T <: AbstractFloat}
    
    # Matches 
    numRows = size(mI, 1);
    numCols = size(mI, 2);
    
    mIx = similar(mI);
    mIy = similar(mI);

    for jj ∈ 1:numCols
        for ii ∈ 1:numRows
            if (ii == 1)
                mIy[ii, jj] = mI[2, jj] - mI[1, jj];
            elseif (ii == numRows)
                mIy[ii, jj] = mI[numRows, jj] - mI[numRows - 1, jj];
            else
                mIy[ii, jj] = 0.5 * (mI[ii + 1, jj] - mI[ii - 1, jj]);
            end
            if (jj == 1)
                mIx[ii, jj] = mI[ii, 2] - mI[ii, 1];
            elseif (jj == numCols)
                mIx[ii, jj] = mI[ii, numCols] - mI[ii, numCols - 1];
            else
                mIx[ii, jj] = 0.5 * (mI[ii, jj + 1] - mI[ii, jj - 1]);
            end
        end
    end

    return mIx, mIy;

end

function GenPyrPatch( mPts :: Matrix{T}, mP :: Matrix{T}, vG, vU :: Vector{T}, vV :: Vector{T} ) where {T <: AbstractFloat}

    numRows, numCols = size(mP);

    numGridPts  = length(vG);
    numPts      = size(mPts, 1);

    vXi = zeros(T, numGridPts * numGridPts * numPts);
    vYi = zeros(T, numGridPts * numGridPts * numPts);

    ii = 0;
    for pp ∈ 1:numPts
        shiftX = mPts[pp, 1] - vU[pp];
        shiftY = mPts[pp, 2] - vV[pp];
        for uu ∈ 1:numGridPts
            for vv ∈ 1:numGridPts
                ii += 1;
                # Clamp can be avoided by setting nearest extrapolation mode
                vXi[ii] = clamp(vG[uu] + shiftX, 1, numCols);
                vYi[ii] = clamp(vG[vv] + shiftY, 1, numRows);
            end
        end
    end

    oIntrp = cubic_spline_interpolation((1:numRows, 1:numCols), mP, extrapolation_bc = Flat());

    # @infiltrate

    return reshape(oIntrp.(vXi, vYi), (numGridPts, numGridPts, numPts));

end

function CalcOptFlSpatialDerivative( tP :: Array{T, 3} ) where {T <: AbstractFloat}

    # The Optical Flow equation will be calculated on half pixel shift grid.  
    #
    # ---------------
    # |x * x * x * x|
    # |x * x * x * x|
    # |x * x * x * x|
    # |x * x * x * x|
    # |x * x * x * x|
    # ---------------
    #
    # The patch `tP` is sampled on `x` while the derivatives are evaluated at 
    # `z` (Half pixel shift).
    
    # TODO: Use views
    tPx = tP[:, 2:end, :] .- tP[:, 1:(end - 1), :]; #<! Horizontal derivative
    tPx = T(0.5) .* (tPx[2:end, :, :] .+ tPx[1:(end - 1), :, :]); #<! Averaging with half pixel shift

    tPy = tP[2:end, :, :] .- tP[1:(end - 1), :, :]; #<! Vertical derivative
    tPy = T(0.5) .* (tPy[:, 2:end, :] .+ tPy[:, 1:(end - 1), :]); #<! Averaging with half pixel shift

    return tPx, tPy;

end

function CalcOptFlDerivative( tP1 :: Array{T, 3}, tP2 :: Array{T, 3} ) where {T <: AbstractFloat}

    # The Optical Flow equation will be calculated on half pixel shift grid.  
    #
    # ---------------
    # |x * x * x * x|
    # |x * x * x * x|
    # |x * x * x * x|
    # |x * x * x * x|
    # |x * x * x * x|
    # ---------------
    #
    # The patches (`tP1` / `tP2`) are sampled on `x` while the derivatives are
    # evaluated at `z` (Half pixel shift).
    
    tPx, tPy = CalcOptFlSpatialDerivative(tP2);

    tTmp = tP2 - tP1;
    # TODO: Use views
    tPt = T(0.25) .* (tTmp[1:(end - 1), 1:(end - 1), :] .+ tTmp[2:end, 1:(end - 1), :] .+
                      tTmp[1:(end - 1), 2:end, :] .+ tTmp[2:end, 2:end, :]);
    
    tPxx = tPx .* tPx;
    tPyy = tPy .* tPy;
    tPxy = tPx .* tPy;    
    tPxt = tPx .* tPt;
    tPyt = tPy .* tPt;

    return tPxx, tPyy, tPxy, tPxt, tPyt;

end

function CalcLocalDerivativeSum( tPxx :: Array{T, 3}, tPyy :: Array{T, 3}, tPxy :: Array{T, 3}, tPxt :: Array{T, 3}, tPyt :: Array{T, 3} ) where {T <: AbstractFloat}

    vPxx = dropdims(sum(tPxx, dims = (1, 2)), dims = (1, 2));
    vPyy = dropdims(sum(tPyy, dims = (1, 2)), dims = (1, 2));
    vPxy = dropdims(sum(tPxy, dims = (1, 2)), dims = (1, 2));
    vPxt = dropdims(sum(tPxt, dims = (1, 2)), dims = (1, 2));
    vPyt = dropdims(sum(tPyt, dims = (1, 2)), dims = (1, 2));

    # Output size matches the number of input points
    return vPxx, vPyy, vPxy, vPxt, vPyt;

end


function SolveLSBnd( vX :: Vector{T}, mAA :: Matrix{T}, vAY :: Vector{T}, K :: N, α :: T, ε :: T, αₘ :: T, M :: N ) where {T <: AbstractFloat, N <: Integer}
    # Solving Bounded Least Squares using Projected Gradient Descent.
    # \arg \minₓ 0.5 * || A x - y ||_2^2 \subject to -1 ≤ xᵢ ≤ 1
    # The objective is equivalent to: x' * A' * A * x - 2 (A' * y)' x + y' y
    # Since `y` is constant, one could optimize x' * A' * A * x - 2 (A' * y)'.
    # Hence the objective value can be evaluated without `y`.
    # Should be used when the given initial value `vX` is close to optimal.  
    # As the method is slow to converge, hence useful when only small number
    # of iterations are required.
    # K - Maximum iterations.
    # α - Initial step size.
    # ε - Stopping threshold.
    # αₘ - Minimum step size.
    # M - Maximum backtracking iterations.

    vT = similar(vX);
    vZ = similar(vX);
    
    for ii ∈ 1:K
        copy!(vT, vX); #<! Buffer
        # TODO: Make the calculation of `vG` non allocating
        vG = mAA * vX - vAY; #<! Gradient
        objVal = dot(vX, mAA, vX) - T(2.0) * dot(vAY, vX);
        vZ .= vX .- α .* vG; #<! Gradient Descent Step
        currVal = dot(vZ, mAA, vZ) - T(2.0) * dot(vAY, vZ);
        jj = 0;
        while ( (currVal > objVal) && (α > αₘ) && (jj < M) )
            α *= T(0.5);
            vZ .= vX .- α .* vG;
            currVal = dot(vZ, mAA, vZ) - T(2.0) * dot(vAY, vZ);
            jj += 1;
        end

        vX .= clamp.(vZ, -one(T), one(T)); #<! Projection
        α = T(2.0) * α;

        vT .-= vX;
        maxChange = maximum(abs, vT);

        if (maxChange < ε)
            break;
        end

    end

    return vX;

end

function EstUV( vPxx :: Vector{T}, vPyy :: Vector{T}, vPxy :: Vector{T}, vPxt :: Vector{T}, vPyt :: Vector{T} ) where {T <: AbstractFloat}

    detThr = 1e-8;
    applyCons = false;
    λ = T(0.0);
    K = 100;
    α = T(1e-2);
    ε = T(1e-6);
    αₘ = T(1e-6);
    M = 10;

    numPts = length(vPxx);

    vU = zeros(T, numPts);
    vV = zeros(T, numPts);

    for ii ∈ 1:numPts
        # Solving a 2x2 linear system of a PD matrix
        detVal = max(detThr, vPxx[ii] * vPyy[ii] - vPxy[ii] * vPxy[ii]); #<! Determinant
        vU[ii] = (vPxy[ii] * vPyt[ii] - vPyy[ii] * vPxt[ii]) / detVal;
        vV[ii] = (vPxy[ii] * vPxt[ii] - vPxx[ii] * vPyt[ii]) / detVal;

        if (applyCons) #<! Apply constraint of |vU[ii], vV[ii]| ≤ 1
            mA = zeros(T, 2, 2);
            vY = zeros(T, 2);
            vX = zeros(T, 2);

            mA[1, 1] = vPxx[ii] + λ;
            mA[1, 2] = vPxy[ii];
            mA[2, 1] = vPxy[ii]; #<! Symmetric
            mA[2, 2] = vPyy[ii] + λ;

            vY[1] = -vPxt[ii];
            vY[2] = -vPyt[ii];
            vX[1] = vU[ii];
            vX[2] = vV[ii];
            
            # Solve:
            # \arg \minₓ || A x - y ||_2^2 subject to -1 ≤ xᵢ ≤ 1
            vX = SolveLSBnd(vX, mA, vY, K, α, ε, αₘ, M);
            
            vU[ii] = vX[1];
            vV[ii] = vX[2];
        end
    end

    return vU, vV;


end

function OpticalFlow( mI1 :: Matrix{T}, mI2 :: Matrix{T}, mPts1 :: Matrix{T}, mPts2 :: Matrix{T}, localPatchRadius :: N, numPyr :: N ) where {T <: AbstractFloat, N <: Integer}

    numRows, numCols = size(mI1);
    
    numPts        = size(mPts1, 1);
    localPatchLen = T(2) * localPatchRadius + one(T);    

    # TODO: Make parameters
    σ = T(1.0);
    radFctr = T(3.0);

    vU  = zeros(T, numPts);
    vV  = zeros(T, numPts);
    vU₀ = zeros(T, numPts);
    vV₀ = zeros(T, numPts);

    mOnes = ones(T, numRows, numCols);
    
    for pp ∈ numPyr:-1:0
        # Per scale

        decFactor = 2 ^ pp;

        if (pp > 0)
            # Build the Gaussian Kernel
            σₚ          = σ * decFactor;
            kerRadius   = ceil(N, radFctr * σₚ);
            mK          = GenGaussianKernel(σₚ, (kerRadius, kerRadius));

            mO  = Conv2D(mOnes, mK; convMode = CONV_MODE_SAME);
            mP1 = Conv2D(mI1, mK; convMode = CONV_MODE_SAME) ./ mO;
            mP2 = Conv2D(mI2, mK; convMode = CONV_MODE_SAME) ./ mO;
        else
            mP1 = mI1;
            mP2 = mI2;
        end

        patchGridRad = decFactor * localPatchRadius;
        vG = -patchGridRad:decFactor:patchGridRad;

        # The multi scale process compensate for the previous scale shift.
        # The idea is each scale has maximum shift of 1 pixel in each direction.  
        # Hence, if the largest scale has 1, the next scale is first shifted by it.
        tP1 = GenPyrPatch(mPts1, mP1, vG, vU, vV); #<! Size: (length(vG) x length(vG) x numPts)
        tP2 = GenPyrPatch(mPts1, mP2, vG, vU₀, vV₀); #<! Size: (length(vG) x length(vG) x numPts)

        tPxx, tPyy, tPxy, tPxt, tPyt = CalcOptFlDerivative(tP1, tP2);
        vPxx, vPyy, vPxy, vPxt, vPyt = CalcLocalDerivativeSum(tPxx, tPyy, tPxy, tPxt, tPyt);

        # Disabled refinements for now
        # for rr ∈ 1:numRef #<! Refinements
        #     vUi, vVi = EstUV(vPxx, vPyy, vPxy, vPxt, vPyt); #<! TODO Estimate du, dv per point
        # end
        vUi, vVi = EstUV(vPxx, vPyy, vPxy, vPxt, vPyt); #<! TODO Estimate du, dv per point

        vU .+= decFactor .* vUi; #<! Shift factored by the scale
        vV .+= decFactor .* vVi;

    end

    # Set the estimated coordinates in next image
    for pp ∈ 1:numPts
        # Per point
        mPts2[pp, 1] = mPts1[pp, 1] + vU[pp]
        mPts2[pp, 2] = mPts1[pp, 2] + vV[pp]
    end

    return mPts2;

end





## Parameters

# Data
imgUrl = raw"https://i.imgur.com/WLROWkE.png"; #<! Lena Image 256x256 Gray

# Point to track
vRefPt = [43.3, 33.7];

# Shift
dU = 0.75;
dV = 0.85;

# Optical Flow
localPatchRadius = 2;
numPyr = 2;



## Generate / Load Data

# Read Image
mI = ConvertJuliaImgArray(load(download(imgUrl)));
mI = copy(mI) ./ 255.0;

numRows = size(mI, 1);
numCols = size(mI, 2);

# Data Index
vX = 1:numCols;
vY = 1:numRows;

# Image Interpolation
imgItrp = interpolate(mI, BSpline(Cubic()));
imgItrp = extrapolate(imgItrp, 0);

numPts  = numRows * numCols;

mII = imgItrp.(vY .+ dV, vX' .+ dU); #<! Interpolations works rows / cols
# mII = reshape(mII, (length(vY), length(vX))); #<! Image for template

# hP = DisplayImage(mI);
# display(hP);
# hP = DisplayImage(mII);
# display(hP);

# Extract 2 Patches
mP1 = mI[100:150, 100:150];
# hP = DisplayImage(mP1);
# display(hP);

mP2 = mII[100:150, 100:150];
# hP = DisplayImage(mP2);
# display(hP);

## Analysis

mPts1 = [126.0 126.0; 127.0 126.0; 125.0 126.0; 126.0 125.0; 126.0 127.0]; #<! (x, y) Integer
mPts2 = similar(mPts1);

# OpticalFlow(mI, mII, mPts1, mPts2, localPatchRadius, numPyr)


## Debug

mI1 = copy(mI);
mI2 = copy(mII);

T = eltype(mI1);
N = typeof(numPyr);

numRowsP, numColsP = size(mI1);
    
numPts        = size(mPts1, 1);
localPatchLen = T(2) * localPatchRadius + one(T);

σ = T(1.0);
radFctr = T(3.0);

vU  = zeros(T, numPts);
vV  = zeros(T, numPts);
vU₀ = zeros(T, numPts);
vV₀ = zeros(T, numPts);

mOnes = ones(T, numRows, numCols);

pp = 0;

decFactor = 2 ^ pp;

σₚ          = σ * decFactor;
kerRadius   = ceil(N, radFctr * σₚ);
mK          = GenGaussianKernel(σₚ, (kerRadius, kerRadius));

mO  = Conv2D(mOnes, mK; convMode = CONV_MODE_SAME);
mP1 = Conv2D(mI1, mK; convMode = CONV_MODE_SAME) ./ mO;
mP2 = Conv2D(mI2, mK; convMode = CONV_MODE_SAME) ./ mO;

patchGridRad = decFactor * localPatchRadius;
vG = -patchGridRad:decFactor:patchGridRad;

tP1 = GenPyrPatch(mPts1, mP1, vG, vU, vV); #<! Size: (length(vG) x length(vG) x numPts)
tP2 = GenPyrPatch(mPts1, mP2, vG, vU₀, vV₀); #<! Size: (length(vG) x length(vG) x numPts)

# tPx, tPy = CalcOptFlSpatialDerivative(tP2);

tPxx, tPyy, tPxy, tPxt, tPyt = CalcOptFlDerivative(tP1, tP2);
vPxx, vPyy, vPxy, vPxt, vPyt = CalcLocalDerivativeSum(tPxx, tPyy, tPxy, tPxt, tPyt);

vUi, vVi = EstUV(vPxx, vPyy, vPxy, vPxt, vPyt); #<! TODO Estimate du, dv per point

## 

# hP = DisplayImage(mII);
# # rect: x0, x1, y0, y1
# shapeRect = rect(vTc[1] - 1, vTc[end] - 1, vTr[1] - 1, vTr[end] - 1, line_color = "green", name = "Template Patch");
# add_shape!(hP, shapeRect);
# display(hP);


# vT = mII[vTr, vTc][:];

# vX = copy(vTc);
# vY = copy(vTr);

# numPtsX = vX[end] - vX[1] + 1;
# numPtsY = vY[end] - vY[1] + 1;

# for ii ∈ 1:numIter
#     vI, ∇I = GetPatch(imgItrp, vX, vY);

    
#     Δp = 1;
# end


## Display Results

# figureIdx += 1;

# vTr = Vector{GenericTrace{Dict{Symbol, Any}}}(undef, length(dSolvers));

# # shapeLine = vline(sOptRes.minimizer, line_color = "green", name = "Optimal Value");
# for (ii, methodName) in enumerate(keys(dSolvers))
#     vTr[ii] = scatter(x = 1:numIterations, y = 20 * log10.(abs.(dSolvers[methodName] .- optVal) ./ abs(optVal)), 
#                mode = "lines", text = methodName, name = methodName, line = attr(width = 3.0))
# end
# oLayout = Layout(title = "Objective Function, Condition Number = $(@sprintf("%0.3f", cond(mA)))", width = 600, height = 600, hovermode = "closest",
#                  xaxis_title = "Iteration", yaxis_title = raw"$$\frac{ \left| {f}^{\star} - {f}_{i} \right| }{ \left| {f}^{\star} \right| }$ [dB]$");

# hP = plot(vTr, oLayout);
# display(hP);

# if (exportFigures)
#     figFileNme = @sprintf("Figure%04d.png", figureIdx);
#     savefig(hP, figFileNme);
# end
