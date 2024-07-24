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
using PlotlyJS;
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

## Functions

function CalcImgGrad( mI :: Matrix{T} ) where {T <: AbstractFloat}
    
    # Matches 
    numRows = size(mI, 1);
    numCols = size(mI, 2);
    
    mDx = similar(mI);
    mDy = similar(mI);

    for jj ∈ 1:numCols
        for ii ∈ 1:numRows
            if (ii == 1)
                mDy[ii, jj] = mI[2, jj] - mI[1, jj];
            elseif (ii == numRows)
                mDy[ii, jj] = mI[numRows, jj] - mI[numRows - 1, jj];
            else
                mDy[ii, jj] = 0.5 * (mI[ii + 1, jj] - mI[ii - 1, jj]);
            end
            if (jj == 1)
                mDx[ii, jj] = mI[ii, 2] - mI[ii, 1];
            elseif (jj == numCols)
                mDx[ii, jj] = mI[ii, numCols] - mI[ii, numCols - 1];
            else
                mDx[ii, jj] = 0.5 * (mI[ii, jj + 1] - mI[ii, jj - 1]);
            end
        end
    end

    return mDx, mDy;

end

function GenGaussKernel( σᵤ :: T, σᵥ :: T, tuRadius :: Tuple{N, N} ) where {T <: AbstractFloat, N <: Integer}

    numRows = 2 * tuRadius[1] + 1;
    numCols = 2 * tuRadius[2] + 1;
    mK = zeros(T, numRows, numCols);

    jj = 0;
    for vv ∈ -tuRadius[2]:tuRadius[2]
        ii = 0;
        jj += 1;
        valV = exp(-(vv * vv) / (2 * σᵥ * σᵥ));
        for uu ∈ -tuRadius[1]:tuRadius[1]
            ii += 1;
            valU = exp(-(uu * uu) / (2 * σᵤ * σᵤ));
            mK[ii, jj] = valU * valV;
        end
    end

    mK ./= sum(mK);

    return mK;

end

function GenGaussKernel( σ :: T, tuRadius :: Tuple{N, N} ) where {T <: AbstractFloat, N <: Integer}

    return GenGaussKernel(σ, σ, tuRadius);

end

function GenGaussKernel( σ :: T, kernelRadius :: N ) where {T <: AbstractFloat, N <: Integer}

    return GenGaussKernel(σ, σ, (kernelRadius, kernelRadius));

end

function GenPyrPatch( mPts :: Matrix{T}, mP :: Matrix{T}, vG, vU :: Vector{T}, vV :: Vector{T} ) where {T <: AbstractFloat}

    numRows = size(mP, 1);
    numCols = size(mP, 2);

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
                vXi[ii] = clamp(vG[uu] + shiftX, 1, N);
                vYi[ii] = clamp(vG[vv] + shiftY, 1, M);
            end
        end
    end

    oIntrp = cubic_spline_interpolation((1:numRows, 1:numCols), mP, extrapolation_bc = Flat());

    return oIntrp.(vXi, vYi);

end



function OpticalFlow( mI1 :: Matrix{T}, mI2 :: Matrix{T}, mPts1 :: Matrix{T}, mPts2 :: Matrix{T}, localPatchRadius :: N, numPyr :: N ) where {T <: AbstractFloat, N <: Integer}

    numPts        = size(mPts1, 1);
    localPatchLen = T(2) * localPatchRadius + one(T);    

    vU  = zeros(T, numPts);
    vV  = zeros(T, numPts);
    vU₀ = zeros(T, numPts);
    vV₀ = zeros(T, numPts);
    
    for pp ∈ numPyr:-1:0

        decFactor = 2 ^ pp;

        if (pp > 0)
            # Build the Gaussian Kernel
            σₚ          = σ * pyrFctr;
            kerRadius   = ceil(radFctr * σₚ);
            mK          = GenGaussKernel(σₚ, kerRadius);

            mO  = Conv2D(mOnes, mK; convMode = CONV_MODE_SAME);
            mP1 = Conv2D(mI1, mK; convMode = CONV_MODE_SAME) ./ mO;
            mP2 = Conv2D(mI2, mK; convMode = CONV_MODE_SAME) ./ mO;
        else
            mP1 = mI1;
            mP2 = mI2;
        end

        patchGridRad = decFactor * localPatchRadius;
        vG = -patchGridRad:decFactor:patchGridRad;

        tP1 = GenPyrPatch(mPts1, mP1, vG, vU, vV);
        tP2 = GenPyrPatch(mPts2, mP1, vG, vU₀, vV₀);

        tD = CalcOptDerivative(tP1, tP2); #<! TODO
        tS = LocalDerivativeSum(tD); #<! TODO

        for rr ∈ 1:numRef #<! Refinements
            vUi, vVi = EstUV(tS); #<! TODO Estimate du, dv per point
        end

        vU .+= decFactor .* vUi;
        vV .+= decFactor .* vVi;

    end

    # Set the estimated coordinates in next image
    for pp ∈ 1:numPts
        mP2[pp, 1] = mP1[pp, 1] + vU[pp]
        mP2[pp, 2] = mP1[pp, 2] + vV[pp]
    end

end





## Parameters

# Data
imgUrl = raw"https://i.imgur.com/WLROWkE.png"; #<! Lena Image 256x256 Gray

# Warp Matrix
αᵤ = 1.0 #<! Scaling
αᵥ = 1.0 #<! Scaling
θ = (5.0 / 180.0) * π; #<! Rotation
tᵤ = 0.95;
tᵥ = 0.25;
mW = [1.0 0.0 tᵤ; 0.0 1.0 tᵥ];
mW = [αᵤ * cos(θ) sin(θ) 0.35; -sin(θ) αᵥ * cos(θ) 0.0];

# Template
vTr = 125:155; #<! Rows
vTc = 110:140; #<! Columns

# Solver
numIter = 100;
mP = [0.0 0.0 vTc[1]; 0.0 0.0 vTr[1]];

## Generate / Load Data
oRng = StableRNG(1234);

# Read Image
mI = ConvertJuliaImgArray(load(download(imgUrl)));
mI = copy(mI) ./ 255.0;

numRows = size(mI, 1);
numCols = size(mI, 2);

vX = 1:numCols;
vY = 1:numRows;

imgItrp = interpolate(mI, BSpline(Cubic()));
imgItrp = extrapolate(imgItrp, 0);

## Analysis
numPts  = numRows * numCols;
mY      = zeros(2, numPts);

ii = 0;
for valX ∈ vX, valY ∈ vY
    # valY is the inner (Running on the rows)
    global ii += 1;
    mY[:, ii] = mW * [valX, valY, 1.0];
end

mII = imgItrp.(mY[2, :], mY[1, :]); #<! Interpolations works rows / cols
mII = reshape(mII, (length(vY), length(vX))); #<! Image for template

DisplayImage(mI)
DisplayImage(mII)

## 

hP = DisplayImage(mII);
# rect: x0, x1, y0, y1
shapeRect = rect(vTc[1] - 1, vTc[end] - 1, vTr[1] - 1, vTr[end] - 1, line_color = "green", name = "Template Patch");
add_shape!(hP, shapeRect);
display(hP);


vT = mII[vTr, vTc][:];

vX = copy(vTc);
vY = copy(vTr);

numPtsX = vX[end] - vX[1] + 1;
numPtsY = vY[end] - vY[1] + 1;

for ii ∈ 1:numIter
    vI, ∇I = GetPatch(imgItrp, vX, vY);

    
    Δp = 1;
end


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
