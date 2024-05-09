# Image Tracking - Lucas Kanade Image Alignment
# Template Tracking via Lucas Kanade's Algorithm.
# References:
#   1.  Simon Baker, Iain Matthews - Lucas Kanade 20 Years On: A Unifying Framework.
# Remarks:
#   1.  Use in Julia as following:
#       -   Move to folder using `cd(raw"<PathToFolder>");`.
#       -   Activate the environment using `] activate .`.
#       -   Instantiate the environment using `] instantiate`.
#   3. 
# TODO:
# 	1.  C
# Release Notes
# - 1.0.000     08/05/2024  Royi Avital RoyiAvital@yahoo.com
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


function ImgWarp( imgInterp, mW :: Matrix{T}, vX :: Vector{T}, vY :: Vector{T} ) where { T <: AbstractFloat }

    numPts  = numRows * numCols;
    mY      = zeros(2, numPts);
    
    ii = 0;
    for valX ∈ vX, valY ∈ vY
        # valY is the inner (Running on the rows)
        ii += 1;
        mY[:, ii] = mW * [valX, valY, 1.0];
    end

    # mY = mW * mX; #<! Calculate the warped coordinates
    mV = imgInterp.(mY[2, :], mY[2, :]); #<! Interpolations works rows / cols
    
    return reshape(mV, size(mX));

end

function GenW( mP :: Matrix{T} ) where { T <: AbstractFloat }

    mW = copy(mP);
    mW[1, 1] += 1;
    mW[2, 2] += 1;

    return mW;

end

function GetPatch( imgInterp, mX :: Matrix{T} ) where { T <: AbstractFloat }

    numPts = size(mX, 2);
    ∇I = zeros(2, numPts);

    for ii ∈ numPts
        ∇I[:, ii] = Interpolations.gradient(imgInterp, mX[2, ii], mX[1, ii]);
    end

    return imgInterp.(mX[2, :], mX[1, :]), ∇I;

end

function GetPatch( imgInterp, vX :: Vector{T}, vY :: Vector{T} ) where { T <: AbstractFloat }

    numPts = length(vX) * length(vY);
    vT = zeros(numPts);
    ∇I = zeros(2, numPts);

    ii = 0;
    for valX ∈ vX, valY ∈ vY
        ii += 1;
        vT[ii]      = imgInterp(valY, valX);
        ∇I[:, ii]   = Interpolations.gradient(imgInterp, valY, valX);
    end

    return vT, ∇I;

end

function ∇I∇W( ∇I :: Matrix{T} vX :: Vector{T}, vY :: Vector{T} ) where { T <: AbstractFloat }

    numPts = size(∇I, 2);
    
    ∇I∇WP   = zeros(1, 6);
    ∇WP     = zeros(2, 6);
    ∇WP[[11, 12]] .= 1.0;
    
    ii = 0;
    for valX ∈ vX, valY ∈ vY
        ii += 1;
        ∇WP[[1, 4]] .= valX;
        ∇WP[[6, 8]] .= valY;

        ∇I∇WP += ∇I[:, ii]' * ∇WP
    end


end

function WeightedErr( vI :: Vector{T}, vT :: Vector{T}, ∇I :: Matrix{T} vX :: Vector{T}, vY :: Vector{T} ) where { T <: AbstractFloat }

    ∇I∇WP   = zeros(1, 6);
    ∇WP     = zeros(2, 6);
    ∇WP[[11, 12]] .= 1.0;
    
    ii = 0;
    for valX ∈ vX, valY ∈ vY
        ii += 1;
        ∇WP[[1, 4]] .= valX;
        ∇WP[[6, 8]] .= valY;

        ∇I∇WP = ∇I[:, ii]' * ∇WP;

        ∇I∇WP += ∇I[:, ii]' * ∇WP * (vI[ii] - vT[ii]);
        mW = 1; #<! Merge into single loop: ∇I∇WP and weighted error
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