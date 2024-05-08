# Image Tracking - Lucas Kanade Image Alignment
# Template Tracking via Lucas Kanade's Algorithm.
# References:
#   1.  Simon Baker, Iain Matthews - Lucas Kanade 20 Years On: A Unifying Framewor
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
using Interpolation;
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


function ImgWarp( imgInterp, mW :: Matrix{T}, mX :: Matrix{T} ) where { T <: AbstractFloat }

    mY = mW * mX; #<! Calculate the warped coordinates
    mV = imgInterp(mY[:, 1], mY[:, 2]);
    
    return reshape(mV, size(mX));

end


## Parameters

# Data
imgUrl = raw"https://i.imgur.com/WLROWkE.png"; #<! Lena Image 256x256 Gray

# Warp Matrix
mW = [1.0 0.0 0.85; 0.0 1.0 0.0];

## Generate / Load Data
oRng = StableRNG(1234);

# Read Image
mI = ConvertJuliaImgArray(load(download(imgUrl)));
mI = copy(mI) ./ 255.0;

numRows = size(mI, 1);
numCols = size(mI, 2);

vX = 0:(numCols - 1);
vY = 0:(numRows - 1);

imgItrp = interpolate(mI);

## Analysis




## Display Results

figureIdx += 1;

vTr = Vector{GenericTrace{Dict{Symbol, Any}}}(undef, length(dSolvers));

# shapeLine = vline(sOptRes.minimizer, line_color = "green", name = "Optimal Value");
for (ii, methodName) in enumerate(keys(dSolvers))
    vTr[ii] = scatter(x = 1:numIterations, y = 20 * log10.(abs.(dSolvers[methodName] .- optVal) ./ abs(optVal)), 
               mode = "lines", text = methodName, name = methodName, line = attr(width = 3.0))
end
oLayout = Layout(title = "Objective Function, Condition Number = $(@sprintf("%0.3f", cond(mA)))", width = 600, height = 600, hovermode = "closest",
                 xaxis_title = "Iteration", yaxis_title = raw"$$\frac{ \left| {f}^{\star} - {f}_{i} \right| }{ \left| {f}^{\star} \right| }$ [dB]$");

hP = plot(vTr, oLayout);
display(hP);

if (exportFigures)
    figFileNme = @sprintf("Figure%04d.png", figureIdx);
    savefig(hP, figFileNme);
end