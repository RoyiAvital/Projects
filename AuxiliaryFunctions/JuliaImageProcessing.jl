# Projects Code - Julia Image Processing
# Set of functions for Image Processing.
# References:
#   1.  
# Remarks:
#   1.  A
# TODO:
# 	1.  Add convolution in frequency domain as in DSP `Q90036`.
#       It should include an auxiliary function: `GenWorkSpace()` for teh buffers.  
#       It should also be optimized for `rfft()`.
#   2.  
# Release Notes
# - 1.4.000     05/07/2023  Royi Avital RoyiAvital@yahoo.com
#   *   Added a convolution function with support for `CONV_MODE_SAME`.
# - 1.3.000     29/06/2023  Royi Avital RoyiAvital@yahoo.com
#   *   Added functions to apply linear color conversion.
# - 1.2.000     29/06/2023  Royi Avital RoyiAvital@yahoo.com
#   *   Added functions to generate Convolution Matrix (Sparse).
#   *   Added functions to generate Filter Matrix (Sparse).
#   *   Added a function to calculate the image Laplace Operator.
#   *   Added a function for in place padding.
# - 1.1.000     23/11/2023  Royi Avital RoyiAvital@yahoo.com
#   *   Added 2D Convolution.
# - 1.0.000     09/07/2023  Royi Avital RoyiAvital@yahoo.com
#   *   First release.

## Packages

# Internal
using SparseArrays;

# External
using ColorTypes;
using LoopVectorization;
using StaticKernels;

## Constants & Configuration

include("./JuliaInit.jl");

## Functions

function ConvertJuliaImgArray(mI :: Matrix{<: Color{T, N}}) where {T, N}
    # According to ColorTypes data always in order RGBA

    # println("RGB");
    
    numRows, numCols = size(mI);
    numChannels = N;
    dataType = T.types[1];

    mO = permutedims(reinterpret(reshape, dataType, mI), (2, 3, 1));

    if numChannels > 1
        mO = permutedims(reinterpret(reshape, dataType, mI), (2, 3, 1));
    else
        mO = reinterpret(reshape, dataType, mI);
    end

    return mO;

end

function ConvertJuliaImgArray(mI :: Matrix{<: Color{T, 1}}) where {T}
    # According to ColorTypes data always in order RGBA
    # Single Channel Image (numChannels = 1)
    
    # println("Gray");
    
    numRows, numCols = size(mI);
    dataType = T.types[1];

    mO = reinterpret(reshape, dataType, mI);

    return mO;

end

function PadArray!( mB :: Matrix{T}, mA :: Matrix{T}, tuPadRadius :: Tuple{N, N}; padMode :: PadMode, padValue :: T = zero(T) ) where {T <: Real, N <: Signed}

    numRowsA, numColsA = size(mA);
    numRowsB, numColsB = (numRowsA, numColsA) .+ (N(2) .* tuPadRadius);
    
    if (padMode == PAD_MODE_CONSTANT)
        mB .= padValue;
        mB[(tuPadRadius[1] + one(N)):(numRowsB - tuPadRadius[1]), (tuPadRadius[2] + one(N)):(numColsB - tuPadRadius[2])] .= mA;
    end

    if (padMode == PAD_MODE_REPLICATE)
        for jj in 1:numColsB
            nn = clamp(jj - tuPadRadius[2], one(N), numColsA);
            for ii in 1:numRowsB
                mm = clamp(ii - tuPadRadius[1], one(N), numRowsA);
                mB[ii, jj] = mA[mm, nn];
            end
        end
    end

    if (padMode == PAD_MODE_SYMMETRIC)
        for jj in 1:numColsB
            nn = ifelse(jj - tuPadRadius[2] < one(N), tuPadRadius[2] - jj + one(N), jj - tuPadRadius[2]);
            nn = ifelse(jj - tuPadRadius[2] > numColsA, numColsA - (jj - (numColsA + tuPadRadius[2]) - 1), nn);
            for ii in 1:numRowsB
                mm = ifelse(ii - tuPadRadius[1] < one(N), tuPadRadius[1] - ii + one(N), ii - tuPadRadius[1]);
                mm = ifelse(ii - tuPadRadius[1] > numRowsA, numRowsA - (ii - (numRowsA + tuPadRadius[1]) - one(N)), mm);
                mB[ii, jj] = mA[mm, nn];
            end
        end
    end

    if (padMode == PAD_MODE_REFLECT) #<! Mirror
        for jj in 1:numColsB
            nn = ifelse(jj - tuPadRadius[2] < one(N), tuPadRadius[2] - jj + N(2), jj - tuPadRadius[2]);
            nn = ifelse(jj - tuPadRadius[2] > numColsA, numColsA - (jj - (numColsA + tuPadRadius[2])), nn);
            for ii in 1:numRowsB
                mm = ifelse(ii - tuPadRadius[1] < one(N), tuPadRadius[1] - ii + N(2), ii - tuPadRadius[1]);
                mm = ifelse(ii - tuPadRadius[1] > numRowsA, numRowsA - (ii - (numRowsA + tuPadRadius[1])), mm);
                mB[ii, jj] = mA[mm, nn];
            end
        end
    end

    if (padMode == PAD_MODE_CIRCULAR)
        for jj in 1:numColsB
            nn = ifelse(jj - tuPadRadius[2] < one(N), numColsA + jj - tuPadRadius[2], jj - tuPadRadius[2]);
            nn = ifelse(jj - tuPadRadius[2] > numColsA, jj - numColsA - tuPadRadius[2], nn);
            for ii in 1:numRowsB
                mm = ifelse(ii - tuPadRadius[1] < 1, numRowsA + ii - tuPadRadius[1], ii - tuPadRadius[1]);
                mm = ifelse(ii - tuPadRadius[1] > numRowsA, ii - numRowsA - tuPadRadius[1], mm);
                mB[ii, jj] = mA[mm, nn];
            end
        end
    end

    return mB;

end

function PadArray( mA :: Matrix{T}, tuPadRadius :: Tuple{N, N}; padMode :: PadMode, padValue :: T = zero(T) ) where {T <: Real, N <: Signed}
    # Works on Matrix
    # TODO: Support padding larger then the input.
    # TODO: Extend ot Array{T, 3}.
    # TODO: Create non allocating variant.

    numRowsA, numColsA = size(mA);
    numRowsB, numColsB = (numRowsA, numColsA) .+ (N(2) .* tuPadRadius);
    mB = Matrix{T}(undef, numRowsB, numColsB);

    mB = PadArray!(mB, mA, tuPadRadius; padMode = padMode, padValue = padValue);

    return mB;
    
    

end

function Conv2D( mI :: Matrix{T}, mK :: Matrix{T}; convMode :: ConvMode = CONV_MODE_FULL ) where {T <: AbstractFloat}
    
    if (convMode == CONV_MODE_FULL)
        mO = Matrix{T}(undef, size(mI) .+ size(mK) .- (1, 1));
    elseif (convMode == CONV_MODE_SAME)
        mO = Matrix{T}(undef, size(mI));
    elseif (convMode == CONV_MODE_VALID)
        mO = Matrix{T}(undef, size(mI) .- size(mK) .+ (1, 1));
    end

    Conv2D!(mO, mI, mK; convMode = convMode);
    
    return mO;

end

function Conv2D!( mO :: Matrix{T}, mI :: Matrix{T}, mK :: Matrix{T}; convMode :: ConvMode = CONV_MODE_FULL ) where {T <: AbstractFloat}

    if (convMode == CONV_MODE_FULL)
        _Conv2D!(mO, mI, mK);
    elseif (convMode == CONV_MODE_SAME)
        _Conv2DSame!(mO, mI, mK);
    elseif (convMode == CONV_MODE_VALID)
        _Conv2DValid!(mO, mI, mK);
    end

    return mO;
    
end

function _Conv2D!( mO :: Matrix{T}, mI :: Matrix{T}, mK :: Matrix{T} ) where {T <: AbstractFloat}
    # Full Convolution

    numRowsI, numColsI = size(mI);
    numRowsK, numColsK = size(mK);

    for jj ∈ 1:(numColsK - 1), ii ∈ 1:(numRowsK - 1) #<! Top Left
        sumVal = zero(T);
        for nn ∈ 1:numColsK, mm ∈ 1:numRowsK
            ib0 = (jj >= nn) && (ii >= mm);
            @inbounds oa = ib0 ? mI[ii - mm + 1, jj - nn + 1] : zero(T);
            @inbounds sumVal += mK[mm, nn] * oa;
        end
        mO[ii, jj] = sumVal;
    end

    for jj ∈ 1:(numColsK - 1), ii ∈ numRowsK:(numRowsI - 1) #<! Middle Left
        sumVal = zero(T);
        for nn ∈ 1:numColsK, mm ∈ 1:numRowsK
            ib0 = (jj >= nn);
            @inbounds oa = ib0 ? mI[ii - mm + 1, jj - nn + 1] : zero(T);
            @inbounds sumVal += mK[mm, nn] * oa;
        end
        mO[ii, jj] = sumVal;
    end

    for jj ∈ 1:(numColsK - 1), ii ∈ numRowsI:(numRowsI + numRowsK - 1) #<! Bottom Left
        sumVal = zero(T);
        for nn ∈ 1:numColsK, mm ∈ 1:numRowsK
            ib0 = (jj >= nn) && (ii < numRowsI + mm);;
            @inbounds oa = ib0 ? mI[ii - mm + 1, jj - nn + 1] : zero(T);
            @inbounds sumVal += mK[mm, nn] * oa;
        end
        mO[ii, jj] = sumVal;
    end

    for jj ∈ numColsK:(numColsI - 1), ii ∈ 1:(numRowsK - 1) #<! Top Middle
        sumVal = zero(T);
        for nn ∈ 1:numColsK, mm ∈ 1:numRowsK
            ib0 = (ii >= mm);
            @inbounds oa = ib0 ? mI[ii - mm + 1, jj - nn + 1] : zero(T);
            @inbounds sumVal += mK[mm, nn] * oa;
        end
        mO[ii, jj] = sumVal;
    end

    for jj ∈ numColsK:(numColsI - 1)
        @turbo for ii ∈ numRowsK:(numRowsI - 1) #<! Middle Middle
            sumVal = zero(T);
            for nn ∈ 1:numColsK, mm ∈ 1:numRowsK
                @inbounds sumVal += mK[mm, nn] * mI[ii - mm + 1, jj - nn + 1];
            end
            mO[ii, jj] = sumVal;
        end        
    end

    for jj ∈ numColsK:(numColsI - 1), ii ∈ numRowsI:(numRowsI + numRowsK - 1) #<! Bottom Middle
        sumVal = zero(T);
        for nn ∈ 1:numColsK, mm ∈ 1:numRowsK
            ib0 = (ii < numRowsI + mm);;
            @inbounds oa = ib0 ? mI[ii - mm + 1, jj - nn + 1] : zero(T);
            @inbounds sumVal += mK[mm, nn] * oa;
        end
        mO[ii, jj] = sumVal;
    end

    for jj ∈ numColsI:(numColsI + numColsK - 1), ii ∈ 1:(numRowsK - 1) #<! Top Right
        sumVal = zero(T);
        for nn ∈ 1:numColsK, mm ∈ 1:numRowsK
            ib0 = (jj < numColsI + nn) && (ii >= mm);
            @inbounds oa = ib0 ? mI[ii - mm + 1, jj - nn + 1] : zero(T);
            @inbounds sumVal += mK[mm, nn] * oa;
        end
        mO[ii, jj] = sumVal;
    end

    for jj ∈ numColsI:(numColsI + numColsK - 1), ii ∈ numRowsK:(numRowsI - 1) #<! Middle Right
        sumVal = zero(T);
        for nn ∈ 1:numColsK, mm ∈ 1:numRowsK
            ib0 = (jj < numColsI + nn);
            @inbounds oa = ib0 ? mI[ii - mm + 1, jj - nn + 1] : zero(T);
            @inbounds sumVal += mK[mm, nn] * oa;
        end
        mO[ii, jj] = sumVal;
    end

    for jj ∈ numColsI:(numColsI + numColsK - 1), ii ∈ numRowsI:(numRowsI + numRowsK - 1) #<! Bottom Right
        sumVal = zero(T);
        for nn ∈ 1:numColsK, mm ∈ 1:numRowsK
            ib0 = (jj < numColsI + nn) && (ii < numRowsI + mm);
            @inbounds oa = ib0 ? mI[ii - mm + 1, jj - nn + 1] : zero(T);
            @inbounds sumVal += mK[mm, nn] * oa;
        end
        mO[ii, jj] = sumVal;
    end

end

function _Conv2DSame!( mO :: Matrix{T}, mI :: Matrix{T}, mK :: Matrix{T} ) where {T <: AbstractFloat}
    # Matches MATLAB
    # Assumes size(mK) <= size(mI)

    numRowsI, numColsI = size(mI);
    numRowsK, numColsK = size(mK);

    ancY = (numRowsK ÷ 2) + 1; #<! Anchor (Aligned with mI[ii, jj])
    ancX = (numColsK ÷ 2) + 1; #<! Anchor (Aligned with mI[ii, jj])

    # Mapping: `mm -> ii + (mm - ancY)` is correlation (Works with `ancY = (numRowsK + 1) ÷ 2;`, mm ∈ 1:numRowsK)
    # Mapping: `mm -> ii + (ancY - mm)` is convolution (For odd size, Symmetric)
    # Mapping: `mm -> ii + (ancY - mm)` is convolution (Works with `ancY = (numRowsK ÷ 2) + 1;`, mm ∈ numRowsK:-1:1)

    for jj ∈ 1:numColsI
        @fastmath @inbounds @simd for ii in 1:numRowsI
            sumVal = zero(T);
            for nn ∈ numColsK:-1:1, mm ∈ numRowsK:-1:1
                inBnd = (jj + ancX <= numColsI + nn) && (jj + ancX > nn) && (ii + ancY <= numRowsI + mm) && (ii + ancY > mm);
                @inbounds valI = inBnd ? mI[ii + ancY - mm, jj + ancX - nn] : zero(T);
                @inbounds sumVal += mK[mm, nn] * valI;
            end
            mO[ii, jj] = sumVal;
        end
    end

end

function _Conv2DValid!( mO :: Matrix{T}, mI :: Matrix{T}, mK :: Matrix{T} ) where {T <: AbstractFloat}

    numRowsI, numColsI = size(mI);
    numRowsK, numColsK = size(mK);

    for jj ∈ 1:(numColsI - numColsK + 1)
        @turbo for ii in 1:(numRowsI - numRowsK + 1)
            sumVal = zero(T);
            for nn ∈ 1:numColsK, mm ∈ 1:numRowsK
                @inbounds sumVal += mK[mm, nn] * mI[ii - mm + numRowsK, jj - nn + numColsK];
            end
            mO[ii, jj] = sumVal;
        end
    end

end


function GenConvMtx( vK :: AbstractVector{T}, numElements :: N; convMode :: ConvMode = CONV_MODE_FULL ) where {T <: AbstractFloat, N <: Integer}
    
    kernelLength = length(vK);

    if (convMode == CONV_MODE_FULL)
        rowIdxFirst = 1;
        rowIdxLast  = numElements + kernelLength - 1;
        outputSize  = numElements + kernelLength - 1;
    elseif (convMode == CONV_MODE_SAME)
        rowIdxFirst = 1 + floor(kernelLength / 2);
        rowIdxLast  = rowIdxFirst + numElements - 1;
        outputSize  = numElements;
    elseif (convMode == CONV_MODE_VALID)
        rowIdxFirst = kernelLength;
        rowIdxLast  = (numElements + kernelLength - 1) - kernelLength + 1;
        outputSize  = numElements - kernelLength + 1;
    end

    mtxIdx = 0;
    vI = ones(N, numElements * kernelLength);
    vJ = ones(N, numElements * kernelLength);
    vV = zeros(T, numElements * kernelLength);

    for jj ∈ 1:numElements
        for ii ∈ 1:kernelLength
            if((ii + jj - 1 >= rowIdxFirst) && (ii + jj - 1 <= rowIdxLast))
                # Valid output matrix row index
                mtxIdx += 1;
                vI[mtxIdx] = ii + jj - rowIdxFirst;
                vJ[mtxIdx] = jj;
                vV[mtxIdx] = vK[ii];
            end
        end
    end

    mK = sparse(vI, vJ, vV, outputSize, numElements);

    return mK;

end


function GenConvMtx( mH :: AbstractMatrix{T}, numRows :: N, numCols :: N; convMode :: ConvMode = CONV_MODE_FULL ) where {T <: AbstractFloat, N <: Integer}
    
    numColsKernel   = size(mH, 2);
    numBlockMtx     = numColsKernel;

    lBlockMtx = Vector{SparseMatrixCSC{T, N}}(undef, numBlockMtx);

    for ii ∈ 1:numBlockMtx
        lBlockMtx[ii] = GenConvMtx(mH[:, ii], numRows; convMode = convMode);
    end

    if (convMode == CONV_MODE_FULL)
        # For convolution shape - 'full' the Doubly Block Toeplitz Matrix
        # has the first column as its main diagonal.
        diagIdx     = 0;
        numRowsKron = numCols + numColsKernel - 1;
    elseif (convMode == CONV_MODE_SAME)
        # For convolution shape - 'same' the Doubly Block Toeplitz Matrix
        # has the first column shifted by the kernel horizontal radius.
        diagIdx     = floor(numColsKernel / 2);
        numRowsKron = numCols;
    elseif (convMode == CONV_MODE_VALID)
        # For convolution shape - 'valid' the Doubly Block Toeplitz Matrix
        # has the first column shifted by the kernel horizontal length.
        diagIdx     = numColsKernel - 1;
        numRowsKron = numCols - numColsKernel + 1;
    end

    vI = ones(N, min(numRowsKron, numCols));
    mK = kron(spdiagm(numRowsKron, numCols, diagIdx => vI), lBlockMtx[1]);

    for ii ∈ 2:numBlockMtx
        diagIdx = diagIdx - 1;
        mK = mK + kron(spdiagm(numRowsKron, numCols, diagIdx => vI), lBlockMtx[ii]);
    end

    return mK;

end

function GenFilterMtx( mH :: AbstractMatrix{T}, numRows :: N, numCols :: N; filterMode :: FilterMode = FILTER_MODE_CONVOLUTION, boundaryMode :: BoundaryMode = BND_MODE_REPLICATE ) where {T <: AbstractFloat, N <: Integer}
    
    numElementsImage    = numRows * numCols;
    numRowsKernel       = size(mH, 1);
    numColsKernel       = size(mH, 2);
    numElementsKernel   = numRowsKernel * numColsKernel;

    vR = repeat(1:numElementsImage, inner = numElementsKernel); #<! `repeat()` returns a Vector from the range
    vC = zeros(Int64, numElementsImage * numElementsKernel);
    vV = zeros(T, numElementsImage * numElementsKernel);

    kernelRadiusV = numRowsKernel ÷ 2;
    kernelRadiusH = numColsKernel ÷ 2;

    if (filterMode == FILTER_MODE_CONVOLUTION)
        mHH = rot180(mH); #<! Code is written for correlation
    elseif (filterMode == FILTER_MODE_CORRELATION)
        mHH = mH; #<! Code is written for correlation
    end

    if (boundaryMode == BND_MODE_CIRCULAR)
        CalcPxShift = _CalcPxShiftCircular;
    elseif (boundaryMode == BND_MODE_REPLICATE)
        CalcPxShift = _CalcPxShiftReplicate;
    elseif (boundaryMode == BND_MODE_SYMMETRIC)
        CalcPxShift = _CalcPxShiftSymmetric;
    elseif (boundaryMode == BND_MODE_ZEROS)
        CalcPxShift = _CalcPxShiftZeros;
    end

    pxIdx    = zero(N);
    elmntIdx = zero(N);

    for jj ∈ 1:numCols
        for ii ∈ 1:numRows
            pxIdx += 1;
            for ll ∈ -kernelRadiusH:kernelRadiusH
                for kk ∈ -kernelRadiusV:kernelRadiusV
                    elmntIdx += 1;
                    
                    # Pixel Index Shift such that pxIdx + pxShift is the linear
                    # index of the pixel in the image
                    pxShift, isValid = CalcPxShift(ii, jj, ll, kk, numRows, numCols, numElementsImage);

                    if (isValid)
                        valH = mHH[kk + kernelRadiusV + one(N), ll + kernelRadiusH + one(N)];
                    else
                        pxShift = zero(N);
                        valH    = zero(T);
                    end
                    
                    vC[elmntIdx] = pxIdx + pxShift;
                    vV[elmntIdx] = valH;
                    
                end
            end
        end
    end
    
    mK = sparse(vR, vC, vV, numElementsImage, numElementsImage);
    dropzeros!(mK); #<! Julia inserts zeros, this removes

    return mK;

end

function _CalcPxShiftCircular( ii :: N, jj :: N, ll :: N, kk :: N, numRows :: N, numCols :: N, numElements :: N ) where {N <: Integer}

    # Pixel Index Shift such that pxIdx + pxShift is the linear
    # index of the pixel in the image
    pxShift = (ll * numRows) + kk;
    
    if((ii + kk) > numRows)
        pxShift -= numRows;
    end
    
    if((ii + kk) < one(N))
        pxShift += numRows;
    end
    
    if((jj + ll) > numCols)
        pxShift -= numElements;
    end
    
    if(jj + ll < one(N))
        pxShift += numElements;
    end

    return pxShift, true;

end

function _CalcPxShiftReplicate( ii :: N, jj :: N, ll :: N, kk :: N, numRows :: N, numCols :: N, numElements :: N ) where {N <: Integer}

    # Pixel Index Shift such that pxIdx + pxShift is the linear
    # index of the pixel in the image
    pxShift = (ll * numRows) + kk;
    
    if((ii + kk) > numRows)
        pxShift = pxShift - (ii + kk - numRows);
    end
    
    if((ii + kk) < one(N))
        pxShift = pxShift + (one(N) - (ii + kk));
    end
    
    if((jj + ll) > numCols)
        pxShift = pxShift - ((jj + ll - numCols) * numRows);
    end
    
    if(jj + ll < one(N))
        pxShift = pxShift + ((one(N) - (jj + ll)) * numRows);
    end

    return pxShift, true;

end

function _CalcPxShiftSymmetric( ii :: N, jj :: N, ll :: N, kk :: N, numRows :: N, numCols :: N, numElements :: N ) where {N <: Integer}

    # Pixel Index Shift such that pxIdx + pxShift is the linear
    # index of the pixel in the image
    pxShift = (ll * numRows) + kk;
    
    if((ii + kk) > numRows)
        pxShift = pxShift - (N(2) * (ii + kk - numRows) - one(N));
    end
    
    if((ii + kk) < one(N))
        pxShift = pxShift + (N(2) * (one(N) - (ii + kk)) - one(N));
    end
    
    if((jj + ll) > numCols)
        pxShift = pxShift - ((N(2) * (jj + ll - numCols) - one(N)) * numRows);
    end
    
    if(jj + ll < one(N))
        pxShift = pxShift + ((N(2) * (one(N) - (jj + ll)) - one(N)) * numRows);
    end

    return pxShift, true;

end

function _CalcPxShiftZeros( ii :: N, jj :: N, ll :: N, kk :: N, numRows :: N, numCols :: N, numElements :: N ) where {N <: Integer}

    # Pixel Index Shift such that pxIdx + pxShift is the linear
    # index of the pixel in the image
    pxShift = (ll * numRows) + kk;
    
    isValid = ((ii + kk <= numRows) && (ii + kk >= one(N)) && (jj + ll <= numCols) && (jj + ll >= one(N)));

    return pxShift, isValid;

end


function CalcImageLaplacian!( mO :: Matrix{T}, mI :: Matrix{T}, mB :: Matrix{T}, mV1 :: Matrix{T}, mV2 :: Matrix{T}, mKₕ :: Matrix{T}, mKᵥ :: Matrix{T} ) where {T <: AbstractFloat}

    _Conv2DValid!(mV1, mI, mKₕ);
    _Conv2D!(mO, mV1, rot180(mKₕ)); #<! TODO: Remove allocation

    _Conv2DValid!(mV2, mI, mKᵥ);
    _Conv2D!(mB, mV2, rot180(mKᵥ)); #<! TODO: Remove allocation

    mO .+= mB;

    return mO;

end

function CalcImageLaplacian( mI :: Matrix{T}; mKₕ :: Matrix{T} = [one(T) -one(T)], mKᵥ :: Matrix{T} = [[one(T), -one(T)];;] ) where {T <: AbstractFloat}

    mO = similar(mI);
    mB = similar(mI); #<! Buffer

    numRows = size(mI, 1);
    numCols = size(mI, 2);

    numRowsKₕ = size(mKₕ, 1);
    numColsKₕ = size(mKₕ, 2);
    numRowsmKᵥ = size(mKᵥ, 1);
    numColsmKᵥ = size(mKᵥ, 2);

    mV1 = Matrix{T}(undef, numRows - numRowsKₕ + 1, numCols - numColsKₕ + 1);
    mV2 = Matrix{T}(undef, numRows - numRowsmKᵥ + 1, numCols - numColsmKᵥ + 1);

    mO = CalcImageLaplacian!(mO, mI, mB, mV1, mV2, mKₕ, mKᵥ);

    return mO;

end


function ConvertColorSpace!( mO :: Array{T, 3}, mI :: Array{T, 3}, mC :: Matrix{T} ) where{T <: AbstractFloat}

    numRows = size(mI, 1);
    numCols = size(mI, 2);
    numPx = numRows * numCols;

    mII = @invoke reshape(mI, numPx, 3);
    mOO = @invoke reshape(mO, numPx, 3);

    mul!(mOO, mII, mC');

    return mO;

end

function ConvertColorSpace( mI :: Array{T, 3}, mC :: Matrix{T} ) where{T <: AbstractFloat}

    mO = similar(mI);

    mO = ConvertColorSpace!(mO, mI, mC);

    return mO;

end

function GenColorConversionMat( colorConvMat :: ColorConvMat )
    
    # RGB <-> YCgCo
    # Y in [0, 1], Cg / Co in [-0.5, 0.5]
    mRgbToYCgCo = [+0.2500 +0.5000 +0.2500; -0.2500 +0.5000 -0.2500; +0.5000 00.0000 -0.5000];
    mYCgCoToRgb = [+1.0000 -1.0000 +1.0000; +1.0000 +1.0000 00.0000; +1.0000 -1.0000 -1.0000];
    
    # RGB <-> YPbPr SD TV
    # Y in [0, 1], Pb / Pr in [-0.5, 0.5]
    mRgbToYPbPrSdTv = [+0.2990 +0.5870 +0.1140; -0.1690 -0.3310 +0.5000; +0.5000 -0.4190 -0.0810];
    mYPbPrToRgbSdTv = [+1.0000 -0.0009 +1.4017; +1.0000 -0.3437 -0.7142; +1.0000 +1.7722 +0.0010];
    
    # RGB <-> YPbPr HD TV
    # Y in [0, 1], Pb / Pr in [-0.5, 0.5]
    mRgbToYPbPrHdTv = [+0.2130 +0.7150 +0.0720; -0.1150 -0.3850 +0.5000; +0.5000 -0.4540 -0.0460];
    mYPbPrToRgbHdTv = [+1.0000 +0.0008 +1.5742; +1.0000 -0.1872 -0.4690; +1.0000 +1.8561 +0.0009];
    
    # RGB <-> YUV
    # Y in [0, 1], U in [-0.436, 0.436], V in [-0.615, 0.615]
    mRgbToYuv = [+0.2990 +0.5870 +0.1140; -0.1470 -0.2890 +0.4360; +0.6150 -0.5150 -0.1000];
    mYuvToRgb = [+1.0000 00.0000 +1.1398; +1.0000 -0.3946 -0.5805; +1.0000 +2.0320 -0.0006];

    # RGB <-> YIQ
    # Y in [0, 1], U in [-0.5957, 0.5957], V in [-0.5226, 0.5226]
    mRgbToYiq = [+0.2990 +0.5870 +0.1140; +0.5959 -0.2746 -0.3213; +0.2115 -0.5227 +0.3112];
    mYiqToRgb = [+1.0000 +0.9560 +0.6190; +1.0000 -0.2720 -0.6470; +1.0000 -1.1060 +0.1703];

    if (colorConvMat == RGB_TO_YCGCO)
        return mRgbToYCgCo;
    elseif (colorConvMat == YCGCO_TO_RGB)
        return mYCgCoToRgb;
    elseif (colorConvMat == RGB_TO_YPBPR_SD)
        return mRgbToYPbPrSdTv;
    elseif (colorConvMat == YPBPR_TO_RGB_SD)
        return mYPbPrToRgbSdTv;
    elseif (colorConvMat == RGB_TO_YPBPR_HD)
        return mRgbToYPbPrHdTv;
    elseif (colorConvMat == YPBPR_TO_RGB_HD)
        return mYPbPrToRgbHdTv;
    elseif (colorConvMat == RGB_TO_YUV)
        return mRgbToYuv;
    elseif (colorConvMat == YUV_TO_RGB)
        return mYuvToRgb;
    elseif (colorConvMat == RGB_TO_YIQ)
        return mRgbToYiq;
    elseif (colorConvMat == YIQ_TO_RGB)
        return mYiqToRgb;
    end


end
