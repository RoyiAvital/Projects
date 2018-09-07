function [ mO ] = BilateralGaussianBlurCompressiveFast( mI, spatialStd, rangeStd, paramK )
% ----------------------------------------------------------------------------------------------- %
% [ mO ] = BilateralGaussianBlurCompressiveFast( mI, spatialStd, rangeStd, paramK )
%   Applies Bilateral Gaussian Blur on an image using Compressive Bilateral
%   Filter approach to have high quality approximation with lower
%   complexity. The Compressive approach uses a Fourier Series to
%   approximate the Range Filter and spearate it from the Spatial Filter.
% Input:
%   - mI                -   Input Image.
%                           Structure: Image Matrix (Single Channel)
%                           Type: 'Single' / 'Double'.
%                           Range: [0, 1].
%   - spatialStd        -   Spatial Standard Deviation.
%                           The STD of Gaussian Kernel for the spatial weights.
%                           Structure: Scalar.
%                           Type: 'Single' / 'Double'.
%                           Range: (0, inf).
%   - rangeStd          -   Range Standard Deviation.
%                           The STD of Gaussian Kernel for the range weights.
%                           Structure: Scalar.
%                           Type: 'Single' / 'Double'.
%                           Range: (0, inf).
%   - paramK            -   Parameter K.
%                           The coefficient to convert from STD to radius f
%                           the filters. The higher the value the better
%                           the quality of the output. Usually value of 3
%                           is good enough. High values means higher
%                           complexity.
%                           Structure: Scalar.
%                           Type: 'Single' / 'Double'.
%                           Range: {3, 4, 5, ...}.
% Output:
%   - mO                -   Output Image.
%                           Structure: Image Matrix (Single Channel)
%                           Type: 'Single' / 'Double'.
%                           Range: [0, 1].
% References:
%   1.  Compressive Bilateral Filtering (https://ieeexplore.ieee.org/document/7120121/).
%   2.  Fast Compressive Bilateral Filter (https://ieeexplore.ieee.org/document/7843844/).
% Remarks:
%   1.  Quality is mainly determined by 'paramL' (Think of it as a
%       generalization of the Range Filter Radius). So either get it higher
%       with 'paramK' that for a given 'rangeStd' will add more
%       coefficients or by higher value of 'rangeStd' that makes the filter
%       easier to approximate (Closer to Spatial Invariant Gaussian Blur).
% TODO:
%   1.  Optimize allocation of arrays.
%   Release Notes:
%   -   1.0.001     06/09/2018  Royi Avital
%       *   Using direct call to 'ApplyGaussianBlur()'.
%       *   Optimized generation of the paramters ('paramL' and 'paramN').
%   -   1.0.000     05/09/2018  Royi Avital
%       *   First release version.
% ----------------------------------------------------------------------------------------------- %

FALSE   = 0;
TRUE    = 1;

OFF = 0;
ON  = 1;

% Setting Parameters
dataType    = class(mI);
numRows     = size(mI, 1);
numCols     = size(mI, 2);

% Parameter L -> Filter Radius
% I adapted this to match images in [0, 1] range.
% Also removing the requirement of 'paramL' to be integer releases 'paramN' to be related to rangeSTD.
% It counter logic to some degree, it also generates good results yet the
% maximum absolute deviation gets much bigger. It seems that the
% approxmation is better as paramL gets higher. So quality is constant for
% given paramL, either raise it by higher value of 'rangeSTD' (Intuitive,
% easier to approxmate something which closer to Gaussian Blur) or by
% higher value of 'paramK'. Maybe add "High Quality Mode" with 'paramL =
% ceil(paramK * rangeStd)' as it seems it yields good enough results.
paramL      = paramK * rangeStd; %<! L in eq. 6
paramTau    = paramK / pi; %<! \tau in the paper
paramN      = ceil(paramK * paramK / pi); %<! Number of Fourier Series Coef (Eq. 7)
paramOmega  = pi / paramL;
vParamD     = 2 * exp(-((1:paramN) .^ 2) / (2 * paramTau * paramTau)); %<! Gamma_n in the paper

% Pre Allocation
mO = zeros(numRows, numCols, dataType); 
mZ = mO;

% Filtering Process - Iteraion number 1
mParamOmega = paramOmega * mI;
mCOmega     = cos(mParamOmega);
mSOmega     = sin(mParamOmega);

mCFiltered  = ApplyGaussianBlur(mCOmega, spatialStd, paramK);
mSFiltered  = ApplyGaussianBlur(mSOmega, spatialStd, paramK);
mO          = mO + vParamD(1) * (mCOmega .* mSFiltered - mSOmega .* mCFiltered);
mZ          = mZ + vParamD(1) * (mCOmega .* mCFiltered + mSOmega .* mSFiltered);

% Iteration Number 2
mS = 2 * mCOmega .* mSOmega; %<! sin(2 * omega)
mC = 2 * mCOmega .* mCOmega - 1; %<! cos(2 * omega);
mCFiltered = ApplyGaussianBlur(mC, spatialStd, paramK);
mSFiltered = ApplyGaussianBlur(mS, spatialStd, paramK);
mO = mO + (2 * vParamD(2)) * (mC .* mSFiltered - mS .* mCFiltered);
mZ = mZ + vParamD(2) *(mC .* mCFiltered + mS .* mSFiltered);

for ii = 3:paramN
    mT = mC .* mSOmega + mS .* mCOmega; %<! sin(ii * omega) [Buffer]
    mC = mC .* mCOmega - mS .* mSOmega; %<! Eq. 17;
    mS = mT;
    mCFiltered = ApplyGaussianBlur(mC, spatialStd, paramK);
    mSFiltered = ApplyGaussianBlur(mS, spatialStd, paramK);
    mO = mO + (ii * vParamD(ii)) * (mC .* mSFiltered - mS .* mCFiltered);
    mZ = mZ + vParamD(ii) * (mC .* mCFiltered + mS .* mSFiltered);
end

% Sum of filter coefficients
% mA = hGaussianBlur(ones(size(mI), 'like', mI), spatialStd); %<! Filter is normalized to 1
% mO = mO ./ (mA + mZ); %<! Filter is normalized to 1
mO = mO ./ (1 + mZ);

mO = mI + (pi * (rangeStd * rangeStd) / paramL) * mO;


end

