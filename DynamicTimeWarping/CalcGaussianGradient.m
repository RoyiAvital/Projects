function [ vGaussianGradXKernel ] = CalcGaussianGradient( gaussianKernelStd, stdToRadiusFactor )

gaussianKernelRadius  = ceil(stdToRadiusFactor * gaussianKernelStd); % Imitating Photoshop - See Reference

vGaussianKernelSupport  = [-gaussianKernelRadius:gaussianKernelRadius];
vGaussianKernel         = exp(-(vGaussianKernelSupport .* vGaussianKernelSupport) / (2 * gaussianKernelStd * gaussianKernelStd));
vGaussianKernel         = vGaussianKernel / sum(vGaussianKernel);
vGaussianGradXKernel    = -vGaussianKernelSupport .* vGaussianKernel;
% vGaussianGradXKernel    = (vGaussianGradXKernel / sum(abs(vGaussianGradXKernel)));

vNegValsIdx = (vGaussianGradXKernel < 0);
vPosValsIdx = (vGaussianGradXKernel > 0);
vGaussianGradXKernel(vPosValsIdx) = vGaussianGradXKernel(vPosValsIdx)/sum(vGaussianGradXKernel(vPosValsIdx));
vGaussianGradXKernel(vNegValsIdx) = vGaussianGradXKernel(vNegValsIdx)/abs(sum(vGaussianGradXKernel(vNegValsIdx)));


end

