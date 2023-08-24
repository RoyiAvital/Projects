
clear();
close('all');

numSamples  = 1000;
valF        = 0.5 * rand();

vN = 0:(numSamples - 1);
vN = vN(:);

vY    = rand(numSamples, 1);
hG    = @(valF) G(valF, vN, vY);
hDerG = @(valF) DerG(valF, vN, vY);

numericGrad = CalcFunGrad(valF, hG);
analyticGrad = hDerG(valF);

disp(['Numeric 1st Derivative: ', num2str(numericGrad), ', Analytic 1st Derivative: ', num2str(analyticGrad)]);

disp(['Numeric 2nd Derivative: ', num2str(CalcFunGrad(valF, hDerG)), ', Analytic 2nd Derivative: ', num2str(analyticGrad)]);


function [ outVal ] = G( valF, vN, vY )

mX = [sin(2 * pi * valF * vN), cos(2 * pi * valF * vN)];

vT = ((mX * ((mX.' * mX) \ (mX.' * vY))) - vY);
outVal = 0.5 * sum(vT .^ 2);

end

function [ derVal ] = DerG( valF, vN, vY )

mX = [sin(2 * pi * valF * vN), cos(2 * pi * valF * vN)];

mT0 = inv(mX.' * mX);
mT1 = mX * mT0 * mX.' * vY;
mT2 = vY.' * mX * mT0;
mT3 = mT1 - vY;
mT4 = ((mT2 * mX') - vY.') * mX * mT0;
mG = (mT3 * mT2) - ((mT1 * mT4) + (mX * mT0 * mX.' * mT3 * mT2)) + (vY * mT4);

mT5 = [2 * pi * (vN .* mX(:, 2)), -2 * pi * (vN .* mX(:, 1))];

derVal = trace(mG.' * mT5);


end