
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

% disp(['Numeric 2nd Derivative: ', num2str(CalcFunGrad(valF, hDerG)), ', Analytic 2nd Derivative: ', num2str(analyticGrad)]);


function [ outVal ] = G( valF, vN, vY )

mX = [sin(2 * pi * valF * vN), cos(2 * pi * valF * vN)];

mZ = mX.' * mX; %<! Keep complex (Not conjugate) for complex step
% mL = chol(inv(mZ), 'lower'); %<! mL * mL.' = mZ; %<! Won't work with complex step
mL = sqrtm(inv(mZ)); %<! mL * mL = mZ;

vT = vY.' * (mX * mL);
outVal = sum(vT .^ 2);

% outVal = (vY' * mX) * inv(mX' * mX) * (mX' * vY);

end

function [ derVal ] = DerG( valF, vN, vY )

mX = [sin(2 * pi * valF * vN), cos(2 * pi * valF * vN)];

mT0 = inv(mX.' * mX);
mT1 = vY.' * (mX * mT0);
mT2 = vY * mT1;
mG  = 2 * (mT2 - (mX * (mT0 * (mX' * mT2))));
mT5 = [2 * pi * (vN .* mX(:, 2)), -2 * pi * (vN .* mX(:, 1))];

derVal = trace(mG.' * mT5);


end