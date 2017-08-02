function [ vX ] = LeastSquaresInterriorPointSolver( vX, mA, vB, mC, tFctr, muFctr, numIterations )
% ----------------------------------------------------------------------------------------------- %
% [ vX ] = ExamAQ3Solver( vA, paramB )
%  Solving the Bill Payment Problem as in Exam A Question 003.
% Input:
%   - vA            -   Constraints Vector.
%                       The upper bound of the payment of the i-th person.
%                       Must be sorted in ascending order.
%                       Structure: Vector (Column).
%                       Type: 'Single' / 'Double'.
%                       Range: (-inf, inf).
%   - paramB        -   The Bill Sum.
%                       The constriants the sum of the payments must add up
%                       to. The summ of the bill.
%                       Structure: Scalar.
%                       Type: 'Single' / 'Double'.
%                       Range: (0, inf).
% Output:
%   - vX            -   Output Vector.
%                       The payment of each person.
%                       Structure: Vector (Column).
%                       Type: 'Single' / 'Double'.
%                       Range: (-inf, inf).
% References
%   1.  a
% Remarks:
%   1.  a
% TODO:
%   1.  U.
% Release Notes:
%   -   1.0.000     27/07/2017  Royi Avital
%       *   First release version.
% ----------------------------------------------------------------------------------------------- %

FALSE   = 0;
TRUE    = 1;

OFF     = 0;
ON      = 1;

hObjFun = @(vX) (tFctr * sum(((mA * vX) - vB) .^ 2)) - ((all((mC * vX) < 0)) * log(-mC * vX)) + ((any((mC * vX) >= 0)) * 1e20);

paramAlpha = 0.75;

for ii = 1:numIterations
    
    convFlag = FALSE;
    
    while(convFlag == FALSE)
        
        vCondVal = mC * vX;
        
        mCondHessian = (1 / (vCondVal(1) * vCondVal(1))) * (mC(1, :).' * mC(1, :));
        for jj = 2:size(mC, 1)
            mCondHessian = mCondHessian + (1 / (vCondVal(jj) * vCondVal(jj))) * (mC(jj, :).' * mC(jj, :));
        end
        mObjHessian     = 2 * tFctr * (mA.' * mA);
        
        vCondGrad   = -mC.' * (1 ./ vCondVal);
        vObjGrad    = 2 * tFctr * mA.' * ((mA * vX) - vB);
        
        vD       = -((mCondHessian + mObjHessian) \ (vCondGrad + vObjGrad));
        stepSize = LineSearchBackTracking(hObjFun, vX, vD, paramAlpha);
        
        vX = vX + (stepSize * vD);
        
        if(norm(stepSize * vD) < 1e-12)
            convFlag = TRUE;
        end
        
    end
    
    tFctr = tFctr * muFctr;
end


end

