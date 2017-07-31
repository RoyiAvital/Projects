function [ vX ] = InterriorPointSolver( vX, hObjFun, hObjGrad, hObjHessian, hCondFun, hCondGrad, hCondHessian, tFctr, muFctr, numIterations )
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

hObjFun = @(vX) (tFctr * hObjFun(vX)) - ((hCondFun(vX) < 0) * log(-hCondFun(vX))) + ((hCondFun(vX) >= 0) * 1e20);

paramAlpha = 0.75;

for ii = 1:numIterations
    
    mCondHessian    = -((hCondGrad(vX) * hCondGrad(vX).') / (hCondFun(vX) .^ 2)) + (hCondHessian(vX) / hCondFun(vX));
    mObjHessian     = tFctr * hObjHessian(vX);
    
    vCondGrad   = hCondGrad(vX) / hCondFun(vX);
    vObjGrad    = tFctr * hObjGrad(vX);
    
    vD       = -((mCondHessian + mObjHessian) \ (vCondGrad + vObjGrad));
    stepSize = LineSearchBackTracking(hObjFun, vX, vD, paramAlpha);
    
    vX = vX + (stepSize * vD);
    
    tFctr = tFctr * muFctr;
end


end

