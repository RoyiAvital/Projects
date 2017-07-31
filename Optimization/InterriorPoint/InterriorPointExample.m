% Interrior Point Example
% Remarks:
%   1.  Using Log Barrier Method for Interrior Point Solver
% TODO:
% 	1.  ds
% Release Notes
% - 1.0.000     31/07/2017  Royi Avital
%   *   First release.


%% General Parameters

run('InitScript.m');


%% Question 003

hObjFun     = @(vX) (2 * (vX .^ 2)) + (3 * vX);
hObjGrad    = @(vX) (4 * vX) + 3;
hObjHessian = @(vX) 4;

% hCondFun        = @(vX) -sqrt(abs(vX)) + ((vX < 0) * 1e20);
% hCondGrad       = @(vX) -0.5 * (vX .^ -0.5);
% hCondHessian    = @(vX) 0.25 * (vX .^ -1.5);

hCondFun        = @(vX) ((vX <= 0) * vX) + ((vX > 0) * 1e20);
hCondGrad       = @(vX) 1;
hCondHessian    = @(vX) 0;


%% Solution by CVX - Exact

cvx_begin('quiet')
    cvx_precision('best');
    variable vX(1)
    minimize( (2 * sum_square_abs(vX)) + (3 * vX) )
    subject to
        -sqrt(vX) <= 0;
        % vX <= 0;
cvx_end

disp([' ']);
disp(['CVX Solution Summary']);
disp(['The CVX Solver Status - ', cvx_status]);
disp(['The Optimal Value Is Given By - ', num2str(cvx_optval)]);
disp(['The Optimal Argument Is Given By - [ ', num2str(vX.'), ' ]']);
disp([' ']);


%% Solution by CVX - Approximate

vXInit          = 5;
tFctr           = 8880;
muFctr          = 5;
numIterations   = 55;

vX = InterriorPointSolver(vXInit, hObjFun, hObjGrad, hObjHessian, hCondFun, hCondGrad, hCondHessian, tFctr, muFctr, numIterations);

vX


%% Restore Defaults

% set(0, 'DefaultFigureWindowStyle', 'normal');
% set(0, 'DefaultAxesLooseInset', defaultLoosInset);

