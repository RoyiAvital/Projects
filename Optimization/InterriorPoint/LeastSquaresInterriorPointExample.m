% Least Squares Interrior Point Example
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

numRows = 4;
numCols = 2;

numConstraints = 4;

mA = randn([numRows, numCols]);
vB = randn([numRows, 1]);

mC = abs(randn([numConstraints, numCols])); %<! For easy initialization



%% Solution by CVX - Exact

cvx_begin('quiet')
    cvx_precision('best');
    variable vX(numCols)
    minimize( sum_square_abs(mA * vX - vB) )
    subject to
        mC * vX <= 0
cvx_end

disp([' ']);
disp(['CVX Solution Summary']);
disp(['The CVX Solver Status - ', cvx_status]);
disp(['The Optimal Value Is Given By - ', num2str(cvx_optval)]);
disp(['The Optimal Argument Is Given By - [ ', num2str(vX.'), ' ]']);
disp([' ']);


%% Solution by CVX - Approximate

vXInit          = -5 * ones([numCols, 1]);
tFctr           = 8;
muFctr          = 5;
numIterations   = 35;

vX = LeastSquaresInterriorPointSolver(vXInit, mA, vB, mC, tFctr, muFctr, numIterations);

vX


%% Restore Defaults

% set(0, 'DefaultFigureWindowStyle', 'normal');
% set(0, 'DefaultAxesLooseInset', defaultLoosInset);

