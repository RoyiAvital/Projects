%
% Demo code running the Frank-Wolfe algorithm with line-search
% for a nuclear norm constrained problem
% (squared loss matrix completion)
%     min ||X-Y||_Obs^2  s.t. ||X||_* <= r
%
% MovieLens100k rating prediction experiment 
load('data/movielens100k/randsplit1','data'); 
% not provided in this code due to license. download directly from
%   http://grouplens.org/datasets/movielens/
Y = data; clear data;

% prepare objective function and gradient
objectiveF = @(X) matrix_completion_squared_error(X,Y);
gradientF =  @(X) matrix_completion_squared_error_gradient(X,Y);
[m,n] = size(Y);

r = 1000; % nuclear norm constraint on X
T = 200;

%% FW algorithm
X = zeros(m,n); % starting point

for t = 0:T
    [u,s,v] = svds(-gradientF(X),1); % LMO is solved by top singular vector pair
    stepSize = 2/(t+2);              % or perform line-search
    X = (1-stepSize)*X + stepSize* r * u*v'; 
end

% NOTE: for larger matrices, don't store the dense iterate X, but just its
% restriction to the train and test entries, or its factorized (rank t)
% representation
fprintf('X after %g iter of FW has error %g and nuclear norm %g \n',T,objectiveF(X),sum(svd(X)));