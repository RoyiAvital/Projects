function g = matrix_completion_squared_error_nabla(X,Y)
% Matrix completion squared error objective:
%  min_X (1/2)(X(Y~=0) - Y(Y~=0))^2

[rows,cols,vals] = find(Y);

% Compute objective
preds = X(Y~=0);
M = preds - vals;
%f = sum(M(:).^2)/2;

% Compute gradient
g = sparse(rows,cols,M);