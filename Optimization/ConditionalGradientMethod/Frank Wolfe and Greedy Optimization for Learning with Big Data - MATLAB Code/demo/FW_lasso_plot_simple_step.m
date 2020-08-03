%
% Demo code running the Frank-Wolfe algorithm with line-search
% for a Lasso problem   min ||Ax-b||^2  s.t. ||x||_1 <= r
%
n=150; % # rows in A
p=20; % dimension of w, # columns of A, variables in x
t=p/4; % number of non-zeros in the true solution xstar
rng('default');
A=randn(n,p);
xstar=[ones(t,1); -ones(t,1); -zeros(p-2*t,1)];
noise=0.1*randn(n,1);
b=A*xstar+noise;
[n,p] = size(A);

T = 15; % number of iterations

r = 1; % the regularization constraint imposed on the l_1-norm


x = zeros(p,1);
Ax = A*x;
progress = [];
figure('Units','normalized','position',[.1 .1 .9 .9])

t=0; i=0;
for t = 0:T
    % call the LMO
    [ res i ] = max( abs(A' *(Ax-b)), [], 1 );
    
	stepSign = sign( A(:,i)' *(b-Ax) ); % = sign of -gradF[i]
    s = zeros(p,1); s(i) = stepSign*r;
    As = stepSign*r*A(:,i); % = the i-th column of the dictionary matrix A
    
    % line-search on the univariate quadratic problem in gamma
    As_minus_Ax = As - Ax;
    stepSize = max(0, min(1, As_minus_Ax' * (b-Ax) / (As_minus_Ax' * As_minus_Ax)));
    %stepSize = 2 / (t+2);
    
    x = (1-stepSize)*x + stepSize*s; % do the FW step
    Ax = (1-stepSize)*Ax + stepSize*As; % update  Ax
    
    f_at_x = sum((Ax-b).^2)/2;
    fprintf('t=%d - f=%f - i=%d - stepsize=%d \n',t,f_at_x,i,stepSize);
    
    progress = [ progress; f_at_x ];
    show_progress(progress,x,s,p,T,r);
    pause(0.2);
    
end

