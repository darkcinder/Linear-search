
% steepest descent method with exact line search 
% for a quadratic and convex cost function 
% f(x) = 0.5*x'*Q*x - b'*x
% where x\in R^n, Q is symmetric and P.D.
clc
clear all

global Q

% define line search parameters for Wolfe conditions
c1 = 0.1;
c2 = 0.9; % using exact line search
a_max = 1e+6;

tol  = 1e-8; % tolerance: iteration stops when |grad| < tol
maximum_iteration = 2000; % muximum number of iterations

n = 100; % dimension of decision variable x
x = zeros(n,maximum_iteration); % x stores approximated solutions at each iterations 
error = zeros(1,maximum_iteration); % value of cost function at each iteration (for convergence analysis)

Q = rand(n,n);
[Q,R] = qr(Q); %QR factorization of Q
lambda_min = 1;
lambda_max = 1000;
lambda = linspace(lambda_min,lambda_max,n); % eigenvalue ranges from 1 to 1000
Q = Q'*diag(lambda)*Q; % generate a PD and symmetric matrix Q with condition number = lambda_max/lambda_min

x_ast = zeros(n,1); % global min

% initial step
x(:,1) = 10*(2*rand(n,1)-1);
error(1) = sqrt((x(:,1)-x_ast)'*Q*(x(:,1)-x_ast));

% search direction is -grad f for steepest descent
p = -grad_fun(x(:,1));

k = 1;
while k < maximum_iteration && norm(p,inf)>tol
    
    % compute the step length along search direction p
    a = (p'*p)/(p'*Q*p); % exact line search formula derived in class
    
%     a = StepLength(p,x(:,k),c1,c2,a_max); % inexact line search
    
    % Definition of line search
    x(:,k+1) = x(:,k) + a*p;
    
    % Compute gradient for next step
    p = -grad_fun(x(:,k+1));
    
    % Q-norm error
    error(k+1) = sqrt((x(:,k+1)-x_ast)'*Q*(x(:,k+1)-x_ast));
    
    k = k+1;
end

% convergence rate: error measured in Q-norm
figure; 
subplot(2,1,1); 
plot(1:1:k,error(1:k),'.'); set(gca,'fontsize', 14);
grid on; hold on;
xlabel('k','fontsize',14); ylabel('Error','fontsize',14)

subplot(2,1,2); 
plot(1:1:k,log10(error(1:k)),'.'); set(gca,'fontsize', 14)
grid on; hold on;
xlabel('k','fontsize',14); ylabel('log10(Error)','fontsize',14)

K = cond(Q);
r = (K-1)/(K+1);
slope = (log10(error(2000))-log10(error(1000)))/(2000-1000)
slope_est = log10(r)

