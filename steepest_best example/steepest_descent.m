function theta = steepest_descent(theta)

tol  = 1e-8; % tolerance: iteration stops when |grad| < tol
maximum_iteration = 2000; % muximum number of iterations
c1 = 0.1;
c2 = 0.9; % using exact line search
alpha_max = 1e+6;
p = -grad_fun(theta);
k = 1;
while k < maximum_iteration && norm(p,inf)>tol
    alpha = StepLength(p, theta, c1, c2, alpha_max);
    % Definition of line search
    theta = theta + alpha*p;
    % Compute gradient for next step
    p = -grad_fun(theta);     
    k = k+1;
end
end
