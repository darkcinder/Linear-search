function theta = newton(theta)
tol  = 1e-8; % tolerance: iteration stops when |grad| < tol
maximum_iteration = 2000; % muximum number of iterations
k = 1;
grad = grad_fun(theta);
while k < maximum_iteration && norm(grad,inf)>tol
    H = hessian(theta);
    p = H\(-grad);
    theta = theta+p;
    grad = grad_fun(theta);
    k = k+1;
end
end