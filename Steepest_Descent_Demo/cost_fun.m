function cost = cost_fun(x)
% defines cost function we want to minimize
% x: Rn
% cost: Rn -> R

global Q

% A quadratic cost with weighting matrix Q
cost = 0.5*x'*Q*x; 

end






