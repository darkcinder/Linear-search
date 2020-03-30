function z = cost_fun(theta)
% TODO: Implement the objective function.
% defines cost function 
% x: Rn
% cost: Rn -> R

global DATA
global LABELS

k = 1;
z = 0.5*norm(theta,inf)^2;
while k <= 100
    z = z-log10(1/(1 + exp(-LABELS(k,1)*theta'*DATA(k, :)')));
    k = k+1;
end
end






