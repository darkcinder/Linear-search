function grad = grad_fun(theta)
% TODO: Implement the gradient of the objective function.
% defines the gradient of the cost function 
% x: Rn
% grad: Rn -> Rn

global DATA
global LABELS


k = 1;
grad1 = norm(theta,inf);
grad2 = norm(theta,inf);
while k <= 100
    grad1 = grad1-log10(2.7183)*exp(-LABELS(k,1)*theta'*DATA(k, :)')*LABELS(k,1)*DATA(k,1)/(1+exp(-LABELS(k,1)*theta'*DATA(k, :)'));
    grad2 = grad2-log10(2.7183)*exp(-LABELS(k,1)*theta'*DATA(k, :)')*LABELS(k,1)*DATA(k,2)/(1+exp(-LABELS(k,1)*theta'*DATA(k, :)'));
    k = k+1;
end
grad = [grad1;grad2];
end






