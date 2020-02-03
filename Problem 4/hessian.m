function [ H ] = hessian( theta )
% TODO: Implement the hessian of the cost function.

global DATA;
global LABELS;

H = eye(2);
h1 = 1;
h2 = 1;
h3 = 1;
h4 = 1;
k = 1;
while k <= 100
    h1 = h1-log10(2.7183)*LABELS(k,1)*DATA(k,1)*(-exp(LABELS(k,1)*theta'*DATA(k, :)'))*LABELS(k,1)*DATA(k,1)/(1+exp(LABELS(k,1)*theta'*DATA(k, :)'))^2;
    h2 = h2-log10(2.7183)*LABELS(k,1)*DATA(k,1)*(-exp(LABELS(k,1)*theta'*DATA(k, :)'))*LABELS(k,1)*DATA(k,2)/(1+exp(LABELS(k,1)*theta'*DATA(k, :)'))^2;
    h3 = h3-log10(2.7183)*LABELS(k,1)*DATA(k,2)*(-exp(LABELS(k,1)*theta'*DATA(k, :)'))*LABELS(k,1)*DATA(k,1)/(1+exp(LABELS(k,1)*theta'*DATA(k, :)'))^2;
    h4 = h4-log10(2.7183)*LABELS(k,1)*DATA(k,2)*(-exp(LABELS(k,1)*theta'*DATA(k, :)'))*LABELS(k,1)*DATA(k,2)/(1+exp(LABELS(k,1)*theta'*DATA(k, :)'))^2;
    k = k+1;
end
H = [h1,h2;h3,h4];
end

