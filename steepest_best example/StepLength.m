function [ alpha ] = StepLength(p, x0, c1, c2, alpha_max)
% Line search step size calculation for strong Wolfe conditions. Requires
% cost_fun(x) and grad_fun(x) to be defined in the same folder as 
% cost_fun.m and grad_fun.m
%
% Input
%   p:          The direction to search
%   x0:         Vector of initial position
%   c1:         Parameter in strong Wolfe conditions
%   c2:         Parameter in strong Wolfe conditions
%   alpha_max:  max step size allowed
%
% Output
%   alpha: step size

MAX_ITER = 3;
alpha_p = 0;
alpha_x = 1;

fx0 = cost_fun(x0);
gx0 = grad_fun(x0);

gx0 = gx0'*p;
fxp = fx0;

i=1;

% Line search algorithm satisfying strong Wolfe conditions
% Algorithms 3.2 on page 59 in Numerical Optimization, by Nocedal and Wright
% alpha_p is alpha_{i-1}
% alpha_x is alpha_i
while true
    xx = x0 + alpha_x*p;
    fxx = cost_fun(xx);
    gxx = grad_fun(xx);
    
    gxx = gxx'*p;
    if (fxx > fx0 + c1*alpha_x*gx0) || ((i > 1) && (fxx >= fxp))
        alpha = Zoom(x0, p, alpha_p, alpha_x, c1, c2);
        return;
    end
    
    if abs(gxx) <= -c2*gx0
        alpha = alpha_x;
        return;
    end
    
    if gxx >= 0
        alpha = Zoom(x0, p, alpha_x, alpha_p, c1, c2);
        return;
    end
    
    alpha_p = alpha_x;
    fxp = fxx;

    if i > MAX_ITER
        alpha = alpha_x;
        return
    end
    % r = rand(1);%randomly choose alpha_x from interval (alpha_p,alpha_max)
    r = 0.8;
    alpha_x = alpha_x + (alpha_max-alpha_x)*r;
    i = i+1;
end
end% end of strongwolfe

function [alpha,fs,gs] = Zoom(x0, p, alpha_low, alpha_high, c1, c2)

i =0;
MAX_ITER = 5;

fx0 = cost_fun(x0);
gx0 = grad_fun(x0)'*p;

while true
    
    alpha_x = 0.5*(alpha_low + alpha_high);
    alpha = alpha_x;
    xx = x0 + alpha_x*p;
    fxx = cost_fun(xx);
    gxx = grad_fun(xx);
    fs = fxx;
    gs = gxx;
    gxx = gxx'*p;
    x_low = x0 + alpha_low*p;
    f_low = cost_fun(x_low);
    
    if ((fxx > fx0 + c1*alpha_x*gx0) || (fxx >= f_low))
        alpha_high = alpha_x;
    else
        if abs(gxx) <= -c2*gx0
            alpha = alpha_x;
            return;
        end
        if gxx*(alpha_high - alpha_low) >= 0
            alpha_high = alpha_low;
        end
        alpha_low = alpha_x;
    end
    
    i = i+1;
    if i > MAX_ITER
        alpha = alpha_x;
        return
    end
end
end