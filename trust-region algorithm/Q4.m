clear;clc
global m
m = 1000;
iter = 0;
x=-1*ones([m,1]);

e = zeros([m,1000]);
h = [18 4;4 2];
H = h;

for i = 1:499
    H = blkdiag(H,h);
    i = 1+1;
end

while norm(grad_fun(x,m),inf) > 0.001 && iter < 1000
    iter = iter + 1;
    p = -H * grad_fun(x,m);
    a = StepLength(p, x, 0.1, 0.9, 1e+6,m);
    x = x + a * p;
    s = a * p;
    y = grad_fun(x,m) - grad_fun(x-a * p,m);
    rho = 1/(y' * s);
    I = eye(m,m);
    H = (I - rho * s * y')*H*(I - rho * y * s') + rho * s * s';
    e(:,iter) = x(:,1);   
end

error = zeros(1,iter);
n = 1;
while n < iter
    error(n+1) = norm(e(:,n+1)-e(:,iter-1),inf);
    n = n+1;
end

figure; 
subplot(2,1,1); 
plot(1:1:n,error(1:n),'.'); set(gca,'fontsize', 14);
grid on; hold on;
xlabel('k','fontsize',14); ylabel('Error','fontsize',14)

subplot(2,1,2); 
plot(1:1:n,log10(error(1:n)),'.'); set(gca,'fontsize', 14)
grid on; hold on;
xlabel('k','fontsize',14); ylabel('log10(Error)','fontsize',14)

function a = grad_fun(x,m)
a = zeros([m,1]);
for i = 1:2:(m-1)
    a(i,1) = 2*(x(i+1,1)-x(i,1)^2)*(-2*x(i,1))-2*(1-x(i,1));
end
for j = 2:2:m
    a(j,1) = 2*(x(j,1)-x(j-1,1)^2);
end
end

function b = cost_fun(x,m)
b = 0;
for i = 1:0.5*m
    b = b + 1*(x(2*i,1) - x(2*i-1,1)^2)^2 + (1 - x(2*i-1,1))^2;
end
end

function [ alpha ] = StepLength(p, x0, c1, c2, alpha_max,m)
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

fx0 = cost_fun(x0,m);
gx0 = grad_fun(x0,m);

gx0 = gx0'*p;
fxp = fx0;

i=1;

% Line search algorithm satisfying strong Wolfe conditions
% Algorithms 3.2 on page 59 in Numerical Optimization, by Nocedal and Wright
% alpha_p is alpha_{i-1}
% alpha_x is alpha_i
while true
    xx = x0 + alpha_x*p;
    fxx = cost_fun(xx,m);
    gxx = grad_fun(xx,m);
    
    gxx = gxx'*p;
    if (fxx > fx0 + c1*alpha_x*gx0) || ((i > 1) && (fxx >= fxp))
        alpha = Zoom(x0, p, alpha_p, alpha_x, c1, c2,m);
        return;
    end
    
    if abs(gxx) <= -c2*gx0
        alpha = alpha_x;
        return;
    end
    
    if gxx >= 0
        alpha = Zoom(x0, p, alpha_x, alpha_p, c1, c2,m);
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

function [alpha,fs,gs] = Zoom(x0, p, alpha_low, alpha_high, c1, c2,m)

i =0;
MAX_ITER = 5;

fx0 = cost_fun(x0,m);
gx0 = grad_fun(x0,m)'*p;

while true
    
    alpha_x = 0.5*(alpha_low + alpha_high);
    alpha = alpha_x;
    xx = x0 + alpha_x*p;
    fxx = cost_fun(xx,m);
    gxx = grad_fun(xx,m);
    fs = fxx;
    gs = gxx;
    gxx = gxx'*p;
    x_low = x0 + alpha_low*p;
    f_low = cost_fun(x_low,m);
    
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