clc
clear all
syms x1 x2;
x1 = 0;
x2 = -1;
c = 10;
k = 1;
p = -grad_fun(c,x1,x2);
f = zeros(1,2000);
while k < 2000 && norm(p,inf) > 1e-8
    p = -grad_fun(c,x1,x2);
    a_ast = steplength(c,x1,x2,p);
    x1 = x1+a_ast*p(1);
    x2 = x2+a_ast*p(2);
    f(k) = (c*(x1+a_ast*p(1))-2)^4+((x2+a_ast*p(2))^2)*(c*(x1+a_ast*p(1))-2)^2+(x2+a_ast*p(2)+1)^2;
    k = k+1;
end
n = 1;
while n < k
    r(n) = log((f(n+1)-f(k))/(f(n)-f(k)));
    n = n+1;
end
plot(r)

function grad = grad_fun(c,x1,x2)
grad = [4*c*(c*x1-2)^3+2*c*(x2^2)*(c*x1-2),2*x2*(c*x1-2)^2+2*(x2+1)];
end

function a_ast = steplength(c,x1,x2,p)
a0 = 0;
a_ast = 0.5*1e+6; %a_max
phi0 = (c*x1-2)^4+(x2^2)*(c*x1-2)^2+x2^2;
phi_pre = phi0;
dphi0 = 4*c*p(1)*(c*(x1)-2)^3+2*(x2)*p(2)*(c*(x1)-2)^2+((x2)^2)*2*(c*(x1)-2)*c*p(1)+2*(x2)*p(2);
i = 1;
while i < 1000
    phi = (c*(x1+a_ast*p(1))-2)^4+((x2+a_ast*p(2))^2)*(c*(x1+a_ast*p(1))-2)^2+(x2+a_ast*p(2)+1)^2;
    if phi > phi0 + 0.1*a_ast*dphi0 || (phi >= phi_pre && i >= 2) %parameters in the Strong Wolfe Conditions
        a_ast=zoom(a0,a_ast,x1,x2,p,phi0,dphi0,c);
    end
    dphi = 4*c*p(1)*(c*(x1+a_ast*p(1))-2)^3+2*(x2+a_ast*p(2))*p(2)*(c*(x1+a_ast*p(1))-2)^2+((x2+a_ast*p(2))^2)*2*(c*(x1+a_ast*p(1))-2)*c*p(1)+2*(x2+a_ast*p(2)+1)*p(2);
    if abs(dphi) <= -0.9*dphi0 %parameters in the Strong Wolfe Conditions
        break
    end
    if dphi >= 0
        a_ast=zoom(a_ast,a0,x1,x2,p,phi0,dphi0,c);
    end
    a0 = a_ast;
    a_ast =0.5* (a_ast+(0.001));
    phi_pre = phi;
    i = i+1;
end
end

function a = zoom(a_lo,a_hi,x1,x2,p,phi0,dphi0,c)
phi_lo = (c*(x1+a_lo*p(1))-2)^4+((x2+a_lo*p(2))^2)*(c*(x1+a_lo*p(1))-2)^2+(x2+a_lo*p(2)+1)^2;
j=1;
while j < 1000
    a = 0.5*(a_lo+a_hi);
    phi = (c*(x1+a*p(1))-2)^4+((x2+a*p(2))^2)*(c*(x1+a*p(1))-2)^2+(x2+a*p(2)+1)^2;
    if phi > phi0 + 0.1*a*dphi0 || phi >= phi_lo
        a_hi = a;
    else
        dphi = 4*c*p(1)*(c*(x1+a*p(1))-2)^3+2*(x2+a*p(2))*p(2)*(c*(x1+a*p(1))-2)^2+((x2+a*p(2))^2)*2*(c*(x1+a*p(1))-2)*c*p(1)+2*(x2+a*p(2)+1)*p(2);
        if abs(dphi) <= -0.9*dphi0
            break
        end
        if dphi*(a_hi-a_lo) >= 0
            a_hi = a_lo;
        end
        a_lo = a;
        phi_lo = phi;
    end
    j = j+1;
end
end
