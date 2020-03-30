clc
clear all
x = randi([-2 2],100,1);
p = -grad_fun(x);
ff1 = grad_fun(x);
ff2 = grad_fun(x);
k = 1;
e = zeros(100,1000000);

while k < 1000000 && norm(p,inf) > 1e-8
    f1 = grad_fun(x);
    a_ast = steplength(x,p);
    x = x+a_ast*p;
    f2 = grad_fun(x);
    
    beta = f2'*f2/(f1'*f1);    %FR
    
%    if ff2'*ff1/norm(ff2,inf)^2 >= 0.1   %restart beta
 %       beta = 0;
 %   end
    
 %   beta = f2'*(f2-f1)/norm(f1,2)^2;    %PR
    
    p = -f2+beta*p;
    e(:,k) = x(:,1);
    ff1 = f1;
    ff2 = f2;
    k = k+1;
end
error = zeros(1,k);
n = 1;
while n < k
    error(n+1) = norm(e(:,n+1)-e(:,k-1),inf);
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


function grad = grad_fun(x)
grad = zeros(100,1);
for i=2:99
    grad(i,1) = (400*x(i,1)*(x(i,1)^2-x(i+1,1))+2*x(i,1)-2)-200*(x(i-1,1)^2-x(i,1));
end
grad(1,1) = 400*x(1,1)*(x(1,1)^2-x(2,1))+2*x(1,1)-2;
grad(100,1) = -200*(x(99,1)^2-x(100,1));
end

function f = PHI(a,p,x)
f = 0;
for i=1:99
    f = f+100*( (x(i,1)+a*p(i,1))^2-x(i+1,1)-a*p(i+1,1))^2+(x(i,1)+a*p(i,1)-1)^2;
end
end

function g = DPHI(a,p,x)
g = 0;
for i=1:99
    g = g+200*( (x(i,1)+a*p(i,1))^2-x(i+1,1)-a*p(i+1,1))*(2*x(i,1)*p(i,1)+2*a*p(i,1)^2-p(i+1,1))+2*(x(i,1)+a*p(i,1)-1)*p(i,1);
end
end

function a_ast = steplength(x,p)
a0 = 0;
a_ast = 0.5*1e+6; %a_max
phi0 = PHI(0,p,x);
phi_pre = phi0;
dphi0 = DPHI(0,p,x);
i = 1;
while i < 1000
    phi = PHI(a_ast,p,x);
    if phi > phi0 + 0.09*a_ast*dphi0 || (phi >= phi_pre && i >= 2) %parameters in the Strong Wolfe Conditions
        a_ast=zoom(a0,a_ast,x,p,phi0,dphi0);
    end
    dphi = DPHI(a_ast,p,x);
    if abs(dphi) <= -0.4*dphi0 %parameters in the Strong Wolfe Conditions
        break
    end
    if dphi >= 0
        a_ast=zoom(a_ast,a0,x,p,phi0,dphi0);
    end
    a0 = a_ast;
    a_ast =0.5* (a_ast+(0.001));
    phi_pre = phi;
    i = i+1;
end
end

function a = zoom(a_lo,a_hi,x,p,phi0,dphi0)
phi_lo = PHI(a_lo,p,x);
j=1;
while j < 1000
    a = 0.5*(a_lo+a_hi);
    phi = PHI(a,p,x);
    if phi > phi0 + 0.09*a*dphi0 || phi >= phi_lo%parameters in the Strong Wolfe Conditions
        a_hi = a;
    else
        dphi = DPHI(a,p,x);
        if abs(dphi) <= -0.4*dphi0%parameters in the Strong Wolfe Conditions
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
