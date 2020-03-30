%clear;clc
global Q 
n = 3;
%Q = generateSPDmatrix(n);
x=[0.1 0 0.1]';
maxiter = 1000;
iter = 1;
trustRegionBound = 10;
method = 1;
det(Q)
e = zeros([3,1000]);
while iter < maxiter && sum(abs(g(x))) > 0.00001
    switch method
        case 1
            iter = iter + 1;
            pU = -g(x)' * g(x) * g(x) / (g(x)' * B(x) * g(x));
            pB = -B(x)^-1 * g(x);
            tau = 2;
            if pB'*pB <= trustRegionBound*trustRegionBound
                p = pB;
            end
            if pB'*pB > trustRegionBound*trustRegionBound && pU'*pU >= trustRegionBound*trustRegionBound
                    p = -trustRegionBound * g(x)./norm(g(x),2);
            end
            if pB'*pB > trustRegionBound*trustRegionBound && pU'*pU < trustRegionBound*trustRegionBound
                a = (pB-pU)'*(pB-pU);
                b = 2 * (pB-pU)' * pU;
                c = pU' * pU-trustRegionBound^2;
                tau = 1 + (-b + sqrt(b^2-4*a*c))/(2*a);
            end       
            if tau < 1
                p = tau * pU;
            else
                p = pU + (tau - 1) * (pB - pU);
            end
            e(:,iter) = x(:,1);
        case 2
            iter = iter + 1;
            pU = -g(x)' * g(x) * g(x) / (g(x)' * B(x) * g(x));
            pB = -B(x)^-1 * g(x);
            if g(x)' * B(x) * g(x) <= 0
                p = -trustRegionBound * g(x)./norm(g(x),2);
            else p = -min(1,norm(g(x),2)^3/(trustRegionBound * g(x)' * B(x) *g(x))) * g(x)./norm(g(x),2);
            end
            e(:,iter) = x(:,1);
    end
    %Evaluate rho
    rho = (f(x)-f(x+p))/(Mk(x,zeros(n,1))-Mk(x,p));  %zeros(3,1)
    %update trust region
    if rho > 0.75 && sqrt(abs(p'*p) - trustRegionBound) < 0.001
        trustRegionBound = min(2 * trustRegionBound, 50);
    else if rho < 0.25
            trustRegionBound = sqrt(abs(p'*p)) * 0.25;
        end
    end
    if rho > 0.0001%Mk(zeros(2,1),x) > Mk(p, x)
        x = x + p;
    end
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

function y = Mk(x,p)
y = f(x) + p'*g(x) + 1/2*p'*B(x)*p;
end

function y = f(x)
global Q 
y = log2(1 + x'*Q*x);
end

function y = B(x)
global Q 
y = Q;
end

function y = g(x)
global Q 
y = 2*Q*x./(1+x'*Q*x);
end

function A = generateSPDmatrix(n)
% Generate a dense n x n symmetric, positive definite matrix

A = rand(n,n); % generate a random n x n matrix

% construct a symmetric matrix using either
A = 0.5*(A+A'); 
%A = A*A';

% The first is significantly faster: O(n^2) compared to O(n^3)

% since A(i,j) < 1 by construction and a symmetric diagonally dominant matrix
%   is symmetric positive definite, which can be ensured by adding nI
A = A + n*eye(n);

end
