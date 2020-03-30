%{
[Q,R] = qr(rand(1000,1000));
v = randi(990,1000,1)'+ 10*ones([1,1000]);
D = diag(v);
A = Q'*D*Q;
b = randi(10,1000,1);
x = ones([1000,1]);
%}

[Q,R] = qr(rand(1000,1000));
v1 = 0.1*randi(20,700,1)'+ 9*ones([1,700]);
v2 = 0.1*randi(20,300,1)'+ 999*ones([1,300]);
v = [v1 v2];
D = diag(v);
A = Q'*D*Q;
b = randi(10,1000,1);
x = ones([1000,1]);

r = b - A * x;
p = r;
rsold = r' * r;
error = zeros([1,1000]);
    for i = 1:length(b)
        Ap = A * p;
        alpha = rsold / (p' * Ap);
        x = x + alpha * p;
        r = r - alpha * Ap;
        rsnew = r' * r;
        error(i) = log10(x'*A*x);
        if sqrt(rsnew) < 1e-10
              break;
        end
        p = r + (rsnew / rsold) * p;
        rsold = rsnew;
    end
plot(1:1:i,error(1:i),'.'); set(gca,'fontsize', 14);
grid on; hold on;
xlabel('k','fontsize',14); ylabel('Error','fontsize',14)