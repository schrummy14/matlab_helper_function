clear
clc
n = 200;
k = 1000;
Q = zeros(n,n,k);
d = zeros(n,k);
A = Q;
x = d;
b = d;
tic;
parfor i = 1:k
    Q(:,:,i) = orth(randn(n,n));
    d(:,i) = logspace(0,-10,n);
    A(:,:,i) = Q(:,:,i)*diag(d(:,i))*Q(:,:,i)';
    x(:,i) = randn(n,1);
    b(:,i) = A(:,:,i)*x(:,i);
end
toc

errbad = zeros(k,1);
resbad = errbad;
errgood = errbad;
resgood = errbad;

tic;
parfor i = 1:k
    ybad(:,i) = pinv(A(:,:,i))*b(:,i);
    errbad(i) = norm(ybad(:,i)-x(:,i));
    resbad(i) = norm(A(:,:,i)*ybad(:,i)-b(:,i));
    ygood(:,i)= A(:,:,i)\b(:,i);
    errgood(i)= norm(ygood(:,i)-x(:,i));
    resgood(i)= norm(A(:,:,i)*ygood(:,i)-b(:,i));
end
toc

figure(1)
subplot(2,1,1)
hist(errbad,100);
subplot(2,1,2)
hist(resbad,100);

figure(2)
subplot(2,1,1)
hist(errgood,100);
subplot(2,1,2)
hist(resgood,100);
% tic, y = inv(A)*b; toc
% err = norm(y-x)
% res = norm(A*y-b)
% tic, z = A\b; toc
% err = norm(z-x)
% res = norm(A*z-b)