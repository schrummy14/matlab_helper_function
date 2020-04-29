clc
clear

M = 10001;
N = 2001;
g = linspace(3.5,4,M);

x = zeros(N,M);
x(1,:) = (g-1)./g + 0.1;
for k = 1:N-1
    x(k+1,:) = g.*x(k,:).*(1-x(k,:));
end

plot(g.',x(end-50:end,:).','b.')