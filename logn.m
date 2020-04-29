function y = logn(x,n)
%% 
% y = logn(x,n)
% Produces the log of x with respect to base n.
% examples
% y = logn(1000,10) = log10(1000) = 3
% y = logn(exp(4),exp(1)) = log(exp(4)) = 4

y = log10(x)./log10(n);
