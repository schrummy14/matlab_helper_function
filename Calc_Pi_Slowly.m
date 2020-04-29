% clear
% clc

TOL = 1e-8;
err = 1;
k = 0;
tot = 0;
method = 2;
tic;

switch method
    case 1
        while err>TOL
            tot = tot + (-1)^k / (2*k+1);
            err = abs(4*tot - pi);
            k = k + 1;
        end
        tot = 4*tot;
    case 2
        while err>TOL
            k = k + 1;
            K = 2*k;
            tot = tot + (-1)^(k-1)/(K*(K+1)*(K+2));
            err = abs(3+4*tot-pi);
        end
        tot = 4*tot + 3;
end
toc;