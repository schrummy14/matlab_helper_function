% Find eigen values of a matrix
function [w1,u1,its] = findEigenValue(A)
    
    u0 = ones(length(A),1);
    u1 = A*u0;
    err = 10;
    w1 = 1;
    its = 0;
    while err>eps||its<5
        its = its + 1;
        w0 = w1;
        u0 = u1;
        u1 = A*u0;
        w1 = u1(1)/u0(1);
        u1 = u1/(u1(1));
        err = norm(u1-u0)/norm(w0,inf);
    end
end
    