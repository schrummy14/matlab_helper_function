function [x_new,delta,count] = lin_solve_jacobi(A,b)

    x = zeros(size(b));
    x_new = x;
    TOL = 1e-10;
    
    count = 0;
    while true
        for k = 1:length(b)
            sum = 0;
            for m = 1:length(b)
                if k~=m
                    sum = sum + A(k,m)*x(m);
                end
            end
            x_new(k,1) = (b(k) - sum)/A(k,k);
        end
        delta = norm(x-x_new);
        if delta < TOL
            break
        end
        x = x_new;
        count = count + 1;
    end
    
end