function [x,res,iters] = gauss_seidel(A, b, TOL, x)

    if nargin < 4
        x = zeros(size(b));
        if nargin < 3
            TOL = 1e-10;
        end
    end
    n = length(A);
    iters = 0;
    err1 = 0;
    MaxIts = 5000;
    done = false;
    ff = 0;
    Sig_Change = zeros(length(x),2);
    while ~done
        ff = ff+1;
        iters = iters+1;
        err2 = err1;
        for m = 1:n
            sig = 0;
            for k =1:n
                if(m~=k)
                    sig = sig + A(m,k)*x(k,:);
                end
            end
            Sig_Change(m,1) = sig;
            x(m,:) = (b(m,:)-sig)/A(m,m);
        end
        res = norm(A*x-b);
        err1 = norm(Sig_Change(:,1)-Sig_Change(:,2),inf);
        Sig_Change(:,2) = Sig_Change(:,1);
        errNorm = norm(err1-err2);
        if res <= TOL || iters>MaxIts
%         if err1 < TOL || errNorm <=TOL ||iters>MaxIts
            done = true;
        end
    end
    res = err1;
end