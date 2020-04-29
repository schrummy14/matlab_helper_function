function M = SSOR_PreCond(A,w)

    if nargin<2
        w = 1;
    end
    if w == 0
        M = eye(size(A));
    else
        D = diag(A);
        invD = 1./D;
        D = diag(D);
        invD = diag(invD);
        E = tril(A);

        M = w/(2-w)*(1/w * D + E)*invD*(1/w * D + E');
    end
    
end