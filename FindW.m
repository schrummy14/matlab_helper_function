function tau = FindW(w,A,b)
    [~,~,tau] = Jacobi(A,b,'weight',w);
end
    