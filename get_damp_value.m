function nlm = get_damp_value(t,x,do_normalize)

    if nargin < 3
        do_normalize = false;
    end

    fun = @(b,t)(b(1) + b(2).*exp(b(3).*t).*cos(b(4).*t));
        
    [nt,nx,c,d] = normalize_freq(t,x);

    b0 = fminunc(@(b)(norm(nx-(fun(b,nt)))),[mean(nx) 1 -1 10]);
    
    if do_normalize
        t = nt;
        x = nx;
    else
        b0(1) = b0(1)*c;
        b0(2) = b0(2)*c;
        b0(3) = b0(3)/abs(d);
        b0(4) = b0(4)/abs(d);
        t = nt*c;
        x = nx*d;
    end
    
    nlm = fitnlm(t,x,fun,b0);

end