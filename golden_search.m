function [x, its] = golden_search(f,a,b,tol)

    if nargin == 0
        [x, its] = golden_search(@(x)(x-2)^2,1,5,1e-10);
        return
    end

    if nargin < 4
        tol = sqrt(eps);
        if nargin < 3
            error('Must provide function and two starting values')
        end
    end
    
    gr = 0.5*(sqrt(5)+1);
    
    c = b - (b-a)/gr;
    d = a + (b-a)/gr;
    
    its = 0;
    err = abs(c-d);
    
    its_plot = its;
    err_plot = abs(2-0.5*(a+b));
    while err > tol
        its = its + 1;
        if f(c) < f(d)
            b = d;
        else
            a = c;
        end
        c = b - (b-a)/gr;
        d = a + (b-a)/gr;
        err = abs(c-d);
        
        its_plot(end+1,1) = its;
        err_plot(end+1,1) = abs(2-0.5*(a+b));
    end
    
    x = 0.5*(a+b);
    
    semilogy(its_plot,err_plot)
    
end