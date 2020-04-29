function [x,fx,k] = funmin_grad_desc(fun,xguess,TOL)


    if nargin < 3
        TOL = 1e-6;
    end

    if nargin < 1
        f = @(x,b)b(1)+b(2).*x+b(3).*x.^2;
        x = linspace(0,1,101);
        y = sin(pi*x);
        fun = @(b)sqrt(sum((f(x,b)-y).^2));
        xguess = [0,1,0].';
    end
    
    MAXITS = 500;
    lam = 0.5;
    grad_f = num_deriv(fun,xguess);
    grad_f_old = grad_f;
    xguess_old = xguess;
    for k = 1:MAXITS
        
        xguess = xguess - lam * grad_f;
        grad_f = num_deriv(fun,xguess);
        lam = (xguess - xguess_old)' * (grad_f - grad_f_old);
        lam = lam/norm(grad_f - grad_f_old)^2;
        lam = min(lam,1);
        
        err_x = norm(xguess - xguess_old);
        err_f = norm(grad_f - grad_f_old);
        
        if max(err_x,err_f) < TOL
            break
        end
        
        xguess_old = xguess;
        grad_f_old = grad_f;
        
        plot(x,y,x,f(x,xguess))
        
    end
    
    x = xguess;
    fx = grad_f;
    
end



function grad_f = num_deriv(fun,x)

    grad_f = zeros(length(x),1);
    f = fun(x);
    x0 = x;
    dx = 1e-6*norm(x);
    for k = 1:length(x)
        x(k) = x(k)+dx;
        grad_f(k) = (fun(x)-f)/dx;
        x = x0;
    end
    
end
    