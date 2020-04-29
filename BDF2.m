function [t,x] = BDF2(fun,Tspan,x0,N,Jac,Tol,MaxIt)
    
    k = diff(Tspan)/N;
    x = zeros(N+1,length(x0));
    t = zeros(N+1,1);
    x(1,:) = x0;
    t(1,1) = Tspan(1);
    for n = 1:N
        mstop = 0;              
        j = 0;                  
        xguess = x(n,:);        
        if n == 1
            while ~mstop && j < MaxIt   
                j = j + 1;
                f = x(n,:)' + k*fun(t(n,1)+k,xguess) - xguess';
                g = k*Jac(t(n,1)+k,xguess)-eye(length(xguess));
                delta = g\f;
                if max(abs(delta)) <= Tol
                    mstop = 1;
                    x(n+1,:) = xguess - delta';
                    t(n+1,1) = Tspan(1) + n*k;
                else
                    xguess = xguess - delta';
                end
            end
        else
            while ~mstop && j < MaxIt
                j = j + 1;
                f = 2*k*fun(t(n,1)+k,xguess)-x(n-1,:)'+4*x(n,:)'-3*xguess';
                g = 2*k*Jac(t(n,1)+k,xguess)-3*eye(length(xguess));
                delta = g\f;
                if max(abs(delta)) <= Tol
                    mstop = 1;
                    x(n+1,:) = xguess - delta';
                    t(n+1,1) = Tspan(1) + n*k;
                else
                    xguess = xguess - delta';
                end
            end
        end
        if mstop == 0,error('Did not find an answer in time...');end
    end
end