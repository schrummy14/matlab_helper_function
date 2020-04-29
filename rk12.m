function [t,x,stepSize] = rk12(fun,tspan,x0)

    x(1,:) = x0;
    t(1,1) = tspan(1);
    k = diff(tspan)/1000;
    n = 1;
    done = 0;
    stepSize(1,1) = k;
    
    TOL = 1e-4;
    nsteps = 0;
    nfailed = 0;
    nfevals = 0;
    
    while ~done
        if t(n,1)+k> tspan(2)
            k = tspan(2)-t(n,1);
           done = 1; 
        end
        
        y1 = fun(t(n,1),x(n,:));
        y2 = fun((t(n,1) + k/2),(x(n,:) + k/2*y1'));
        
        nfevals = nfevals + 2;
        
        x1 = x(n,:) + k*y1';
        x2 = x(n,:) + k*y2';
        
        tau = max(abs(x2-x1));
        
        if tau < TOL || k < sqrt(eps)
            x(n+1,:) = x2;
            t(n+1,1) = t(n,1) + k;
            stepSize(n+1,1) = k;
            n = n + 1;
            nsteps = nsteps + 1;
        else
            nfailed = nfailed + 1;
        end
        
        k = 0.95*k*(TOL/(tau+eps))^(1/2);
        
    end
    
    if nargout==1
        Sol.x = t.';
        Sol.y = x.';
        Sol.solver = 'rk12';
        Sol.stats.nsteps = nsteps;
        Sol.stats.nfailed = nfailed;
        Sol.stats.nfevals = nfevals;
        t = Sol;
    end
end