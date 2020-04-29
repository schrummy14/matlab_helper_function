function varargout = TR_BDF2(fun,tspan,x0,options)
    
    if nargin<4
        options = [];        
        if nargin == 0 % Test of function
            c = -10e10;
            options = {...
                'refine',1,...
                'Rel_TOL',1e-3,...
                'Abs_TOL',1e-6,...
                'Do_Plots',true,...
                'Do_Waitbar',true,...
                'dfun',@(t,x)c};
            fun = @(t,x)(c*(x-cos(t))-sin(t));
            tspan = [0 2*pi];
            x0 = 3;
        end
    end
    
    [mu,min_dt,max_dt,Newton_TOL,Rel_TOL,Abs_TOL,dfun,refine_n,do_plots,...
        do_waitbar] = ...
        getOptions(options);
    
    num_eqs = length(x0);
    k_mu = (-3*mu^2 + 4*mu - 2)/(12*(2-mu));
    pow = 1/(2+1);
    nEval = 0;
    nFail = 0;
    waitbar_k = 0;
    
    x_out(1,:) = x0;
    t_out(1,1) = tspan(1);
    x = x0(:);
    t = tspan(1);
    dt = 16*(min_dt);
    done = false;
    Max_Newton_Its = 500;
    N_Tries = 0;
    while ~done
        
        min_dt = max(eps(t),min_dt);
        
        N_Tries = N_Tries + 1;
        if t+dt>tspan(end)
            dt = tspan(end)-t;
            done = true;
        end
        
        f = feval(fun,t,x);
        nEval = nEval + 1;
        delta = ones(num_eqs,1);
        x_guess = x;
        a1 = mu*dt/2;
        
        num_its = 0;
        N_Failed = false;
        while norm(delta,inf)>Newton_TOL && num_its < Max_Newton_Its
            num_its = num_its + 1;
            f_mu = feval(fun,t+mu*dt,x_guess);
            [df,jac_evals] = jac(fun,dfun,t+mu*dt,x_guess);
            delta = (a1*df-eye(num_eqs))\(x + a1*(f+f_mu)-x_guess);
            x_guess = x_guess - delta;
            
            nEval = nEval + 1 + jac_evals;
            
        end
        if num_its >= Max_Newton_Its
            N_Failed = true;
        end
        
        delta = ones(num_eqs,1);
        a1 = (1-mu)*dt/(2-mu);
        a2 = 1/(mu*(2-mu));
        a3 = a2*(1-mu)^2;
        
        f_mu = feval(fun,t+mu*dt,x_guess);
        nEval = nEval + 1;
        x_mu = x_guess;
        
        num_its = 0;
        while norm(delta,inf)>Newton_TOL && num_its < Max_Newton_Its && ~N_Failed
            num_its = num_its + 1;
            f_n = feval(fun,t+dt,x_guess);
            [df,jac_evals] = jac(fun,dfun,t+dt,x_guess);
            delta = (a1*df-eye(num_eqs))\(a1*f_n + a2*x_mu - a3*x - x_guess);
            x_guess = x_guess - delta;
            
            nEval = nEval + 1 + jac_evals;
            
        end
        if num_its >= Max_Newton_Its
            N_Failed = true;
        end        
        
        f_n = feval(fun,t+dt,x_guess);
        f_error = 2*k_mu*dt*(f/mu - f_mu/(mu*(1-mu)) + f_n/(1-mu));
        nEval = nEval + 1;
        if N_Failed 
            r = 10;
        else
            r = norm(f_error,inf)/(Rel_TOL*norm(x_guess,inf)+Abs_TOL);
        end
        
        if r<=2
            [t_new,x_new] = refine_grid(t,x,t+dt,x_guess,f,f_n,refine_n);
            x_out((end+1):(end+refine_n+1),:) = x_new;
            t_out((end+1):(end+refine_n+1),1) = t_new;
            
            if(do_plots)
                plot(t_out,x_out,'-o');
                xlim([tspan(1),tspan(end)])
                drawnow;
            end
            
            if(do_waitbar)
                if waitbar_k==0
                    wait_vals = linspace(tspan(1),tspan(end),500);
                    wait_h = waitbar(0,'Solving System...');
                    waitbar_k = 1;
                end
                if t+dt>wait_vals(waitbar_k)
                    waitbar_k = find(wait_vals>t+dt,1);
                    if isempty(waitbar_k)
                        waitbar(1,wait_h);
                    else
                        waitbar(wait_vals(waitbar_k-1)/tspan(end),wait_h);
                    end
                end
            end                
            
            t = t + dt;
            x = x_guess;
            N_Tries = 0;
        else
            done = false;
            nFail = nFail + 1;
        end
        
        if r < 0
            dt = 2*dt;
        else
            dt = 0.95*dt/(r^pow+eps);
        end
        
        dt = min(max_dt,dt);
        dt = max(min_dt,dt);
        
        if N_Tries > 500
            break
        end
        
    end
    
    if nargout > 1
        varargout{1} = t_out;
        varargout{2} = x_out;
    elseif nargout > 0
        SOL.t = t_out;
        SOL.x = x_out;
        SOL.n_fevals = nEval;
        SOL.n_Steps  = length(t_out)-1;
        SOL.n_Failed = nFail;
        varargout{1} = SOL;
    else
        plot(t_out,x_out,'-o')
        xlim([tspan(1),tspan(end)])
    end
    
    if(do_waitbar)
        delete(wait_h)
    end
    
end

function [mu,min_dt,max_dt,Newton_TOL,Rel_TOL,Abs_TOL,dfun,refine_n,...
    do_plots,do_waitbar] = ...
        getOptions(options)

    mu = 2-sqrt(2);
    min_dt = 16*eps;
    max_dt = inf;
    Newton_TOL = 1e-10;
    Rel_TOL = 1e-3;
    Abs_TOL = 1e-6;
    dfun = [];
    refine_n = 0;
    do_plots = false;
    do_waitbar = false;
    
    m = 1;
    while m < length(options)
        switch options{m}
            case 'dfun'
                dfun = options{m+1};
            case 'mu'
                mu = options{m+1};
                if mu==0 || mu==1
                    error 'Valid values for mu are 0<mu<1'
                end
            case 'refine'
                refine_n = options{m+1};
                if refine_n<0||(floor(refine_n)-refine_n)~=0
                    error('refine must be a positive integer')
                end
            case 'min_dt'
                min_dt = options{m+1};
            case 'max_dt'
                max_dt = options{m+1};
            case 'Newton_TOL'
                Newton_TOL = options{m+1};
            case 'Rel_TOL' 
                Rel_TOL = options{m+1};
            case 'Abs_TOL'
                Abs_TOL = options{m+1};
            case 'Do_Plots'
                do_plots = options{m+1};
            case 'Do_Waitbar'
                do_waitbar = options{m+1};
        end
        m = m + 2;
    end
end

function [t_new,x_new] = refine_grid(t0,x0,t1,x1,dx0,dx1,n)

    if n > 0
        
        x0 = x0';
        x1 = x1';
        
        dx0 = dx0';
        dx1 = dx1';    
        
        h = t1-t0;

        t_new(:,1) = t0:(h/(n+1)):t1;
        t = (t_new-t0)./h;

        h00 = my_horn([2 -3 0 1],t);
        h10 = h*my_horn([1 -2 1 0],t);
        h01 = my_horn([-2 3 0 0],t);
        h11 = h*my_horn([1 -1 0 0],t);

        x_new = h00*x0 + h10*dx0 + h01*x1 + h11*dx1;
    
        t_new = t_new(2:end);
        x_new = x_new(2:end,:);
        
    else
        
        t_new = t1;
        x_new = x1;
        
    end
    
end

function [jac,evals] = jac(fun,dfun,t,x)

    if isempty(dfun)
        n = length(x);
        fx = fun(t,x);
        xn = x;
        jac = zeros(length(x));
        step = 1e-6;
        for m = 1:n
            xn(m) = xn(m) + step;
            jac(:,m) = (fun(t,xn)-fx)/step;
            xn(m) = x(m);
        end
        evals = n + 1;
    else
        if isa(dfun,'function_handle')
            jac = feval(dfun,t,x);
            evals = 1;
        else
            jac = dfun;
            evals = 0;
        end
    end
    
end

function polval = my_horn(P,t)
    polval = zeros(length(t),1);
    polval(:) = P(1);
    for k = 2:length(P)
        polval = t.*polval + P(k);
    end
end