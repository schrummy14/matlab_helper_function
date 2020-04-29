function varargout = velocity_verlet(fun,tspan,x0,v0,opts)

    if nargin < 4
        opts = [];
    end
    
    opts = set_options(opts);
    
    if opts.FixedStep
        NumSteps = ceil(diff(tspan)/opts.dt);
        t = nan(NumSteps+1,1);
        x = nan(NumSteps+1,numel(x0));
        v = x;
        t(1) = tspan(1);
        x(1,:) = x0;
        v(1,:) = v0;
        
        fun_val = fun(t(1),x(1,:),v(1,:));
        for k = 1:NumSteps
            v(k+1,:) = v(k,:) + 0.5*opts.dt*fun_val;
            x(k+1,:) = x(k,:) + opts.dt*v(k+1,:);
            fun_val = fun(t(k),x(k,:),v(k,:));
            v(k+1,:) = v(k,:) + 0.5*opts.dt*fun_val;
            t(k+1) = t(k) + opts.dt;
        end
        
    else
        error 'Not yet implemented'
        % Assume spring dampener as a starting point
        w = sqrt(fun(0,1,0));
        dt_guess = 1/w;
        CHUNK = 200;
        t = nan(CHUNK,1);
        x = nan(CHUNK,numel(x0));
        v = nan(CHUNK,numel(v0));e
    end
    
    switch nargout
        case 0
            plot(t,x,t,v)
        case 1
            Sol.t = t;
            Sol.x = x;
            Sol.v = v;
            varargout{1} = Sol;
        case 2
            varargout{1} = t;
            varargout{2} = [x,v];
        case 3
            varargout{1} = t;
            varargout{2} = x;
            varargout{3} = v;
    end
    
end

function opts = set_options(opts)
    
    if isempty(opts)
        opts.FixedStep = false;
    end
end