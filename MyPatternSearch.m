function [x, fx, stats] = MyPatternSearch(Fun,x0,opts)

    if nargin < 3
        opts = [];
        if nargin == 0
            Fun = @(x)x(1)*exp(-(x(1)^2 + x(2)^2)) + (x(1)^2 + x(2)^2)/20;
            x0 = [1 1];
            opts.Tol = 1e-10;
            opts.lb = -inf(numel(x0),1);
            opts.ub = inf(numel(x0),1);
            opts.display = 'debug';
        end
    end
    
    % Initialization
    [M,N] = size(x0);
    x0 = x0(:);
    opts = get_defaults(x0,opts);
    N_vars = length(x0);
    fbest = Fun(x0);
    opts.fevals = 1;
    x1 = x0;
    
    if strcmpi(opts.display,'debug')
        fprintf("   Iteration\tFunc-count\t\tf(x)\t\tMeshSize\t\tMethod\n");
        opts.its = 0;
    end
    if opts.do_plot
        opts.x_path(1,:) = x0;
        opts.f_path(1,1) = fbest;
    end
    while true
        % Start Restart
        [x1,fbest,opts] = start_restart(Fun,x1,fbest,N_vars,opts);
        if any(opts.P < opts.T)
            break
        end

        % Pattern Move
        [x1,x0,fbest,opts] = pattern_move(Fun,x1,x0,fbest,opts);
    end
    if opts.do_plot
        for k = 1:N_vars
            opts.x_path(:,k) = opts.x_path(:,k)/(abs(opts.x_path(end,k)));
        end
        plot(opts.f_path,opts.x_path,'x-')
    end
    x = reshape(x1,M,N);
    fx = fbest;
    stats = opts;
    
end

function [x1,x0,fbest,opts] = pattern_move(Fun,x1,x0,fbest,opts)
    x2 = x0 + opts.a*(x1 - x0);
    meets_conds = true;
    for k = 1:length(x0)
        if x2(k) < opts.ub(k)
            if x2(k) > opts.lb(k)
                continue;
            else
                meets_conds = false;
            end
        else
            meets_conds = false;
        end
    end
    CurF = Fun(x2);
    opts.fevals = opts.fevals + 1;
    if strcmpi(opts.display,'debug')
        opts.its = opts.its + 1;
        if mod(opts.its,opts.display_k) == 0
            print_debug(opts,fbest,'Moving');
        end
    end 
    if CurF <= fbest && meets_conds
        fbest = CurF;
        x0 = x1;
        x1 = x2;
        if opts.do_plot
            opts.x_path(end+1,:) = x0;
            opts.f_path(end+1,1) = fbest;
        end
        [x1,x0,fbest,opts] = pattern_move(Fun,x1,x0,fbest,opts);
    else
        x0 = x1;
    end
end
    

function [x1,fbest,opts] = start_restart(Fun,x1,fbest,N_vars,opts)
    has_improved = false;
    AtBound = false;
    x1_old = x1;
    for k = 1:N_vars
        x1(k) = x1(k) + opts.P(k);
        if x1(k) > opts.ub(k)
            x1(k) = opts.ub(k);
            AtBound = true;
        end
        CurF = Fun(x1);
        opts.fevals = opts.fevals + 1;
        if CurF < fbest
            fbest = CurF;
            has_improved = true;
        else
            x1(k) = x1(k) - 2*opts.P(k);
            if x1(k) < opts.lb(k)
                x1(k) = opts.lb(k);
                AtBound = true;
            end
            CurF = Fun(x1);
            opts.fevals = opts.fevals + 1;
            if CurF < fbest
                fbest = CurF;
                has_improved = true;
            elseif ~AtBound
                x1(k) = x1(k) + opts.P(k);
            end
        end
    end
    if strcmpi(opts.display,'debug')
        opts.its = opts.its + 1;
        if mod(opts.its,opts.display_k) == 0
            print_debug(opts,fbest,'Search');
        end
    end
    if ~has_improved
        x1 = x1_old;
        opts.P = 0.5 * opts.P;
        if ~any(opts.P < opts.T)
            [x1,fbest,opts] = start_restart(Fun,x1,fbest,N_vars,opts);
        end
    end
    if opts.do_plot
        opts.x_path(end+1,:) = x1;
        opts.f_path(end+1,1) = fbest;
    end
end

function print_debug(opts,fbest,method)
    fprintf("%12g%14g%18g%20g\t\t%s\n",...
        opts.its,opts.fevals,fbest,min(opts.P),method);
end

function opts = get_defaults(x0,opts_in)

    opts.do_plot = false;
    opts.P = 0.5*ones(size(x0));
    opts.lb = -inf(size(x0));
    opts.ub =  inf(size(x0));
    opts.T = 1e-6*x0;
    opts.a = 2.0;
    opts.display = 'final';
    opts.display_k = 10;
    
    if isempty(opts_in)
        return
    end
    
    Passed_Names = fieldnames(opts_in);
    
    for k = 1:length(Passed_Names)
        switch Passed_Names{k}
            case 'Tol'
                opts.T = opts_in.Tol;
            otherwise
                opts.(Passed_Names{k}) = opts_in.(Passed_Names{k});
        end
    end

end