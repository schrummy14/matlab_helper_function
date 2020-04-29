function [x, fx, stats] = Mylsqnonlin(Fun,x0,opts)

    if nargin < 3
        opts = [];
        if nargin == 0
%             Fun = @(x)x(1)*exp(-(x(1)^2 + x(2)^2)) + (x(1)^2 + x(2)^2)/20;
            Fun = @(x)((x(1)-sqrt(2)).^2 + (x(2)+sqrt(5)).^2);
            x0 = [-1,3];
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
    
    
end

function opts = get_defaults(x0,opts_in)

    opts.do_plot = false;
    opts.P = 0.5*ones(size(x0));
    opts.lb = -inf(size(x0));
    opts.ub =  inf(size(x0));
    opts.T = 1e-6*x0;
    opts.display = 'final';
    
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