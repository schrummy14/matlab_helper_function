function varargout = rk43(fun,tspan,x0,varargin)
    
    if nargout > 2
        error('rk43 can only output either the solution structure or [t,x]')
    end
    TOL = 1e-3;
    kmin = eps;
    kmax = inf;
    kStart = 16*sqrt(eps);
    autoRefine = 1;
    minPoints = 35;
    
    if nargin > 3
        for m = 1:2:nargin-3
            switch varargin{m}
                case 'TOL'
                    TOL = varargin{m+1};
                case 'kStart'
                    kStart = varargin{m+1};
                case 'kMin'
                    kmin = varargin{m+1};
                case 'kMax'
                    kmax = varargin{m+1};
                case 'Refine'
                    Refine = varargin{m+1};
                    if strcmp(Refine,'auto')
                        autoRefine = 1;
                    else
                        autoRefine = 0;
                    end
                case 'minPoints'
                    minPoints = varargin{m+1};
                otherwise
                    warning([varargin{m+1}, ' is a bad option... Skipping'])
            end
        end
    end

    [A,b,c] = getCoefs;
    x(1,:) = x0;
    t(1,1) = tspan(1);
    done = false;
    y = zeros(length(c),length(x0));
    y(1,:) = fun(tspan(1),x0);
    funEvals = 1;
    failedStep = 0;
    k = kStart;
    steps = 0;
    n = 1;
    while ~done
        if t(n,1) + k > tspan(2)
            k = tspan(2)-t(n,1);
            done = 1;
        end
        
        for m = 2:length(c)
            y(m,:) = fun(t(n,1)+c(m)*k, x(n,:)+k*A(m,:)*y);
        end
        funEvals = funEvals + length(c) - 1;
        Err = norm(k*(b(1,:) - b(2,:))*y,inf);
        
        if Err<TOL
            if autoRefine
                Refine = ceil(k*minPoints/diff(tspan));
            end
            [tnew,xnew] = ...
            denseOut([t(n,1);t(n,1)+k],x(n,:),k*y,Refine);
            t(n+1:n+Refine,1) = tnew(2:end,1);
            x(n+1:n+Refine,:) = xnew(2:end,:);
            y(1,:) = y(end,:);
            n = n + Refine;
            steps = steps + 1;
        else
            done = false;
            failedStep = failedStep + 1;
        end
        k =  min(kmax,max(kmin,0.9*k*(TOL/(Err + eps))^(1/4)));
            
    end
    
    if nargout > 1
        varargout{1} = t;
        varargout{2} = x;
    elseif nargout > 0
        Sol.x = x;
        Sol.t = t;
        Sol.funEvals = funEvals;
        Sol.steps = steps;
        Sol.failedStep = failedStep;
        varargout{1} = Sol;
    else
        figure
        plot(t,x,'-o')
        ylabel('x(t)')
        xlabel('t')
    end
end

function [A,b,c] = getCoefs()

    A = [
        0       0       0       0       0
        1/3     0       0       0       0
       -1/3     1       0       0       0
        1      -1       1       0       0
        1/8     3/8     3/8     1/8     0
        ];
    
    c = sum(A,2);
    b = [
        1/8     3/8     3/8     1/8     0
        1/24    5/8     1/8    -1/8     1/3
        ];
end

function [tnew,xnew] = denseOut(t,x,f,Refine)

    tnew = linspace(t(1),t(2),Refine+1).';
    dift = diff(t);
    tt = (tnew-t(1))/dift;

    tt2 = tt.^2;
    tt3 = tt.^3;
    
    A = 24-45*tt+22*tt2;
    B = 3-2*tt;   

    xnew = ones(size(tt))*x(1,:) + ...
        tt/24.*A*f(1,:) + 5*tt2/8.*B*f(2,:) + tt2/8.*B*(f(3,:)-f(4,:)) + tt3/3*f(5,:);
    
end




















