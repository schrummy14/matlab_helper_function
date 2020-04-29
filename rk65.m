function varargout = rk65(fun,tspan,x0,varargin)
    
    SolverName = 'rk65';
    if nargout > 2
        error('rk65 can only output either the solution structure or [t,x]')
    end
    TOL = 1e-4;
    kmin = eps;
    kmax = inf;
    kStart = diff(abs(tspan))*sqrt(eps);
    autoRefine = 1;
%     minPoints = 25;
    userSetKmin = false;
    
    if nargin > 3
        for m = 1:2:nargin-3
            switch varargin{m}
                case 'TOL'
                    TOL = varargin{m+1};
                case 'kStart'
                    kStart = varargin{m+1};
                case 'kMin'
                    kmin = varargin{m+1};
                    userSetKmin = true;
                case 'kMax'
                    kmax = varargin{m+1};
                case 'Refine'
                    Refine = varargin{m+1};
                    if strcmp(Refine,'auto')
                        autoRefine = 1;
                    else
                        autoRefine = 0;
                    end
%                 case 'minPoints'
%                     minPoints = varargin{m+1};
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
    nfevals = 1;
    nfailed = 0;
    k = kStart;
    nsteps = 0;
    n = 1;
    failedOneStep = 0;
    while ~done
        if t(n,1) + k > tspan(2)
            k = tspan(2)-t(n,1);
            done = 1;
        end
        
        for m = 2:length(c)
            y(m,:) = fun(t(n,1)+c(m)*k, x(n,:)+k*A(m,:)*y);
        end
        nfevals = nfevals + length(c) - 1;
        Err = norm(k*(b(1,:) - b(2,:))*y,inf);
        
        if Err<TOL
            if autoRefine
%                 Refine = min(max(1,ceil(k*minPoints/diff(tspan))),6);
                c1 = ceil(10*k/diff(tspan));
                c2 = ceil(10/max(abs(x0))*(max(abs(k*b(2,:)*y))));
                Refine = min(50,max(c1,c2));
            end
            [tnew,xnew] = ...
            denseOut([t(n,1);t(n,1)+k],x(n,:),k*y,Refine);
            t(n+1:n+Refine,1) = tnew(2:end,1);
            x(n+1:n+Refine,:) = xnew(2:end,:);
            y(1,:) = y(end,:);
            n = n + Refine;
            nsteps = nsteps + 1;
            failedOneStep = 0;
        else
            done = false;
            nfailed = nfailed + 1;
            failedOneStep = failedOneStep + 1;
        end
        if ~userSetKmin
            kmin = eps(t(n,1));
        end
        k = min(kmax,max(kmin,0.9*k*(TOL/(Err + eps))^(1/6)));
        
        if failedOneStep > 10
            disp('Failed to integrate the solution, please try a "Stiff" solver...')
            figure
            plot(t,x);
            return
        end
        
    end
    
    if nargout > 1
        varargout{1} = t;
        varargout{2} = x;
    elseif nargout > 0
        Sol.solver = SolverName;
        Sol.extdata.odefun = fun;
        Sol.extdata.options = [];
        Sol.extdata.varargin = varargin;
        Sol.t = t;
        Sol.x = x;
        Sol.stats.nsteps = nsteps;
        Sol.stats.nfailed = nfailed;
        Sol.stats.nfevals = nfevals;
        Sol.idata = [];
        varargout{1} = Sol;
    else
        figure
        if length(t)<100
            plot(t,x,'-o')
        else
            plot(t,x)
        end
        ylabel('x(t)')
        xlabel('t')
    end
end

function [A,b,c] = getCoefs()

    A = zeros(10);
    A(2,1) = 1/32;
    A(3,1:2) = 1/72*[1 2];
    A(4,1:3) = 1/64 * [1 0 3];
    A(5,1:4) = 1/125 * [53 0 -204 176];
    A(6,1:5) = [1/96 0 0 4/33 125/1056];
    A(7,1:6) = [-19/24 0 0 64/33 -875/264 8/3];
    A(8,1:7) = [987/320 0 0 -36/5 975/64 -459/40 177/160];
    A(9,1:8) = [-1238/105 0 0 33536/1155 -14900/231 1796/35 -21/5 8/7];
    A(10,1:9)= 1/90*[7 0 0 0 0 32 12 32 7];
    
    c = sum(A,2);
    b = 1/450*[
        35 0 0 0 0 160 60 160 28  7
        35 0 0 0 0 160 60 160 42 -7
        ];
end

function [tnew,xnew] = denseOut(t,x,f,Refine)

    tnew = linspace(t(1),t(2),Refine+1).';
    dift = diff(t);
    tt = (tnew-t(1))/dift;

    C = [
          1 0 0 0 0    0   0    0  0  0  
        -25 0 0 0 0   48 -36   16 -1 -2
         70 0 0 0 0 -208 228 -112  1 21
        -10 0 0 0 0   36 -48   28  1 -7
          4 0 0 0 0  -16  24  -16 -1  5
          ];
      Cf = C*f;

      tt2 = tt.*tt;
      tt3 = tt.*tt2;
      tt4 = tt.*tt3;
      tt5 = tt.*tt4;

      xnew = ones(size(tt))*x(1,:) + ...
          (tt*Cf(1,:) + 1/6*tt2*Cf(2,:) + 1/9*tt3*Cf(3,:) + 2/3*tt4*Cf(4,:) + 8/15*tt5*Cf(5,:));
end






