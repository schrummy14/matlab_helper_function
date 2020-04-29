function [t,x] = rk4(fun,tspan,x0,varargin)
    
    method = 'default';
    TOL = 1e-3;
    autoRefine = true;
    if nargin>3
        gg = 1;
        if isnumeric(varargin{1})
            N = varargin{1};
            fixStep = 1;
            gg = 2;
        else
            fixStep = 0;
        end
        for m=gg:2:length(varargin)
            switch varargin{m}
                case 'method'
                    method = varargin{m+1};
                case 'fixStep'
                    fixStep = varargin{m+1};
                case 'TOL'
                    TOL = varargin{m+1};
                case 'Refine'
                    Refine = varargin{m+1};
                    autoRefine = false;
                    if strcmp('auto',Refine)
                        autoRefine = true;
                    end
                otherwise
                    warning([varargin{m},' is not an option...skipping'])
            end
        end
    elseif nargin>2
        fixStep = false;
    end
    if numel(tspan)<2
        tspan(2) = tspan;
        tspan(1) = 0;
    end
    if size(x0,1)>1,x0=x0.';end
    
    
    [A,b,c] = getCoefs(method);
    y = zeros(length(c),length(x0));
    
    if fixStep
        k = diff(tspan)/N;
        x = zeros(N,length(x0));
        t = zeros(N,1);
        x(1,:) = x0;
        t(1,1) = tspan(1);
        for n = 1:N
            for m = 1:length(c)
                y(m,:) = fun(t(n,1)+k*c(m),x(n,:)+k*A(m,:)*y);
            end
            x(n+1,:) = x(n,:) + k*b*y;
            t(n+1,1) = tspan(1) + n*k;
        end
    else
        x(1,:) = x0;
        t(1,1) = tspan(1);
        n = 1;
        k = diff(tspan)/100;
        kmax = diff(tspan)/10;
        done = false;
        while ~done
            if t(n,1)+k>tspan(2)
                k = tspan(2)-t(n,1);
                done = true;
            end
            kh = k/2;
            for m = 1:length(c)
                y(m,:) = fun(t(n,1)+kh*c(m),x(n,:)+kh*A(m,:)*y);
            end
            x1(1,:) = x(n,:) + kh*b*y;
            for m = 1:length(c)
                y(m,:) = fun(t(n,1)+kh+kh*c(m),x1+kh*A(m,:)*y);
            end
            x2 = x1 + kh*b*y;
            for m = 1:length(c)
                y(m,:) = fun(t(n,1)+k*c(m),x(n,:)+k*A(m,:)*y);
            end
            xnew = x(n,:) + k*b*y;
            
            err = norm(xnew-x2,inf);
            if err<TOL
                tnew = t(n,1)+k;
                if autoRefine
                    Refine = max(1,ceil(k*35/diff(tspan)));
                end
                if Refine > 1
                    [tnew,xnew] = denseOutput([t(n,1);tnew],...
                        [x(n,:);xnew],[y(1,:);y(end,:)],Refine);
                end
                x(n+1:n+Refine,:) = xnew;
                t(n+1:n+Refine,1) = tnew;
                n = n + Refine;
            else
                done = false;
            end
            k = min(kmax,0.95*k*(TOL/(err+eps))^(1/4));
        end
    end
    
    if nargout == 0
        plot(t,x,'o-')
    end
    
end

function [A,b,c] = getCoefs(method)

    switch method
        case 'default'
            A = [
                0    0    0    0
                1/2  0    0    0
                3/10 1/5  0    0
                0   -3/2  5/2  0
                ];
            b = [1/6 -1/6 5/6 1/6];
            c = [0   1/2  1/2  1];
        case 'classic'
            A = [
                0   0   0   0
                1/2 0   0   0
                0   1/2 0   0
                0   0   1   0
                ];
            b = [1/6 1/3 1/3 1/6];
            c = [0  1/2 1/2 1];
        otherwise
            warning('Bad method passed.... using default...')
            [A,b,c] = getCoefs('default');
    end
end

function [tnew,xnew] = denseOutput(t,x,y,Refine)

    tnew1 = linspace(t(1),t(end),Refine+1).';
    
    tdiff = t(end)-t(1);
    tt(:,1) = (tnew1-t(1))./tdiff;
    
    tt3 = tt.^3;
    tt2 = tt.^2;
    h00 = 2*tt3 - 3*tt2 + 1;
    h10 = tt3 - 2*tt2 + tt;
    h01 = -2*tt3 + 3*tt2;
    h11 = tt3 - tt2;
    
    p0 = x(1,:);
    p1 = x(end,:);
    m0 = y(1,:);
    m1 = y(end,:);
    
    xnew1 = h00*p0 + h10*tdiff*m0 + h01*p1 + h11*tdiff*m1;
    
    xnew = xnew1(2:end,:);
    tnew = tnew1(2:end,1);
    
end
                