function varargout = rk56(fun,tspan,x0,varargin)

    if size(x0,1)~=1
        x0 = x0';
    end
    if numel(tspan) == 1
        tspan(2) = tspan;
        tspan(1) = 0;
    end
    
    refined = 3;
    kmin = sqrt(eps);
    kmax = inf;
    TOL = 1e-3;
    nStep = 0;
    
    if nargin > 3
        for n = 1:2:numel(varargin)
            switch varargin{n}
                case 'TOL'
                    TOL = varargin{n+1};
                case 'kMin'
                    kmin = varargin{n+1};
                case 'kMax'
                    kmax = varargin{n+1};
                case 'refined'
                    refined = varargin{n+1};
                otherwise
                    error([varargin{n},' was not found. Current options are: ',...
                        'TOL, kMin, kMax.'])
            end
        end
    end
    
    k = kmin;
    x(1,:) = x0;
    t(1,1) = tspan(1);
    n = 1;
    done = false;
    
    if kmax < kmin
        error('The min step size is larger than the max step size...')
    end
    
    [A,b,c] = getCoefs();
    stages = length(c);
    
    y = zeros(stages,length(x0));
    y(1,:) = fun(tspan(1),x0);
    funEvals = 1;
    stepFail = 0;
    while ~done
        if t(n,1) + k > tspan(2)
            k = tspan(2)-t(n,1);
            done = true;
        end
        
        for m = 2:stages
            y(m,:) = fun(t(n,1)+c(m)*k, x(n,:)+k*A(m,:)*y);
        end
        funEvals = funEvals + 6;
        x6 = x(n,:) + k*b(1,:)*y;
        x5 = x(n,:) + k*b(2,:)*y;
        err = norm(x6-x5,inf);
        
        if err<TOL || k<=kmin
            x(n+1,:) = x6;
            t(n+1,1) = t(n,1) + k;
            if refined>1
                [tnew,xnew] = denseOutput(t(n:n+1,1),x(n:n+1,:),[y(1,:);y(end,:)],refined);
                t(n+1:n+length(tnew)-1,1) = tnew(2:end,1);
                x(n+1:n+length(xnew)-1,:) = xnew(2:end,:);
                n = n + refined;
            else
                n = n + 1;
            end
            y(1,:) = y(end,:);
            nStep = nStep + 1;
        else
            done = false;
            stepFail = stepFail + 1;
        end
        
        k = min(kmax,max(kmin,0.9*k*(TOL/(err + eps))^(1/6)));
    end
    if nargout>1
        if nargout>4,varargout{5} = nStep;end
        if nargout>3,varargout{4} = stepFail;end
        if nargout>2,varargout{3} = funEvals;end
        if nargout>1,varargout{1} = t;varargout{2} = x;end
    elseif nargout > 0
        Sol.t=t;
        Sol.x=x;
        Sol.funEvals=funEvals;
        Sol.stepFail=stepFail;
        Sol.nStep = nStep;
        varargout{1} = Sol;
    else
        plot(t,x,'-o')
    end
end

function [A,b,c] = getCoefs()
    
    An = [
           0           0           0           0           0           0           0
     2507295           0           0           0           0           0           0
   -63776695    28685351           0           0           0           0           0
     3533952   -15757219    25824805           0           0           0           0
   142856474   -27679525    -8953421    39475178           0           0           0
  -102080007   263098293    42259935   -61207707    21727210           0           0
     2283975    12787441     5280135    35387115    17062480     8523395           0
     ];

    Ad = [
           1           1           1           1           1           1           1
    16027766           1           1           1           1           1           1
   249919681    40448660           1           1           1           1           1
    17171693   118540766    40718209           1           1           1           1
   147722655    38795022    76268989    52297341           1           1           1
    51834440   118101608    23319751    37058539    38200917           1           1
    50504276    50938375    17610394   164314183   118677213   190365469           1
    ];

    bn = [
     2283975    12787441     5280135    35387115    17062480     8523395           0
     1519879    22723713     4341101     8524261     4702885     2590840        1137
     ];
    bd = [
    50504276    50938375    17610394   164314183   118677213   190365469           1
    29412040    96181777    13702901    42022552    31869831    57799907    44379772
    ];

    cn = [
           0     2507295    38810536    38613965    65532703   156327905           1
           ];
    cd = [
           1    16027766    85487551    54608393    73549072   158276552           1
           ];
       
    A = An./Ad;
    b = bn./bd;
    c = cn./cd;
end

function [tnew,xnew] = denseOutput(t,x,f,refined)

    xx = linspace(t(1),t(end),(2+(refined-1))).';
    xV = xx(2:end-1);
    tdiff = t(end)-t(1);
    tt(:,1) = (xV-t(1))./tdiff;
    
    tt3 = tt.^3;
    tt2 = tt.^2;
    h00 = 2*tt3 - 3*tt2 + 1;
    h10 = tt3 - 2*tt2 + tt;
    h01 = -2*tt3 + 3*tt2;
    h11 = tt3 - tt2;
    
    p0 = x(1,:);
    p1 = x(end,:);
    m0 = f(1,:);
    m1 = f(end,:);
    
    xnew1 = h00*p0 + h10*tdiff*m0 + h01*p1 + h11*tdiff*m1;
    tnew1 = xV;
    
    tnew = [t(1);tnew1;t(end)];
    xnew = [x(1,:);xnew1;x(end,:)];  
    
end












