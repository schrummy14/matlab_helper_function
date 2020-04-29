function [t,x] = rk78(fun,tspan,x0,varargin)

%   [t,x] = rk78(fun,tspan,x0,varargin)
%   rk78 is an embeded Runge Kutta method of orders 7 and 8 as described by
%   Dormand and Prince outlined in their paper High order embedded Runge-
%   Kutta formulae in Journal of Computational and Applied Mathematics.
%   
%   rk78 is a ode initial value systems of ODEs solver which uses the 8th
%   order method to bound the error of the 7th order. This allows the use
%   of an adaptive step.
%
%   rk78's inputs are as follows:
%       fun is the system of odes you would like to solve. They must be
%           function handels of the form sys(t,x). fun must return a column
%           vector.
%       tspan is the time of which you would like to solve the system.
%           tspan takes the form tspan(StartTime EndTime). If tspan is
%           given a sigle value, it will be treated as the endtime with the
%           start time being 0.
%       x0 is the row vector containing the initial conditions.
%       varargin is a list of options that take the form ('name',value).
%           The following options are currently available;
%             TOL       :The tolerence used in estimating the bounded error 
%                           of the 7th order method.
%             kMin      :The minimum step size that the solver will use.
%                           Note: This can cause instabilities and it is
%                           suggested to first change TOL to a lower 
%                           tolerence.
%             kmax      :The maximum step size the solver will use. 
%             fixStep   :The solver will use a fixed step. This fixed step
%                           is determined by the difference of tspan 
%                           divieded by the number of steps given. 
%             exitCon   :This is a value that the solver will be looking 
%                           for to satisfy. Once this condition is found to
%                           be true, the solver will stop. This does not 
%                           increase the tspan value to allow longer 
%                           interation time. 

    if size(x0,1)~=1 % If x0 is given as a column vector, flip it.
        x0 = x0.';
    end
    if numel(tspan) == 1    % If Tspan only holds one value, assume it's
        tspan(2) = tspan;   % the end value and start at zero.
        tspan(1) = 0;
    end
    
    % Stock values
    kmin = eps(1+tspan(1));
    kmax = inf;
    TOL = 1e-10;
    useFixStep = 0;
    exitCon = 0;
    exCond = @(x)0;
    
    if nargout == 0 
        do_plot = true;
    else
        do_plot = false;
    end
    
    % If options are passed, check them
    if nargin > 3
        for n = 1:2:numel(varargin)
            switch varargin{n}
                case 'TOL'
                    TOL = varargin{n+1};
                case 'kMin'
                    kmin = varargin{n+1};
                case 'kMax'
                    kmax = varargin{n+1};
                case 'fixStep'
                    useFixStep = 1;
                    N = varargin{n+1};
                    kmin = diff(tspan)/N;
                case 'exitCon'
                    exitCon = 1;
                    exCond = varargin{n+1};
                    conTry = 0;
                otherwise
                    error([varargin{n},' was not found. Current options are: ',...
                        'TOL, kMin, kMax, fixStep, exitCon.'])
            end
        end
    end
    
    % Initialize the method
    k = 16*kmin;   % Start with a small step.
    x(1,:) = x0;
    t(1,1) = tspan(1);
    n = 1;    
    done = 0;
    kErrorLoc = inf;
    
    if kmax < kmin  % Make sure kmax is not smaller than kmin
        error('The min step size is larger than the max step size...')
    end
    
    
    % Set up the runge kutta 8(7) method
    [A,b,c] = getCoefs;
    
    y = zeros(13,length(x0));
    kMinError = 0;
    if useFixStep
        for n = 1:N
            y = zeros(13,length(x0));
            for m = 1:size(y,1)
                y(m,:) = fun(t(n,1)+c(m)*k, x(n,:)+k*A(m,:)*y);
            end
            x(n+1,:) = x(n,:) + k*b(1,:)*y;
            t(n+1,:) = tspan(1) + n*k;
        end
    else
        while ~done

            if t(n,1) + k > tspan(2)
                k = tspan(2) - t(n,1);
                done = 1;
            end

            for m = 1:size(y,1)
                y(m,:) = fun(t(n,1)+c(m)*k, x(n,:)+k*A(m,:)*y);
            end
            if sum(sum(isnan(y)))>0 || sum(sum(isinf(y)))>0
                warning('Answer contains NAN or Inf Values. Stopping Solver.')
                done = 1;
                tau = 0;
            else
                x8 = x(n,:) + k*b(1,:)*y;
                tau = norm(k*(b(1,:)-b(2,:))*y,2);
            end
            
            if tau < TOL || k <= kmin
                kmin = eps(max(x8(:)));
                x(n+1,:) = x8;
                t(n+1,1) = t(n,1) + k;
                
                if do_plot
                    plot(t,x,'-o')
                    xlim([tspan(1),tspan(end)])                    
                    drawnow
                end
                
                if k<=kmin && ~done && n>1
                    kMinError = 1;
                    kErrorLoc = n;
                end
                if exitCon && exCond(x8)
                    if kErrorLoc == n,kMinError = 0;end
                    if conTry > 15
                        done = 1;
                        exitCon = 0;
                    else
                        conTry = conTry + 1;
                        n = n - 1;
                        tau = 10^6*TOL;
                    end
                else
                    conTry = 0;
                end
                n = n +1;
            elseif done
                done = 0;            
            end
            k = min(kmax,max(kmin,0.95*k*(TOL/(tau+eps))^(1/8))); 
        end
    end
    if kMinError
        warning('Minimum k value has been reached, solution may not be correct...')
    end
    if exitCon
        warning('Solver did not reach the exit condition...')
    end
    switch nargout
        case 0
            figure
            plot(t,x,'-o')
            t = [t,x];
        case 1
            t = [t,x];
    end
end

function [A,b,c] = getCoefs

    Atop = [ % The top half of the A matrix from the butcher table
        0           0       0       0               0           0               0               0               0               0               0           0               0
        1           0       0       0               0           0               0               0               0               0               0           0               0
        1           1       0       0               0           0               0               0               0               0               0           0               0
        1           0       3       0               0           0               0               0               0               0               0           0               0
        5           0      -75      75              0           0               0               0               0               0               0           0               0
        3           0       0       3               3           0               0               0               0               0               0           0               0
        29443841    0       0       77736538       -28693883    23124283        0               0               0               0               0           0               0
        16016141    0       0       61564180        22789713    545815736      -180193667       0               0               0               0           0               0
        39632708    0       0      -433636366      -421739975   100302831       790204164       800635310       0               0               0           0               0
        246121993   0       0      -37695042795    -309121744  -12992083        6005943493      393006217       123872331       0               0           0               0
       -1028468189  0       0       8478235783      1311729495 -10304129995    -48777925059     15336726248    -45442868181     3065993473      0           0               0
        185892177   0       0      -3185094517     -477755414  -703635378       5731566787      5232866602     -4093664535      3962137247      65686358    0               0
        403863854   0       0      -5068492393     -411421997   652783627       11173962825    -13158990841     3936647629     -160528059       248638103   0               0
        ];
    
    Abot = [ % Bottom of the A matrix
        1           1       1       1               1           1               1               1               1               1               1           1               1
        18          1       1       1               1           1               1               1               1               1               1           1               1
        48          16      1       1               1           1               1               1               1               1               1           1               1
        32          1       32      1               1           1               1               1               1               1               1           1               1
        16          1       64      64              1           1               1               1               1               1               1           1               1
        80          1       1       16              20          1               1               1               1               1               1           1               1
        614563906   1       1       692538347       1125000000  1800000000      1               1               1               1               1           1               1
        946692911   1       1       158732637       633445777   2771057229      1043307555      1               1               1               1           1               1
        573591083   1       1       683701615       2616292301  723423059       839813087       3783071287      1               1               1           1               1
        1340847787  1       1       15268766246     1061227803  490766935       2108947869      1396673457      1001029789      1               1           1               1
        846180014   1       1       508512852       1432422823  1701304382      3047939560      1032824649      3398467696      597172653       1           1               1
        718116043   1       1       667107341       1098053517  230739211       1027545527      850066563       808688257       1805957418      487910083   1               1
        491063109   1       1       434740067       543043805   914296604       925320556       6184727034      1978049680      685178525       1413531060  1               1
        ];
    
    A = Atop./Abot; % Combine
        
    
    btop = [ % Start of the b matrix. The 8th order method is the 1st row.
        13451932    0       0       0               0           -808719846      1757004468      656045339       -3867574721     465885868       53011238     2              0
        14005451    0       0       0               0           -59238493       181606767       561292985       -1041891430     760417239       118820643    -528747749     1
        ];
    
    bbot = [
        455176623   1       1       1               1           976000145       5645159321      265891186       1518517206      322736535       667516719   45              1
        335480064   1       1       1               1           1068277825      758867731       797845732       1371343529      1151165299      751138087   2220607170      4  
        ];
    
    b = btop./bbot;
    
    ctop = [ % Start of the c vector
        0           1       1       1               5           3               59              93              5490023248      13              1201146811  1               1
        ];
    cbot = [
        1           18      12      8               16          8               400             200             9719169821      20              1299019798  1               1
        ];
    
    c = ctop./cbot;
    
end

