function [t,x] = TRBDF2s(fun,tspan,x0,varargin)

% [t,x] = TRBDF2s(fun,tspan,x0,varargin)
% TRBDF2s is an implicit embedded ode solver. The solver takes a half step
% using the trapezoidal rule and finishes the step using a 2nd order
% backward difference method. Durring this, the midpoint method is also
% calculated. This allows the error of the midpoint methoded to be bounded
% by the trapezoidal and backward difference methods. By default, the
% solver computes the jacobian numericaly using forward difference.
%
% fun is the system of odes you would like to solve. fun must return a
%   column vector
%
% tspan is the time interval that the simulation will cover. tspan must
%   either be a vector containing the start and stop times or must just
%   contain the stop time. The start time is assumed to be 0.
%
% x0 is the initial starting points for the solver. x0 must be a row
%   vector.
%
% varargin is used to pass additional options in the form 'name', value.
%   Current options include:
%
%       errorTol:   The upper bound on error that must be met for the 
%                   solver to take another step.
%
%       solTol:     The bound on error when using the newton solver.
%
%       maxIts:     The maximum number of tries to solve the newton system.
%
%       kMin:       The minimum step size the solver will use.
%
%       kMax:       The maximum step size the solver will use.
%
%       fixStep:    Will make the solver use N uniformally spaced points
%                   to solve the system.
%
%       jac:        Is the function handld containing the jaccobian matrix
%                   of the system.
%       
%       exitCon:    This is a value that the solver will be looking for to 
%                   satisfy. Once the condition has been met, the solver
%                   will terminate. This will not increase the tspan value
%                   to allow longer interation time.
%
%       Note: Using kMin and kMax is discouraged. The user should adjust 
%       the tolerances to help with speed.

    TOL1 = 1e-3;
    TOL2 = 1e-10;
    MaxIt = 500;
    kmax = inf;
    kmin = 1e-12;
    useFixStep = 0;
    exitCon = 0;
    exCond = @(x)0;
    
    if nargin > 3
        for n = 1:2:numel(varargin)
            switch varargin{n}
                case 'errorTol'
                    TOL1 = varargin{n+1};
                case 'solTol'
                    TOL2 = varargin{n+1};
                case 'maxIts'
                    MaxIt = varargin{n+1};
                case 'kMin'
                    kmin = varargin{n+1};
                case 'kMax'
                    kmax = varargin{n+1};
                case 'fixStep'
                    N = varargin{n+1};
                    useFixStep = 1;
                case 'Jac'
                    if isa(varargin{n+1},'function_handle')
                        Jac = varargin{n+1};
                    else
                        warning(['Jac must be a function handle,',...
                            'continuing with a numerical jacobian'])
                    end
                case 'exitCon'
                    exitCon = 1;
                    exCond = varargin{n+1};
                    conTry=0;
                otherwise
                    error([varargin{n},' was not found'])
            end
        end
    end
    
    x(1,:) = x0;
    t(1,1) = tspan(1);
    k = kmin;
    done = 0;
    n = 1;
    if useFixStep
        k = diff(tspan)/N;
        for n = 1:N;
            
            fFound = 0;
            j = 0;
            xguess = x(n,:);
            
            while ~fFound && j <MaxIt
                j = j + 1;
                f = x(n,:)' + k/4 * (fun(t(n,1),x(n,:)) + fun(t(n,1)+k/2,xguess)) - xguess';

                try
                    g = k/4 * Jac(t(n,1)+k/2,xguess) - eye(length(xguess));
                catch
                    g = k/4 * jac(fun,t(n,1) + k/2,xguess) - eye(length(xguess));
                end

                delta = g\f;

                if max(abs(delta)) < TOL2
                    fFound = 1;
                    xmid = xguess - delta';
                else
                    xguess = xguess - delta';
                end
            end

            if j >= MaxIt
                error('Max Iterations on Solver has been Exceeded')
            end
            
            fFound = 0;
            j = 0;
            xguess = x(n,:);

            while ~fFound && j < MaxIt

                j = j + 1;
                f = k*fun(t(n,1)+k,xguess)-x(n,:)'+4*xmid'-3*xguess';

                try
                    g = k*Jac(t(n,1)+k,xguess) - 3*eye(length(xguess));
                catch
                    g = k*jac(fun,t(n,1)+k,xguess)-3*eye(length(xguess));
                end

                delta = g\f;

                if max(abs(delta)) < TOL2
                    fFound = 1;
                    x(n+1,:) = xguess - delta';
                    t(n+1,1) = tspan(1)+n*k;
                else
                    xguess = xguess - delta';
                end
            end

            if j >= MaxIt
                error('Max Iterations on Solver has been Exceeded')
            end
            
        end              
    else
        while ~done

            if t(n,1) + k > tspan(2)
                k = tspan(2) - t(n,1);
                done = 1;
            end

            fFound = 0;
            j = 0;
            xguess = x(n,:);

            while ~fFound && j < MaxIt

                j = j + 1;
                f = x(n,:)'+k/4*(fun(t(n,1),x(n,:))+fun(t(n,1)+k/2,xguess))-xguess';

                try
                    g = k/4*Jac(t(n,1)+k/2,xguess)-eye(length(xguess));
                catch
                    g = k/4*jac(fun,t(n,1)+k/2,xguess)-eye(length(xguess));
                end

                delta = g\f;

                if max(abs(delta)) < TOL2
                    fFound = 1;
                    xmid = xguess - delta';
                else
                    xguess = xguess - delta';
                end
            end

            if j >= MaxIt
                error('Max Iterations on Solver has been Exceeded')
            end

            x1 = 2*xmid - x(n,:);

            fFound = 0;
            j = 0;
            xguess = x(n,:);

            while ~fFound && j < MaxIt

                j = j + 1;
                f = k*fun(t(n,1)+k,xguess)-x(n,:)'+4*xmid'-3*xguess';

                try
                    g = k*Jac(t(n,1)+k,xguess) - 3*eye(length(xguess));
                catch
                    g = k*jac(fun,t(n,1)+k,xguess)-3*eye(length(xguess));
                end

                delta = g\f;

                if max(abs(delta)) < TOL2
                    fFound = 1;
                    x2 = xguess - delta';
                else
                    xguess = xguess - delta';
                end
            end

            if j >= MaxIt
                error('Max Iterations on Solver has been Exceeded')
            end

            tau = norm(x2-x1,2);
            if tau < TOL1
                x(n+1,:) = x2;
                t(n+1,:) = t(n,1) + k;
                if exitCon && exCond(x2)
                    if conTry > 25
                        done = 1;
                        exitCon = 0;
                    else
                        conTry = conTry + 1;
                        n = n - 1;
                        tau = 100*TOL1;
                    end
                else
                    conTry = 0;
                end                
                n = n + 1;
            else
                done = 0;
            end
            k = min(kmax,max(kmin,0.95*k*(TOL1/(tau+eps()))^(1/2)));        
        end
    end
    if exitCon
        warning('Solver did not reach the exit condition...')
    end
end

function jac = jac(fun,t,x)
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
end