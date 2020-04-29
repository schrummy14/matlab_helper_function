function varargout = scaledEuler(fun,tspan,x0)

    if nargout > 2
        error('Output can only handle 2 outputs at max')
    end
    h_old = 16*sqrt(eps);
    TOL = 1e-5;
    M = eye(length(x0));

    lamba = 1.4;
    alpha = 0.95;
    rho = zeros(length(x0),1);

    n = 1;
    t(n,1) = tspan(1);
    x(n,:) = x0;
    nStep = 0;
    fEval = 0;
    nFail = 0;

    h = 2*lamba*h_old;
    done = false;
    while ~done

        if t(n,1) + h > tspan(2)
            h = tspan(2) - t(n,1);
            done = true;
        end

        x11 = tryStep22(fun,t(n,1),x(n,:),M,h);
        x22 = tryStep12(fun,t(n,1),x(n,:),M,h);
        fEval = fEval + 3;
        yError = x22-x11;
        hprime = h/(2*norm(yError,inf)/TOL)^(1/2);

        if h > 10*hprime
            h = 2*hprime;
            nFail = nFail + 1;
        else

            mprime = M*lamba;

            xm11 = tryStep22(fun,t(n,1),x(n,:),mprime,h);
            xm22 = tryStep12(fun,t(n,1),x(n,:),mprime,h);
            mError = xm22 - xm11;
            fEval = fEval + 3;

            for m = 1:length(rho)
                rho(m) = (h^2*alpha^2*M(m,m) + h*alpha*M(m,m) - 1 + alpha - h + h*alpha^2)/...
                    (h*alpha*M(m,m)*(1+h));

                if abs(mError(m)) < abs(yError(m))
                    M(m,m) = M(m,m)*lamba;
                elseif abs(mError(m)) > abs(yError(m))
                    M(m,m) = max(1,M(m,m)*rho(m));
                end

            end

            x(n+1,:) = x11;
            t(n+1,1) = t(n,1) + h;
            n = n + 1;
            h = 2*lamba*h;
            nStep = nStep + 1;
        end
    end
    if nargout > 1
        varargout{1} = t;
        varargout{2} = x;
    else
        Sol.t = t;
        Sol.x = x;
        Sol.fEvals = fEval;
        Sol.nSteps = nStep;
        Sol.nFail  = nFail;
        varargout{1} = Sol;
    end
end

function x = tryStep12(fun,t,x,M,h)
    A = eye(length(x))+h/2*M;
    x1 = x + (h/2*(1+h/2)*(A\fun(t,x))).';
    x  = x1+ (h/2*(1+h/2)*(A\fun(t+h/2,x1))).';    
end

function x = tryStep22(fun,t,x0,M,h)
    A = eye(length(x0))+h*M;
    x = x0 + (h*(1+h)*(A\fun(t,x0))).';
end   
  