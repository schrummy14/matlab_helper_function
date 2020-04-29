function varargout = bvpSol(fun,tspan,bc,varargin)

    guess = zeros(1,length(bc)/2);
    TOL = 1e-8;
    maxEval = 500;
    refine = 'auto';
    derMeth = 'forward';
    
    for m = 1:2:length(varargin)
        switch varargin{m}
            case 'guess'
                guess = varargin{m+1};
                if size(guess,1)>1 % Make sure its a row vector
                    guess = guess.';
                end
                if length(guess) ~= length(bc)/2 % Make sure its the right size
                    error('Guess vector must be the same size as half of the BC...')
                end
            case 'TOL'
                TOL = varargin{m+1};
            case 'maxEval'
                maxEval = varargin{m+1};
            case 'Refine'
                refine = varargin{m+1};
            case 'derMeth'
                derMeth = varargin{m+1};
            otherwise
                warning([varargin{m},' is not a valid argument...skipping...'])
        end
    end
    
    gg = @(s)solFind(fun,tspan,bc,s,TOL);
    [ss,res,fevals] = zeroFind(gg,guess,TOL,maxEval,derMeth);
    SOL = rk65(fun,tspan,[bc(1:length(bc)/2),ss],'TOL',TOL,'Refine',refine);
    fevals = fevals + 1;
    res = abs(res);
   
    if nargout < 2
        if nargout < 1
            plot(SOL.t,SOL.x,'-o')
        end
        sol.t = SOL.t;
        sol.x = SOL.x;
        sol.ss = ss;
        sol.res = res;
        sol.fevals = fevals*SOL.stats.nfevals;
        varargout{1} = sol;
    elseif nargout == 2
        varargout{1} = SOL.t;
        varargout{2} = SOL.x;
    end
    
end

function res = solFind(fun,tspan,bc,s,TOL)
    [~,x] = rk65(fun,tspan,[bc(1:length(bc)/2),s],'TOL',TOL,'Refine',1);
    res =(bc(length(bc)/2+1:end)-x(end,1:length(bc)/2));
end

function [ss,res,fevals] = zeroFind(gg,sguess,TOL,maxEval,derMeth)
    
    xold = sguess;    
    fold = gg(xold).';
    fevals = 1;
    done = 0;
    if TOL == 0
        warning('TOL can not be zero... using machine epsilon insted...')
        if TOL == 0
            TOL = eps(numel(sguess));
        end
    end
    m = 0;
    while ~done
        jac = nJacZeroFind(gg,xold,derMeth);
        xnew = xold - (jac\fold).';
        fnew = gg(xnew).';
        fevals = fevals + 2 + length(xold);
        m = m + 1;
        if norm(fnew)<TOL(1) || norm(xnew-xold)<eps || fevals > maxEval
            done = 1;
            if fevals > maxEval
                warning('Solution was not found given the supplied guess...')
            end
%         elseif norm(fnew)>norm(fold) && m > 10
%             done = 1;
%             xnew = xold;
%             fnew = fold;
        else
            xold = xnew;
            fold = fnew;
        end
    end
    res = norm(fnew);
    ss = xnew;
end

function Jac = nJacZeroFind(fun,xold,derMeth)
    switch derMeth
        case 'forward'
            n = length(xold);
            fold = fun(xold).';
            xoldn = xold;
            Jac = zeros(n);
            step = 1e-8;
            for m = 1:n
                xoldn(m) = xoldn(m) + step;
                Jac(:,m) = (fun(xoldn).'-fold)/step;
                xoldn(m) = xold(m);
            end
        case 'backward'
            n = length(xold);
            fold = fun(xold).';
            xoldn = xold;
            Jac = zeros(n);
            step = 1e-8;
            for m = 1:n
                xoldn(m) = xoldn(m) - step;
                Jac(:,m) = (fold - fun(xoldn).')/step;
                xoldn(m) = xold(m);
            end
        case 'central'
            n = length(xold);
            fold = fun(xold).';
            xoldn = xold;
            Jac1 = zeros(n);
            step = 2e-8;
            for m = 1:n
                xoldn(m) = xoldn(m) + step;
                Jac1(:,m) = (fun(xoldn).'-fold)/step;
                xoldn(m) = xold(m);
            end
            n = length(xold);
            fold = fun(xold).';
            xoldn = xold;
            Jac2 = zeros(n);
            for m = 1:n
                xoldn(m) = xoldn(m) - step;
                Jac2(:,m) = (fold - fun(xoldn).')/step;
                xoldn(m) = xold(m);
            end
            Jac = 1/2*(Jac1+Jac2);            
    end
            
end