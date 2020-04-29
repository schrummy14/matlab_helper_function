function [x,tau,n] = Jacobi(A,b,varargin)

    x0 = zeros(size(A,1),1);
    TOL = 0;
    maxIt = 5000;
    w = 1;
    for n = 1:2:length(varargin)
        switch varargin{n}
            case 'x0'
                x0 = varargin{n+1};
            case 'TOL'
                TOL = varargin{n+1};
            case 'maxIt'
                maxIt = varargin{n+1};
            case 'weight'
                w = varargin{n+1};
            otherwise
                error('Bad Argument Passed')
        end
    end    
    
    diag_val = abs(diag(A));
    row_sum = sum(abs(A-diag(diag(A))),2);
     if any(row_sum>diag_val)
         error('Jacobi Method will not converge, please try a different method...')
     end
    
    d = diag(diag(A));
    upT = d - triu(A);
    loT = d - tril(A);
    R = upT + loT;
    invd = diag(1./diag(A));
    
    done = 0;
    n = 0;
    xold = x0;
    while ~done
        xnew = w*invd*R*xold + w*invd*b +(1-w)*xold;
        n = n+1;
        tau = norm(A*xnew - b,2);
        tauDif = norm(xnew-xold,inf);
        if tau<TOL || n>=maxIt || tauDif < eps
            done = 1;
            x = xnew;
            if TOL<tau && TOL > 0
                msg1 = 'Solver has stopped before Tol value was reached. ';
                if n >maxIt
                    msg2 = 'Max Iterations allowed has been reached.';
                elseif tauDif < eps
                    msg2 = 'The difference between the old and new iteration values is less than eps.';
                else
                    msg2 = 'Something else happened...';
                end
                warning([msg1,msg2]);
            end
        else
            xold = xnew;
        end
    end
end


