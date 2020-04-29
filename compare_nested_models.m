function [h, Ftest, Fcrit] = compare_nested_models(SSEr, SSEf, n, p, q, alpha)
    
    if (isa(SSEr,'LinearModel') || isa(SSEr,'NonLinearModel')) ...
            && (isa(SSEf, 'LinearModel') || isa(SSEf, 'NonLinearModel'))
        
        if nargin < 3
            n = 0.05;
        end
        
        nr = SSEr.NumObservations;
        nf = SSEf.NumObservations;
        
        alpha = n;
        
        if nr ~= nf
            error('Linear Models must come from the same data set')
        end
        
        n = nr;
        p = n - SSEf.DFE;
        q = n - SSEr.DFE;
        
        SSEr = SSEr.SSE;
        SSEf = SSEf.SSE;
        
    else

        if nargin < 6
            alpha = 0.10;
        end
        
    end

    if isempty(SSEr) && isempty(SSEf)
        if nargout > 1
            error 'Can only get the Fcrit if SSEs are not given'
        end
        h = finv(1-(alpha),n-p,p-q);
        return;
    else
        MSE = SSEf./(n-p);
        Ftest = abs(((SSEr - SSEf)./(p-q))./MSE);
    end
    
    Fcrit = finv(1-(alpha),n-p,p-q);
    h = Ftest > Fcrit;
    
end 