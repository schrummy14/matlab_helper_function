function [meanData, CI_Data] = meanCI(data, alpha)

    if min(size(data)) > 1
        error 'Data must be a vector'
    end
    
    if nargin < 2
        alpha = 0.05;
    end
    
    data(isnan(data)) = [];
    
    alpha = alpha/2;

    m = mean(data);
    s = std(data);
    N = length(data);
    
    val = tinv(1-alpha,N)*s/sqrt(N);
    
    if nargout > 1
        meanData = m;
        CI_Data = val;
    else
        meanData = m + [-1,1]*val;
    end
    
end