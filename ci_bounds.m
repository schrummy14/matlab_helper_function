function [ci_low, ci_high] = ci_bounds(data,alpha)

    if nargin < 2
        alpha = 0.05;
    end

    N = size(data,1);
    std_data = std(data);
    
    ci = std_data/sqrt(N) * tinv(1-alpha/2,N-1);
    
    if nargout < 2
        ci_low = ci;
    else
        mean_data = mean(data);
        ci_low = mean_data - ci;
        ci_high = mean_data + ci;
    end
    
end
    