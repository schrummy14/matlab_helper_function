function y = msmooth(x,n)
    n = abs(n);
    if abs(round(n)-n)>0
        error('n must be an odd integer')
    end
    if n/2 == floor(n/2)
        error('n must an odd integer')
    end
    y = zeros(size(x));
    n2 = floor(n/2);
    if n == 1
        y = x;
        return
    end 
    for m = 1:length(x)
        if m <= n2
            y(m) = meanFind(x(m:n+m-1));
        elseif m >= length(x)-n2
            y(m) = meanFind(x(m-n+1:m));
        else
            y(m) = meanFind(x(m-n2:m+n2));
        end
    end
end
function y = meanFind(x)
    g = sort(x);
    g = g(~isnan(g));
    y = mean(g(2:end-1));
end