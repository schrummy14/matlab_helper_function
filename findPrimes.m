function p_list = FindPrimes(n,p_list)

    
    p_start = 1;
    xx = 0;
    if nargin < 2
        p_list = zeros(floor(sqrt(n)),1);
        if nargin < 1
            error 'Function must have at least one input'
        end
    else
        if length(p_list) < n
            gg = n-length(p_list);
            p_start = p_list(end)+1;
            p_list(end+1:end+gg,:) = 0;
            xx = length(p_list(p_list>0));
        elseif length(p_list) == floor(sqrt(n))
            return
        end
    end
    for m = p_start:n
        if m > 2
            found_prime = true;
            for j = 1:xx
                if mod(m,p_list(j)) == 0
                    found_prime = false;
                    break
                end
            end
            if found_prime
                p_list(xx+1) = m;
                xx = xx + 1;
            end
        else
            if m > 1
                p_list(xx+1) = 2;
                xx = xx + 1;
            end
        end
    end
            
    p_list = p_list(p_list > 0);
    
end