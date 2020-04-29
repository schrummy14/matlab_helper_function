function [avex, turx, a2tx] = average(x,n)
%%
% This function will calculate the average value of x every n points along
% with the difference from the mean for every entry in x.
% If n is not entered, avex = average(x) = mean(x), and turx = mean(x)-x.
% a2tx are the x points in which the turbulence was calculated for.

    if nargin == 1,n = length(x);end

    if n == 0
        error('Can not average zero entries...')
    end

    if n > 1
        s1 = size(x,1);
        M = s1 - mod(s1,n);
        y = reshape(x(1:M),n,[]);
        avex = transpose(sum(y,1)/n);
    else
        avex = x;
    end

    if nargout>=2
        
        k = n;
        for i = length(avex):-1:1
            for j = n:-1:1
                turx((i-1)*k+j) = avex(i)-x((i-1)*k+j);
            end
        end
        
        if nargout>=3

            a2tx(:,1) = x(1:length(turx));
            
        end
    end
    
end

