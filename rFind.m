function [x,nTry] = rFind(xn,n)
    if xn<0 && iseven(n)
        xguess = 1i;
    else
        xguess = 1;
    end
    done = false;
    nTry = 0;
    oldRelErr = inf;
    while ~done      
        
        x = 1/n .* ((n-1).*xguess + xn./xguess.^(n-1)); 
        
        RelErr = abs((x - xguess)./x);
        if RelErr == 0
            done = true;
        elseif oldRelErr <= RelErr && nTry > 5
            done = true;
            x = xguess;
        else
            oldRelErr = RelErr;
            xguess = x;
        end
        if nTry > 9 && RelErr<=eps
            done = true;
        end
        nTry = nTry + 1;
    end
    if abs(imag(x))>0 && RelErr>0
        rx = real(x);
        ix = imag(x);
        if abs(rx) < eps(abs(x))
            rx = 0;
        end
        if abs(ix) < eps(abs(x))
            ix = 0;
        end
        x = rx+ix*1i;
    end
end

function f = iseven(x)
    if x/2==floor(x/2)
        f = true;
    else
        f = false;
    end
end