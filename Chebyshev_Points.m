function x = Chebyshev_Points(a,b,n)

    if nargin<3
        n = 100;
    end
    
    k = n:-1:1;
    x = 0.5.*(a+b)+0.5.*(b-a).*cos((k-0.5)./n.*pi);
    
end