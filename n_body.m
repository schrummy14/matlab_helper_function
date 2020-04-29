function df = n_body(~,x,mass,G,dim)


    N = length(mass);
    
    F = zeros(N,N,dim);
    for m = 1:N
        for n = m:N
            if m~=n
                xm = x((m-1)*dim+1:m*dim);
                xn = x((n-1)*dim+1:n*dim);
                num = G*(xm-xn);
                den = norm(xm-xn,2).^3;
                F(n,m,:) = num./den;
            end
        end
    end
    
    for m = 1:dim
        F(:,:,m) = F(:,:,m)-F(:,:,m)';
        F(:,:,m) = F(:,:,m).*repmat(mass,N,1);
    end
    
    F = sum(F,2);
    F = F(:);
    F = reshape(F,N,dim);
    F = F';
    F = F(:);
    
    df = zeros(size(x));
    
    df(1:N*dim) = x(N*dim+1:end);
    df(N*dim+1:end) = F;
    
end