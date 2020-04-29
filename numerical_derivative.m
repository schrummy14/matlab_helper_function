function J = numerical_derivative(fun,x)

    n = length(x);
    fx = fun(x);
    xn = x;
    J = zeros(n);
    dx = 1e-6*norm(x);
    for m = 1:n
        xn(m) = xn(m) + dx;
        J(:,m) = (fun(xn)-fx)/dx;
        xn(m) = x(m);
    end

end