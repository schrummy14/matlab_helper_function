function [t,x] = RungeKutta4(fun,tspan,x0,N)

    if size(x0,1) ~= 1
        x0 = x0';
    end
    if numel(tspan) == 1, tspan(2) = tspan; tspan(1) = 0; end

    k = diff(tspan)/N;

    t(1,1) = tspan(1);
    x(1,:) = x0;

    for n = 1:N
        t(n+1,1) = tspan(1) + n*k;

        k1 = fun(t(n,1),x(n,:));

        y = x(n,:) + k/2 * k1';
        k2 = fun(t(n,1)+k/2,y);

        y = x(n,:) + k/2 * k2';
        k3 = fun(t(n,1)+k/2,y);

        y = x(n,:) + k*k3';
        k4 = fun(t(n,1)+k,y);

        x(n+1,:) = x(n,:) + k/6 * (k1 + 2*k2 + 2*k3 + k4)';
    end
    
end