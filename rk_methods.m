clc 
clear

c = -0.5;
fun = @(t,x)c.*x;

x0 = 5;
t0 = 0;
tN = 10;
N = 1000;

dt = (tN-t0)/N;

x_euler = zeros(N+1,1);
x_buler = zeros(N+1,1);
x_mod_e = zeros(N+1,1);
x_rk3   = zeros(N+1,1);
x_rk4   = zeros(N+1,1);
t       = zeros(N+1,1);

x_euler(1) = x0;
x_buler(1) = x0;
x_mod_e(1) = x0;
x_rk3(1)   = x0;
x_rk4(1)   = x0;
t(1)       = t0;

for m = 1:N
    
    % Euler Method
    x_euler(m+1) = x_euler(m) + dt * fun(t(m),x_euler(m));
    
    % Backward Euler (Linear Function...) 
    % x_(n+1) = x_(n) + dt*fun(t_(n+1),x_(n+1))
    % x_(n+1) = x_(n) + dt*(c*x_(n+1))
    % (x_(n+1)-dt*c*x_(n+1)) = x_(n)
    % x_(n+1)*(1-dt*c) = x_(n)
    % x_buler(m+1) = x_buler(m)/(1-dt*c);
    
    % Mod Euler (Heun's Method)
    x_mid = x_mod_e(m) + dt * fun(t(m),x_mod_e(m));
    x_mod_e(m+1) = x_mod_e(m) + dt/2 * (fun(t(m),x_mod_e(m)) + fun(t(m)+dt,x_mid));
    
    % Kutta third order
    k1 = fun(t(m),x_rk3(m));
    k2 = fun(t(m)+0.5*dt,x_rk3(m)+0.5*dt*k1);
    k3 = fun(t(m)+1.0*dt,x_rk3(m)+dt*(2*k2-k1));
    x_rk3(m+1) = x_rk3(m) + dt/6 * (k1 + 4*k2 + k3);
    
    % The Runge Kutta Method
    k1 = fun(t(m),x_rk4(m));
    k2 = fun(t(m)+0.5*dt,x_rk4(m)+dt/2 * k1);
    k3 = fun(t(m)+0.5*dt,x_rk4(m)+dt/2 * k2);
    k4 = fun(t(m)+1.0*dt,x_rk4(m)+dt   * k3);
    x_rk4(m+1) = x_rk4(m) + dt/6 * (k1 + 2*k2 + 2*k3 + k4);
    
    t(m+1) = t(m) + dt;
    
end

x_true = x0*exp(c*t);

h = semilogy(t,abs(repmat(x_true,1,4) - [x_euler,x_mod_e,x_rk3,x_rk4]));
% h = plot(t,[x_euler,x_mod_e,x_rk3,x_rk4]);
for m = 1:4
    h(m).LineWidth = 2.5;
end
xlim([0 10])
ylim([1e-20 1])