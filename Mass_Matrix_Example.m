function Mass_Matrix_Example(N)
%FEM1ODE  Stiff problem with a time-dependent mass matrix 

if nargin < 1
  N = 5001;
end
h = pi/(N+1);
y0 = 1*sin(h*(1:N)');
tspan = linspace(0, pi,N).';

% The Jacobian is constant.
e = repmat(1/h,N,1);    %  e=[(1/h) ... (1/h)];
d = repmat(-2/h,N,1);   %  d=[(-2/h) ... (-2/h)]; 
J = spdiags([e d e], -1:1, N, N);

options = odeset('Mass',@mass,'MStateDependence','none', ...
                 'Jacobian',J);

Sol = ode15s(@f,tspan,y0,options,N);
if strcmp(Sol.solver,'ode56')
    y = Sol.y.';
    tspan = Sol.x.';
else
    y = deval(Sol,tspan).';
end
disp(Sol.stats)

surf((1:N)/(N+1),tspan,y,'LineStyle','none');
set(gca,'ZLim',[0 1]);
view(142.5,30);
title(['Finite element problem with time-dependent mass ' ...
       'matrix, solved by ',Sol.solver]);
xlabel('space ( x/\pi )');
ylabel('time');
zlabel('solution');
%---------------------------------------------------------------
function out = f(~,y,N)
h = pi/(N+1);
e = repmat(1/h,N,1);    %  e=[(1/h) ... (1/h)];
d = repmat(-2/h,N,1);   %  d=[(-2/h) ... (-2/h)]; 
J = spdiags([e d e], -1:1, N, N);
out = J*y;
%---------------------------------------------------------------
function M = mass(t,N)
h = pi/(N+1);
e = repmat(exp(-t)*h/6,N,1);  % e(i)=exp(-t)*h/6
e4 = repmat(4*exp(-t)*h/6,N,1); 
M = spdiags([e e4 e], -1:1, N, N);