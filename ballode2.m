function ballode2
%BALLODE  Run a demo of a bouncing ball.
%   This is an example of repeated event location, where the initial
%   conditions are changed after each terminal event.  This demo computes ten
%   bounces with calls to ODE23.  The speed of the ball is attenuated by 0.9
%   after each bounce. The trajectory is plotted using the output function
%   ODEPLOT.
%
%   See also ODE23, ODE45, ODE56 ODESET, ODEPLOT, FUNCTION_HANDLE.

%   Mark W. Reichelt and Lawrence F. Shampine, 1/3/95
%   Copyright 1984-2014 The MathWorks, Inc.

tstart = 0;
tfinal = 30;
y0 = [0; 20];
refine = 10;
RelTol = 1e-3;
AbsTol = 1e-4;
options = odeset('Events',@events,'OutputFcn',@odeplot,'OutputSel',1,...
   'Refine',refine,'RelTol',RelTol,'AbsTol',AbsTol);

fig = figure;
ax = axes;
ax.XLim = [0 30];
ax.YLim = [0 25];
box on
hold on;

tout = tstart;
yout = y0.';
teout = [];
yeout = [];
ieout = [];
for i = 1:10
   % Solve until the first terminal event.
   [t,y,te,ye,ie] = ode56(@f,[tstart tfinal],y0,options);
   if ~ishold
      hold on
   end
   % Accumulate output.  This could be passed out as output arguments.
   nt = length(t);
   tout = [tout; t(2:nt)];
   yout = [yout; y(2:nt,:)];
   teout = [teout; te];          % Events at tstart are never reported.
   yeout = [yeout; ye];
   ieout = [ieout; ie];
   
   ud = fig.UserData;
   if ud.stop
      break;
   end
   
   % Set the new initial conditions, with .9 attenuation.
   y0(1) = 0;
   y0(2) = -.9*y(nt,2);
   
   % A good guess of a valid first timestep is the length of the last valid
   % timestep, so use it for faster computation.  'refine' is 4 by default.
   options = odeset(options,'InitialStep',t(nt)-t(nt-refine),...
      'MaxStep',t(nt)-t(1),'Refine',refine,'RelTol',RelTol,'AbsTol',AbsTol);
   
   tstart = t(nt);
end

plot(teout,yeout(:,1),'ro')
xlabel('time');
ylabel('height');
title('Ball trajectory and the events');
hold off
odeplot([],[],'done');

% --------------------------------------------------------------------------

function dydt = f(~,y)
dy1 = y(2);
dy2 = -9.8;%-.1*y(1) - .1*abs(y(1)).*y(1);
dydt = [
    dy1
    dy2
    ];

% --------------------------------------------------------------------------

function [value,isterminal,direction] = events(~,y)
% Locate the time when height passes through zero in a decreasing direction
% and stop integration.
value = y(1);     % detect height = 0
isterminal = 1;   % stop the integration
direction = -1;   % negative direction