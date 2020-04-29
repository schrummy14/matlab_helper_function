clc
clear
clf

% Global Characteristics
g = 9.81;

% Characteristics of the first pendulum
l1 = 1;
m1 = 1;
T1_0 = pi;
dT1_0= 0;

% Characteristics of the second pendulum
l2 = 0.5;
m2 = 1;
T2_0 = pi/20;
dT2_0= 0;

a = (m1+m2)*l1;
b = @(x)m2*l2*cos(x(1)-x(3));
c = @(x)m2*l1*cos(x(1)-x(3));
d = m2*l2;
e = @(x)-m2*l2*x(4)* x(4)*sin(x(1)-x(3))-g*(m1+m2)*sin(x(1));
f = @(x)m2*l1*x(2)*x(2)*sin(x(1)-x(3))-m2*g*sin(x(3)) ;

fun = @(t,x)[
    x(2)
    
    (e(x)*d-b(x)*f(x))/(a*d-c(x)*b(x))
    
    x(4)
    
    (a*f(x)-c(x)*e(x))/(a*d-c(x)*b(x))
    ];

% Opts = {'refine',5,'Rel_TOL',1e-8,'Abs_TOL',1e-12,'Do_Waitbar',true};
% [t,deg] = TR_BDF2(fun,[0 60],[T1_0,dT1_0,T2_0,dT2_0],Opts);
opts = odeset('RelTol',1e-6,'AbsTol',1e-10);
[t,deg] = ode56(fun,[0 10],[T1_0,dT1_0,T2_0,dT2_0],opts);
%%
deg_1 = deg(:,1);
deg_2 = deg(:,3);

x1 = l1*sin(deg_1);
y1 = -l1*cos(deg_1);

x2 = x1 + l2*sin(deg_2);
y2 = y1 - l2*cos(deg_2);

s = 0.02;
t_old = 0;
tracker_1 = [];
tracker_2 = [];

axis_val = 1.1*[-(l1+l2), (l1+l2), -(l1+l2), (l1+l2)];
x1_1 = 0;
y1_1 = 0;
for k = 1:length(x1)
    if k == 1 || t(k)>t_old+s || k == length(x1)
        clf
        x1_2 = x1(k);
        y1_2 = y1(k);

        x2_1 = x1(k);
        y2_1 = y1(k);
        x2_2 = x2(k);
        y2_2 = y2(k);
        
        tracker_1(end+1,:) = [x1_2 y1_2];
        tracker_2(end+1,:) = [x2_2 y2_2];

        plot( ...
            [x1_1 x1_2],[y1_1 y1_2],'b-', ...
            [x2_1 x2_2],[y2_1 y2_2],'r-',...
            tracker_1(:,1),tracker_1(:,2),'b.',...
            tracker_2(:,1),tracker_2(:,2),'r.')
        title(['Time = ',num2str(t(k)),' seconds'])
        axis(axis_val)
        drawnow
        t_old = t_old + s;
    end
end
    

% plot(x1,y1,x2,y2)
