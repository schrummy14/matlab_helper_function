clear;

u0=[0;1;-1];

tspan = [0 70];

OPTIONS = odeset('RelTol',1e-10,'AbsTol',1e-14);
[t,u] = ode56(@lorenz,tspan,u0,OPTIONS);
% OPTIONS = {'TOL',5e-16};
% [t,u] = rk78(@lorenz,tspan,u0,OPTIONS{:});

figure;
clf;
tmp = 'xyz';
for k=1:3
    subplot(3,1,k)
    p1=plot(t,u(:,k),'r-');
    set(p1,'linewidth',0.5);
    set(gca,'fontsize',16);
    if k==1
        t1=title('Lorenz equation');
        set(t1,'fontsize',16);
    end
    t2=xlabel('t');
    set(t2,'fontsize',16);
    t3=ylabel([tmp(k),'(t)']);
    set(t3,'fontsize',16);
    set(gca,'xtick',0:10:70);
    if k==1
        set(gca,'ytick',-20:10:20);
    end
    if k==2
        set(gca,'ytick',-50:25:50);
    end
    if k==3
        set(gca,'ytick',0:15:60);
    end
end

figure
clf;
p2=plot3(u(:,1),u(:,2),u(:,3),'r-');
set(p2,'linewidth',0.5);
view([50 30]);
set(gca,'fontsize',16);
axis([-20 20 -30 30 0 50]);
set(gca,'xtick',-50:10:50);
set(gca,'ytick',-50:10:50);
set(gca,'ztick',-50:10:50);
grid on;
s1=xlabel('x');
s2=ylabel('y');
s3=zlabel('z');
set(s1,'fontsize',16);
set(s2,'fontsize',16);
set(s3,'fontsize',16);
t1=title('Lorenz equation example');
set(t1,'fontsize',16);