% Lets build a grid generator....
clear
clc
clf
hold on

xb = [0.0 0.4 0.6 1.0 1.0 0.0 0.0];
yb = [0.0 0.0 0.2 0.2 1.0 1.0 0.0];

plot(xb,yb,'-b')

grid_x_start = linspace(min(xb),max(xb),31);
grid_y_start = linspace(min(yb),max(yb),31);
[X,Y] = meshgrid(grid_x_start,grid_y_start);

plot(X,Y,'ro')




hold off