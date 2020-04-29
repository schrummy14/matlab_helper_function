clc
clear
x = linspace(0,1,25);
y = x.^3 + 1/10.*rand(1,length(x));
y = log(y+1);
order = 3;
[ynew,rsq,beta] = LinReg(x,y,order);
plot(x,ynew,x,y,'o')
str = sprintf('Regression of order %i. R^2 value is %.2f%%.',order,rsq*100);
title(str)