function plotFigureSpeedRun
 x = linspace(0,2*pi,1001);
 y1= cos(x);
 y2= sin(x);
 y3= y1-y2;
 y4= y2-y1;
 y5= y1+y2;
 y6= y1.*y2;
 y = [y1;y2;y3;y4;y5;y6];
 numpoints = 31;
 Legends = {
 'Thanks'
 'for'
 'looking'
 'at'
 'my'
 'function!'};
 plotFigure(x,y,numpoints,Legends)
end