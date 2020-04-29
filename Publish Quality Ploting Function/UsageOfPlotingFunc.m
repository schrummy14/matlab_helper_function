clear; clc; close all

PlotData = subplot_class;

% X = 0:0.1:8;
% Y = sin(X);
% for i = 1:20
%     for j = 1:i
%         PlotData(i).X(:,j)= X;
%         PlotData(i).Y(:,j) = exp(j)*Y;
%     end
%     PlotData(i).Axis = [0 8 -(exp(i)+1) exp(i)+1 8 6 1 1];
%     PlotData(i).XLabel = ['\makebox[4in][c]{\textbf {X(',num2str(i),')}}'];
%     PlotData(i).YLabel = ['\makebox[4in][c]{\textbf {Y(',num2str(i),')}}'];
% end
% 
% Legend_Names = {'Sim 1' 'Sim 2' 'Sim 3' 'Sim 4' 'Sim 5' 'Sim 6'...
%                 'Sim 7' 'Sim 8' 'Sim 9' 'Sim 10' 'Sim 11' 'Sim 12'};
% multiple_LinePlot_func(PlotData,'PerRow',[1 1 1],...
%     'CommonLegend',Legend_Names,'CommonLegendLocation','northeast',...
%     'CommonTitle','','Size',1,'FigureHeight',10,'Print','Black' ,...
%     'MarkerFrequency',10,'WhiteSpace',[0.15 0.05 0.1 0.05],'Margin',[0.1 0.1 0.1]);


X = linspace(0,1,101).';
Y1 = X;
Y2 = X.^2;
Y3 = X.^3;
Y4 = X.^10;
for k = 1:3
    PlotData(1).X(:,k) = X;
end
PlotData(2).X = X;

PlotData(1).Y(:,1) = Y1;
PlotData(1).Y(:,2) = Y2;
PlotData(1).Y(:,3) = Y3;
PlotData(2).Y(:,1) = Y4;

for k = 1:2
    PlotData(k).Axis = [0 1 0 1 8 6 1 1];
end

Legend_Names = {'My Sim 1', 'Another Sim', 'Because I wanted To',...
                'One More Time!!!'};

multiple_LinePlot_func(PlotData,'Print','Color','MarkerFrequency',8,...
                                'CommonLegend',Legend_Names,...
                                'CommonLegendLocation','northeast',...
                                'PerRow',[1 1])


