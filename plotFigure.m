% Function plotFigure(x,y,numpoints,Legends,color,LineStyle,Markers)
% This function is used to easily plot multiple datasets quickly. It
% automaticaly sets different linestyles and marker styles for each
% dataset, upto 33 of them. It does this in a black color that is ready for
% input into any paper. It even correctly labels each dataset with its own
% legend entry.
%
% plotFigure(x,y,numpoints,Legends)
% 
% where
% x contains the span of the data. This must be a vector.
%
% y contains the range of the data and can contain upto 33 different
% datasets
%
% numpoints is the number of markes used for the plot. This can be either a
% single value that will be used for all plots, or the user can specify how
% many markers are used per plot. (Default value is 21)
%
% Legends will be the name given on a plot per plot basis. The user will
% need to name each plot seperetly if they choose to use a legend. Legend
% needs to be entered as a class type. Legend = {'mydata1','mydata2',...}.
% The exception to this if the user inputs the keywords 'default', or
% 'none'. These may be entered as text. The default setting will give the
% following plot names (data 1, data 2, data 3, ...). The none setting will
% simply not show a legend.
%
% Example
% x = linspace(0,2*pi,1001);
% y1= cos(x);
% y2= sin(x);
% y3= y1-y2;
% y4= y2-y1;
% y5= y1+y2;
% y6= y1.*y2;
% y = [y1;y2;y3;y4;y5;y6];
% numpoints = 31;
% Legends = {
% 'Thanks'
% 'for'
% 'looking'
% 'at'
% 'my'
% 'function!'};
% plotFigure(x,y,numpoints,Legends)


% My info
% Function written by 
% Matt Schramm 
% of Iowa State University 
% on June 12, 2016

% This 


%%%%%% Main Function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotFigure(x,y,numpoints,Legends,color,LineStyles,Markers)

    if nargin < 7
        if nargin == 6
            error 'Need both LineStyle and Markers'
        else
            LineStyles = 'default';
            Markers   = 'default';
            
            if nargin < 5
                color = {'default'};
                if nargin < 4
                    Legends = {'default'};
                    if nargin < 3
                        numpoints = 21;
                    end
                end
            end
        end
    end
    
    % Check input data
    [x,y,numpoints,num_data_sets,mX,mY,LegStr,color] = ...
        check_and_update_data(x,y,numpoints,Legends,color);
    
    % Generates all combinations of linestyles and markers needed (max 33)
    [LineStyles,Markers] = getStyles(num_data_sets,LineStyles,Markers);
    
    % Clear open figure, if any
    clf
    hold on
    for k = num_data_sets:-1:1
        SkipLines = round(linspace(1,length(x),numpoints(k)));
        LS = linestyle(k,LineStyles,color);
        plot(x,y(:,k),LS{:});
        MS = markerstyle(k,LineStyles,Markers,color);
        plot(x(SkipLines),y(SkipLines,k),MS{:});
        LegMark = getLegendMarker(LS{2},MS{2});
        h(k) = plot(mX,mY,LegMark,'Color',color{k},'visible','off');
    end
    hold off
    if ~strcmp(LegStr{1},'none')
        h_leg = legend(h,LegStr{:},'Location','Eastoutside');
        set(h_leg,'FontSize',14);
    end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% Helper Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function LegMark = getLegendMarker(LS,MS)

    if strcmp(MS,'none')
        LegMark = LS;
    else
        LegMark = [LS,MS];
    end

end

function [LS,MS] = getStyles(num_data_sets,LineStyles,Markers)

    if strcmp(LineStyles,'default')
        posLS = getLS;
        LS = cell(num_data_sets,1);
        for mn = 1:num_data_sets
           LS{mn} = posLS{mod(mn-1,length(posLS))+1};
        end
    else
        if iscell(LineStyles)
            if strcmp(LineStyles{1},'default')
                posLS = getLS;
                LS = cell(num_data_sets,1);
                for mn = 1:num_data_sets
                   LS{mn} = posLS{mod(mn-1,length(posLS))+1};
                end
            elseif length(LineStyles)<num_data_sets
                error 'There must be a LineStyle for each dataset...'
            else
                LS = LineStyles;
            end
        end
    end
    
    if strcmp(Markers,'default')
        posMS = getMS;
        MS = cell(num_data_sets,1);
        for mn = 1:num_data_sets
           MS{mn} = posMS{mod(mn-1,length(posMS))+1};
        end
    else
        if iscell(Markers)
            if strcmp(Markers{1},'default')
                posMS = getMS;
                MS = cell(num_data_sets,1);
                for mn = 1:num_data_sets
                   MS{mn} = posMS{mod(mn-1,length(posMS))+1};
                end
            elseif length(Markers)<num_data_sets
                error 'There must be a Marker for each dataset...'
            else
                MS = Markers;
            end
        end
    end

end

function out=getLS
    out = {'-','--','-.'};
end

function out=getMS
    out = {'none','o','*','x','s','d','v','^','<','>','pentagram'};
end

function out = linestyle(k,LS,color)
    out = {'LineStyle',LS{k},'Color',color{k}};
end

function out = markerstyle(k,LS,MS,color)

    out1 = {'Marker',MS{k},'MarkerSize',8};
    out2 = linestyle(k,LS,color);
    out2{2} = 'none';
    out = [out1,out2];
    
end

function [x,y,numpoints,num_data_sets,mX,mY,LegStr,color] = ...
        check_and_update_data(x,y,numpoints,Legends,color)

    if ~iscell(Legends)
        if strcmp(Legends,'none')||strcmp(Legends,'default')
            Legends = {Legends};
        else
            error(['Legends must either be set to ''default'',',...
                ' ''none'', or be of type cell'])
        end
    end
    
    if ~iscell(color)
        if strcmp(color,'default')
            color = {color};
        elseif size(color,2)==num
            error 'There must be a color for each plot inside a cell...'
        end
    end
    
    if numel(numpoints)>max(size(numpoints))
        error('The number of points must be a single row/column vector')
    end
    
    if size(numpoints,1)>1,numpoints = numpoints';end
    
    [m,n] = size(x);
    [M,N] = size(y);
    if min(m,n)>1
        error('The x axis must be consistant')
    end
    
    if max(numpoints)>length(x)
        error('Can''t have more numpoints then there are data points...')
    end
    
    if m>1
        num_data_sets = N;
    else
        num_data_sets = M;
        x = x';
        y = y';
    end
    
    if length(numpoints)>1&&length(numpoints)~=num_data_sets
        error('If numpoints>1, then numpoints==num_data_sets')
    end
    
    if num_data_sets>length(numpoints)
        numpoints = repmat(numpoints,num_data_sets,1);
    end
    
    mX = range(x)/2;
    mY = range(y(:))/2;
    
    LegStr = getLegends(Legends,num_data_sets);
    
    if length(color)<num_data_sets && length(color)>1
        error 'There must be a color for each plot...'
    elseif strcmp(color{1},'default')
        for m=num_data_sets:-1:1
            color{m} = [0 0 0];
        end
    end
    
end

function LegStr = getLegends(Legends,num_data)

    if strcmp(Legends{1},'default')
        txt = 'data ';
        LegStr = cell(num_data,1);
        for m = 1:num_data
            LegStr{m} = [txt,num2str(m)];
        end
    elseif strcmp(Legends{1},'none');
        LegStr = Legends;
    else
        Legends = Legends(:);
        if length(Legends)~=num_data
            error('There must be a legend title per data set')
        else
            LegStr = Legends;
        end
    end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%