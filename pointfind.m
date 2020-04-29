function [per,clevels,points] = pointfind(x,y,levels,lognum,varargin)

%% Function [points,per] = pointfind(x,y,levels,lognum)
%
%  This function takes two data sets and finds the frequency in which
%  the data intersects and plots a contour plot of the given data. It will
%  then provide how many points fall within an evenly spaced level.
%
%  pointfind(x,y) returns a plot of original data over 5 levels.
%
%  pointfind(x,y,n) returns a plot of original data over n levels. n must
%  be an integer greater than 1.
%
%  pointfind(x,y,0) returns a plot of the original data over 5 levels.
%  Same as pointfind(x,y)
%
%  pointfind(x,y,1) returns a plot of the narural log of the frequency data
%  to see smaller changes in the frequency ploted over 5 levels.
%
%  pointfind(x,y,n,1) returns a plot of the natural log of the frequency
%  data to see smaller changes in the frequency ploted over n levels.
%
%  Example:
%  len = 1000;
%  x = floor(25*rand(len,1));
%  y = floor(25*rand(len,1));
%  [per,clevles,points] = pointfind(x,y,10,0);

%% Checking Input
%  Check to see which conditions need to be satisfied
switch nargin
    case 2
        levels = 5;
        lognum = 0;
    case 3
        if levels == 1
            lognum=1;
            levels = 5;
        elseif levels == 0
            lognum=0;
            levels = 5;
        else
            lognum=0;  
        end
end

if length(x)~=length(y)
    error('X and Y must be the same length');
elseif levels < 1
    error('Levels must be greater than one.');
elseif lognum ~= 0 && lognum ~= 1
    error('Log flag must be either 0 == off, or 1 == on')
end
       
%% Find Frequency data
%figure()
[edgesX,edgesY,N]=ndhist(x,y,varargin);
%  If we are using a ln, we must "fix" the frequency data.

if lognum == 1
    N = log(N);
    N=N+1; % Since all bins that had a value of 1 were mapped to zero, We 
           % will add one to the bins to account for the original loss.
           % This will not effect the NaN values.
    N(N == -inf)=NaN; % Set -inf values to Not a Number values. 
    g=N;
    
    if levels > 1
        temp = linspace(min(min(N)),max(max(N)),levels+1);
        clevels=temp(2:end-1);
    elseif levels == 1
        clevels = 1;
    end

    contourf(edgesX,edgesY,g,[1,clevels],'LineStyle','none')
%    g(g<log(5))=NaN;
else
    g=N;
    
    if levels > 1
        temp = linspace(min(min(N)),max(max(N)),levels+1);
        clevels=temp(2:end-1);
    elseif levels == 1
        clevels = 1;
    end
    
    contourf(edgesX,edgesY,g,[1,clevels],'LineStyle','none')
%    g(g<5)=NaN;
end

points= zeros(numel(clevels),1);

for i = 1 : numel(clevels)+1
    if i == 1 && numel(clevels)>1
        points(i) = sum(N(N<clevels(1)));
    elseif i == numel(clevels)+1
        points(i) = sum(N(N>=clevels(i-1)));
    else
        if numel(clevels)>1
            points(i) = sum(N(N>=clevels(i-1) & N<clevels(i)));
        end
    end
end

tot = sum(points);

per = points*100/tot;        