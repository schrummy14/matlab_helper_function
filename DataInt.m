function int = DataInt(x,y,bound,tol)
%% int = DataInt(x,y,bound,tol)
% DataInt is a numerical integrater used for integrating data. DataInt does
% this by creating a cubic spline over the data and integrating the
% corresponding polynomial.
%
% x == Is the location of the Data Points
% y == Data Values
% bound == Where in the data you would like to integrate
% tol == The tolerence used in calculating the integral

if numel(x)~=numel(y)
    error('x and y must contain the same number of elements')
end
g = size(x);
if min(g) > 1 
    error('DataInt only integrates row vectors or column vectos at the moment')
end

if size(x,1)>1
    x = x';
end
if size(y,1)>1
    y = y';
end

switch nargin()
    case 2
        a = x(1);
        b = x(end);
        tol = 1e-8;
    case 3
        if numel(bound) == 2
            a = bound(1);
            b = bound(2);
        elseif numel(bound) == 1
            a = x(1);
            b = bound;
        else error('Bound must contain the start and end of your integration')
        end
        tol = 1e-8;
    case 4
        if isempty(tol),tol = 1e-8;end
        if isempty(bound)
            a=x(1); b=x(end);
        else
            a=bound(1); b=bound(2);
        end
end

if a<x(1) || b>x(end)
    error('You are integrating outside your data set.');
end

dy1 = gradient(y(1:2),diff(x(1:2)));
dy2 = gradient(y(end-1:end),diff(x(end-1:end)));
pp = spline(x,[dy1(1),y, dy2(end)]);

int = integral(@(x)ppval(pp,x),a,b,'AbsTol',tol);

