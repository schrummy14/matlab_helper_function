% surface_fit to apply a fit to a n-dimensional dataset such that the fitted
% line best fits the data while being strictly equal to or greater than the
% data set. 
% surface_fit solves a linear program to findthe best solution with the
% conditions
% minimize E such that Ax > b and Ax - E = b
% where 
%   E is the errors from the fit to the dataset
%   A is a matrix consisting of the linear terms plus a constant
%   x are the coefficients that we would like to find
%   b is the dataset that we would like to fit
%
% Usage
% bnew = surface_fit(X,y)
%   Applies the algorithm to the dataset y using the linear coefficients in
%   X
%
% [bnew, res] = surface_fit(X,y)
%   Same as above but also returns the raw residules
% [bnew, res, fval] = surface_fit(X,y)
%   Same as above but also returns the function evaluation at bnew
%
% ~ = surface_fit(X,y,over_under)
%   over_under is a flag that tell the program if you would like the fitted
%   model to be strictly greater than the dataset ('over') (default), or
%   strictly less than the dataset ('under')

function [bnew,res,fval] = surface_fit(X, y, over_under, do_const)
    
    % Check to see if we should just run the demo
    if nargin < 1
        over_under = 'over';
        run_demo(over_under);
        return
    end
    
    % Check if the under_over flag was passed, if not use default values
    if nargin < 4
        do_const = true;
        if nargin < 3
            over_under = 'over';
        end
    end
    
    % If over_under is set to 'under', make changes and rerun program
    if strcmpi(over_under,'under')
        [bnew,res] = surface_fit(X, -y,'over',do_const);
        
        % To get the correct outcome, we must negate the results
        res = -res;
        bnew = -bnew;
        return
    end

    % Get the size of X
    [m,n] = size(X);

    % Build the linear matrix (To return the least squares result == A1\b)
    if do_const
        A1 = [ones(m,1),X];
    else
        A1 = X;
        n = n-1;
    end
    
    % Build equality matrix and conditional matracies
    Aeq = [A1, -speye(m)];
    A   = [A1, sparse(m,m)];
    b   = y(:);
    f   = [sparse(1,n+1), ones(1,m)];
    
    % Solve linear program
    [bb, fval] = linprog(f,-A,-b,Aeq,b);
    
    if isempty(bb)
        error(' ')
    end
        
    % Get coefficients and residules
    bnew = bb(1:n+1);
    res  = bb(n+2:end);

end

function run_demo(over_under)

    x = linspace(-3*pi,3*pi,20).';
    y = cos(x);

    subplot(2,1,1)
    A = [ones(length(x),1),x];
    y1 = A*(A\y);
    A = [A,x.^2];
    y2 = A*(A\y);
    A = [A,x.^4];
    y3 = A*(A\y);
    plot( ...
        x,y,'ko', ...
        x,y1,'b-', ...
        x,y2,'r--', ...
        x,y3,'g.-')
    legend('Data','Constant','Quadratic Model', '4^{th} Order Model')
    axis([-10 10 -1 2])
    
    
    subplot(2,1,2)
    X = x;
    [~, res] = surface_fit(X,y,over_under);
    y1 = y + res;
    X = [X,x.^2];
    [~, res] = surface_fit(X,y,over_under);
    y2 = y + res;
    X = [X,x.^4];
    [~, res] = surface_fit(X,y,over_under);
    y3 = y + res;
    plot( ...
        x,y,'ko', ...
        x,y1,'b-', ...
        x,y2,'r--', ...
        x,y3,'g.-')
    legend('Data','Constant','Quadratic Model', '4^{th} Order Model')
    axis([-10 10 -1 2])
end