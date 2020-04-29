function ydot = threeBody(~,y,m0,m1,m2,dim)

%*****************************************************************************
%
%% SIMPLE_F returns the right hand side of the three body ODE system.
%
%  Modified:
%
%    03 April 2011
%
%  Author:
%
%    Dominik Gruntz, Joerg Waldvogel
%
%  Reference:
%
%    Dominik Gruntz, Joerg Waldvogel,
%    "Orbits in the Planar Three-Body Problem",
%    Walter Gander, Jiri Hrebicek,
%    Solving Problems in Scientific Computing using Maple and Matlab,
%    Springer, 1997,
%    ISBN: 3-540-61793-0,
%    LC: Q183.9.G36.
%
%  Parameters:
%
%    Input, real T, the current time.
%
%    Input, real Y(12), the current solution.
%
%    Output, real YDOT(12), the derivatives of the current solution.
%  
%  Example
%    [t,x] = rk78(@(t,x)threeBody(t,x,5,3,4,2),[0 50],[1 -1 0 0 1 3 0 0 -2 -1 0 0],'TOL',1e-12);

    x0 = y(0*dim+1:1*dim);
    x1 = y(1*dim+1:2*dim);
    x2 = y(2*dim+1:3*dim);
  
    d0 = ( x2 - x1 ) / norm ( x2 - x1 )^3;
    d1 = ( x0 - x2 ) / norm ( x0 - x2 )^3;
    d2 = ( x1 - x0 ) / norm ( x1 - x0 )^3;
    
    ydot(1:3*dim) = y(3*dim+1:6*dim);
    ydot(3*dim+1:4*dim) = m1 * d2 - m2 * d1;
    ydot(4*dim+1:5*dim) = m2 * d0 - m0 * d2;
    ydot(5*dim+1:6*dim) = m0 * d1 - m1 * d0;
  
    ydot = ydot(:);

    return
end