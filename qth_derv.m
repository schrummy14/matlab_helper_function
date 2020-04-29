function qx = qth_derv(q,f,h)

%     qx = 0;
%     qx_old = 0;
%     done = false;
%     m = 0;
%     while ~done
%         qx = qx + 1/h^q*(-1)^m.*qchoosem(q,m).*fx(x+(q-m).*h);
%         err = qx - qx_old;
%         if norm(err,inf)<1e-4
%             done = true;
%         end
%         qx_old = qx;
%         m = m + 1;
%     end
% 
% end
% function out = qchoosem(q,m)
% 
%     out = gamma(q+1)./(gamma(m+1).*gamma(q-m+1));
% 
% end

% fgl_deriv
%
%   Computes the fractional derivative of order alpha (a) for the function
%   y sampled on a regular grid with spacing h, using the Grunwald-Letnikov
%   formulation.
%
%   Inputs:
%   a : fractional order
%   y : sampled function
%   h : period of the sampling lattice
%
%   Output:
%   Y : fractional derivative estimate given y
%
%   Note that this implementation is similar to that of Bayat 2007
%   (FileExchange: 13858-fractional-differentiator), and takes the exact
%   same inputs, but uses a vectorized formulation to accelerates the
%   computation in Matlab.
%
%   Copyright 2014 Jonathan Hadida
%   Contact: Jonathan dot hadida [a] dtc.ox.ac.uk

    n  = numel(f);
    J  = 0:(n-1);
    G1 = gamma( J+1 );
    G2 = gamma( q+1-J );
    s  = (-1) .^ J;

    M  = tril( ones(n) );
    R  = toeplitz( f(:)' );
    T  = meshgrid( (gamma(q+1)/(h^q)) * s ./ (G1.*G2) );
    qx  = reshape(sum( R .* M .* T, 2 ), size(f));
    
end