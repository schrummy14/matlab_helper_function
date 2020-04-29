function [B]=LSE(Y,X)
% -----------------------------------------------------------
% General Usage: to estimate unknown parameters of a linear equation
%                [B]=LSE(Y,X)
%
% Inputs:     Y: observation matrix
%             X: inputs
%
% Outputs:    B: unknown parameters
%
% Author: Farhad Sedaghati, PhD student at The University of Memphis
% Contact: <farhad_seda@yahoo.com>
% Written: 09/30/2015
%
% References:
%     1) Ljung, L., System Identification: Theory for the User, Second ed.,
%        Prentice Hall, 1999
% -----------------------------------------------------------





    [n,m]=size(X);

    % define a diagonal matrix with large number elements
    s=10^5*eye(m,m);

    B=zeros(m,1);

    for k=1:n
        s=s-(s*(X(k,:))'*X(k,:)*s)/(1+X(k,:)*s*(X(k,:))');
        B=B+s*(X(k,:))'*(Y(k,1)-X(k,:)*B);
    end


end