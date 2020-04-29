function [ynew,rsq,beta] = LinReg(datax,datay,order,xx)
    flip = 0;
    % datax and datay must be column vectors
    if size(datax,2) > 1, datax = datax';flip=1; end
    if size(datay,2) > 1, datay = datay';flip=1; end
    if exist('xx','var')
        if size(xx,2)>1, xx = xx.';flip = 1;end
    end
    sM = size(datax);
    % Construct the A matrix
    A = ones(length(datax),order+1);
    for m = 1:order
        A(:,m+1) = datax.^m;
    end
    
    % Solve A'y = A'Ab to give beta values
%     beta = (A'*A)\(A'*datay);
    beta = A\datay;
    
    % Give ynew approx values
    % if xx is passed, use its values
    if exist('xx','var')
        xxM = size(xx);
        ycheck = beta(1)*ones(sM(1),sM(2));
        ynew = beta(1)*ones(xxM(1),xxM(2));
        for n = 2:length(beta)
            ycheck = ycheck + beta(n)*datax.^(n-1);
            ynew = ynew + beta(n)*xx.^(n-1);
        end
    else
        ynew = beta(1)*ones(sM(1),sM(2));
        for n = 2:length(beta)
            ynew = ynew + beta(n)*datax.^(n-1);
        end
        ycheck = ynew;
    end
    
    % Calculate r squared 
    rsq = corr(datay,ycheck)^2;
    
    if flip
        ynew = ynew.';
    end
end