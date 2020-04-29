function [D,V] = Jacobi_eig(A,epsilon,show)

    if nargin==2, show = 0; end
    D = A;
    [m,n] = size(A);
    if m~=n,error('A must be a square matrix');end
    V = eye(n);
    cnts = 0;
    cntr = 0;
    done = 0;
    working = 1;
    stat = working;
    while (stat==working),
        cnts = cnts+1;
        t = sum(diag(D));
        stat = done;
        for p = 1:(n-1),
            for q = (p+1):n,
                if ((abs(D(p,q))/t)>epsilon),
                    t = D(p,q)/(D(q,q) - D(p,p));
                    c = 1/sqrt(t*t+1);
                    s = c*t;
                    R = [c s; -s c];
                    D([p q],:) = R'*D([p q],:);
                    D(:,[p q]) = D(:,[p q])*R;
                    V(:,[p q]) = V(:,[p q])*R;
                    cntr = cntr+1;
                    if show==1,
                        home; if cntr==1, clc; end; 
                        disp(['Jacobi iteration No. ',int2str(cntr)]),disp(''),...
                        disp(['Zeroed out the element  D(',num2str(p),...
                              ',',num2str(q),') = ']),disp(D(p,q)),...
                        disp('New transformed matrix  D = '),disp(D)
                    end
                    stat = working;
                end
            end
        end
    end
    ev = diag(D);
    [ev,idx] = sort(-ev);
    ev = -ev;
    if nargout > 1
        D = diag(ev);
        V = V(:,idx);
    else
        D = ev;
    end
%     D = diag(diag(D));
    
end