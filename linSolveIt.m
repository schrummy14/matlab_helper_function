function [x,tau,its] = linSolveIt(A,b,TOL)

    if ~exist('TOL','var')
        TOL = 1e-12;
    end
    
    its = 0;
    xold = zeros(size(A,1),1);
    xnew = xold;
    
    done = 0;
    
    while ~done
        n = length(A);
        x = xold;
        relax = 1;
        for i = 1:n
            x(i) = (1/A(i, i))*(b(i) - A(i, 1:n)*x + A(i, i)*x(i));
        end
        if relax == 1
            xnew = x;
        else
            xnew = (1-relax)*xold + relax*x;
        end
        tau = norm(A*xnew - b);
        xdif = norm(xnew-xold);
        
        if tau < TOL 
            done = 1;
        elseif xdif == 0
			done = 1;
			msg1 = 'There is no difference between iterations, so the solver has stopped...';
			msg2 = 'Final error value of ';
			msg3 = num2str(tau);
			msg4 = ' was reached.';
			warning(msg1,msg2,msg3,msg4);
            xold = xnew;
        end
        its = its + 1;
        xold = xnew;
    end
    
    x = xnew;
end
        