function x = Bi_CG_Stab_L(A,b,ell,TOL,MaxIts)

    if nargin == 0
        n = 150;
        for mm = 1:n
            for nn = 1:n
                if nn+1==mm || nn-1==mm
                    A(mm,nn) = -24;
                elseif nn+2==mm || nn-2==mm
                    A(mm,nn) = -120;
                elseif nn+5==mm || nn-5==mm
                    A(mm,nn) = -12;
                elseif nn==mm 
                    A(mm,mm) = 500;
                end
            end
        end
        
        b = sum(A,2);
        ell = 3;
        TOL = 1e-10;
        MaxIts = 100;
    end
    
    n2b = norm(b);
    tolb = TOL*n2b;
    n = length(b);
    x = zeros(n,1);
    x0 = x;
    flag = 1;
    imin = 0;
    
    r = b;
    normr = norm(r);
    normrmin = normr;
    
    ukml = 0;
    rk = r;
    xk = zeros(n,1);
    r0 = rk;
    rho0 = 1;
    omega = 1;
    ut = zeros(n,ell+1);
    rt = ut;
    alpha = 0;
    tau = zeros(ell+1,ell+1);
    sigma = zeros(ell+1,1);
    gamma = sigma;
    gammap = sigma;
    gammapp = sigma;
    
    for kk = 1:MaxIts
        ut(:,1) = ukml;
        rt(:,1) = rk;
        xt = xk;
        rho0 = -omega*rho0;
        for jj = 1:ell
            rho1 = r0'*rt(:,jj);
            beta = alpha*rho1/rho0;
            rho0 = rho1;
            ut(:,1:jj) = rt(:,1:jj)-beta*ut(:,1:jj);
            ut(:,jj+1) = A*ut(:,jj);
            gamma_s = r0'*ut(:,jj+1);
            alpha = rho0/gamma_s;
            rt(:,1:jj) = rt(:,1:jj) - alpha*ut(:,2:jj+1);
            rt(:,jj+1) = A*rt(:,jj);
            normr = norm(rt(:,1));
            xt = xt + alpha*ut(:,1);
            
            if normr<=tolb
                rt(:,1) = b-A*(x0+xt);
                normr_act = norm(rt(:,1));
                if normr_act<=tolb
                    flag = 0;
                    iter = ((kk-1)*ell*2+ell+1)/(2*ell);
                    break
                end
            end
            
            if normr<normrmin
                normrmin = normr;
                xmin = xt;
                imin = ((kk-1)*ell*2+jj)/(2*ell);
            end
        end
        if flag == 0
            break;
        end
        
        for jj = 2:ell+1
            for ii = 2:jj-1
                tau(ii,jj) = (rt(:,jj)'*rt(:,ii))/sigma(ii);
                rt(:,jj) = rt(:,jj) - tau(ii,jj)*rt(:,ii);
            end
            sigma(jj) = rt(:,jj)'*rt(:,jj);
            gammap(jj) = (rt(:,1)'*rt(:,jj))/sigma(jj);
        end
        gamma(ell+1) = gammap(ell+1);
        omega = gamma(ell+1);
        
        for jj = ell:-1:2
            gamma(jj) = gammap(jj) - tau(jj,jj+1:ell+1)*gamma(jj+1:ell+1);
        end
        
        for jj=2:ell
            gammapp(jj) = gamma(jj+1) + tau(jj,jj+1:ell)*gamma(jj+2:ell+1);
        end
        
        xt = xt + gamma(2)*rt(:,1);
        rt(:,1) = rt(:,1) - gammap(ell+1)*rt(:,ell+1);
        ut(:,1) = ut(:,1) - gamma(ell+1)*ut(:,ell+1);
        normr = norm(rt(:,1));
        
        if normr<=tolb
            rt(:,1) = b-A*(x0+xt);
            normr_act = norm(rt(:,1));
            if normr_act<=tolb
                flag = 0;
                iter = ((kk-1)*ell*2+ell+1)/(2*ell);
                break
            end
        end
        
        if normr<normrmin
            normrmin = normr;
            xmin = xt;
            imin = ((kk-1)*ell*2+ell+1)/(2*ell);
        end
        
        for jj = 2:ell
            ut(:,1) = ut(:,1) - gamma(jj)*ut(:,jj);
            xt = xt+gammapp(jj)*rt(:,jj);
            rt(:,1) = rt(:,1) - gammap(jj)*rt(:,jj);
            normr = norm(rt(:,1));
            
            if normr<=tolb
                rt(:,1) = b-A*(x0+xt);
                normr_act = norm(rt(:,1));
                if normr_act<=tolb
                    flag = 0;
                    iter = ((kk-1)*ell*2+ell+1)/(2*ell);
                    break
                end
            end

            if normr<normrmin
                normrmin = normr;
                xmin = xt;
                imin = ((kk-1)*ell*2+ell+1)/(2*ell);
            end
        end
        
        if flag == 0
            break;
        end

        ukml = ut(:,1);
        rk = rt(:,1);
        xk = xt;         
        
    end
    
    
    if flag == 0
        x = x0 + xt;
        relres = normr_act/n2b;
    else
        r1 = b - A*(x0+xmin);
        r2 = b - A*(x0+xt);
        if norm(r1)<norm(r2)
            x = x0 + xmin;
            iter = imin;
            relres = norm(r1)/n2b;
        else
            x = x0 + xt;
            iter = kk;
            relres = norm(r2)/n2b;
        end
    end
    disp(['Number of iterations == ',num2str(iter)])
    disp(['Relative Residule == ',num2str(relres)])
end