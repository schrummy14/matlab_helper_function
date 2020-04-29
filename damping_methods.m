function varargout = damping_methods
    
    n_bonds = 55;
    diam = 2.83e-3;
    Y = 1.4e9;
    rho = 420;
    beta = .0005;
    radi = 0.5*diam;
    K = Y*pi*radi*radi/diam;
    m = rho * 4/3 * pi * radi*radi*radi;
    F_Pull = 1;
    
    dt = 1.0e-7;
    count_out = 0.0005/dt;
    t_end = 0.025;
    
    n_parts = n_bonds+1;
    t = 0;
    x = diam*(0:n_parts-1)';
    v = zeros(size(x));
    bondhist = zeros(n_bonds,1);
    F = zeros(size(x));
    F(end) = F_Pull;
    count = 0;
    count_val = count_out;
    clf
    F_out = [];
    v_out = [];
    x_out = [];
    t_out = [];
    
    F_out(:,end+1) = F;
    v_out(:,end+1) = v;
    x_out(:,end+1) = x;
    t_out(1,end+1) = t;
    
    sqrt_K_m = sqrt(K*m);
    two_beta_sqrt_K_m = 2*beta*sqrt_K_m;
    Kdt = K*dt;
    
    while t < t_end
        
        v = v + 0.5*dt*F/m;        
        x = x + dt*v;
        
        F = 0*F;
        for k = 1:n_bonds
        
            % force caused by bond
            dv = v(k)-v(k+1);
            df = -Kdt*dv;
            bondhist(k) = bondhist(k) + df;
            
            Eke = m*dv*dv*sqrt(K);
            Epe = bondhist(k)*bondhist(k);
            E = 0.5*(Eke + Epe);
            f_damp = -beta*sqrt(E);
%             f_damp = -two_beta_sqrt_K_m*dv;

            F(k) = F(k) + bondhist(k) + f_damp*sign(v(k));
            F(k+1) = F(k+1) + -bondhist(k) - f_damp*sign(-v(k+1));
            
            % force caused by contact
            dx = abs(x(k)-x(k+1)) - diam;
            if dx < 0
                F_contact = -K*dx - 0.001*dv;
                F(k) = F(k) + F_contact;
                F(k+1) = F(k+1) - F_contact;
            end
            
        
        end
        
        F(end) = F(end) + F_Pull;
        F(1) = 0;
        
        if any(isnan(F))
            disp(['nan values present at step ', num2str(count)])
            break
        end
        
        v = v + 0.5*dt*F/m;
        t = t + dt;
        
        count = count + 1;
        if count >= count_val || t >= t_end
            count_val = count_val + count_out;
            F_out(:,end+1) = F;
            v_out(:,end+1) = v;
            x_out(:,end+1) = x;
            t_out(1,end+1) = t;
            subplot(3,1,1)
            plot(t_out,F_out(end,:),'o')
            title(['time == ',num2str(t)])
            ylabel Force
            subplot(3,1,2)
            plot(t_out,v_out(end,:),'o')
            ylabel Vel
            subplot(3,1,3)
            plot(t_out,x_out(end,:),'o')
            ylabel Pos
            drawnow
        end
        
    end
    
    switch nargout
        case 1
            varargout{1} = [t_out.',x_out.',v_out.',F_out.'];
        case 2
            varargout{1} = t_out.';
            varargout{2} = x_out.';
        case 3
            varargout{1} = t_out.';
            varargout{2} = x_out.';
            varargout{3} = v_out.';
        case 4
            varargout{1} = t_out.';
            varargout{2} = x_out.';
            varargout{3} = v_out.';
            varargout{4} = F_out.';
    end
    
end
    
    