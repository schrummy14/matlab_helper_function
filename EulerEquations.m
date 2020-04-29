function EulerEquations
    
    global gamma
    gamma = 1.4;
    m = 101;
    n = 101;
    U = zeros(m,n,4);
    rho1 = 1.4;
    rho2 = 1*1.4;
    P1   = 10;
    P2   = 11;
    u1   = 0;
    u2   = 0;
    v1   = 0;
    v2   = 0;
    U1 = [rho1,rho1*u1,rho1*v1,P1/(gamma-1)+rho1*(u1^2+v1^2)/2];
    U2 = [rho2,rho2*u2,rho2*v2,P2/(gamma-1)+rho2*(u2^2+v2^2)/2];
    
    U(:,:,1) = U1(1);
    U(:,:,2) = U1(2);
    U(:,:,3) = U1(3);
    U(:,:,4) = U1(4);
    
    U(:,1:10,1) = U2(1);
    U(:,1:10,2) = U2(2);
    U(:,1:10,3) = U2(3);
    U(:,1:10,4) = U2(4);
    
    dx = 0.1;
    dy = 0.1;
    figure(1) 
    surf(U(:,:,1),'LineStyle','none')
    view(2)
    drawnow
    t_out = 0.001;
    t_count = 0;
    t = 0;
    while t < 0.5
        dt = getdt(dx,dy,U);
        U = updateU(dx,dy,dt,U);
        t = t + dt;
        if imag(t) ~= 0
            break
        end
        if t > t_count
            t_count = t_count + t_out;
            figure(1) 
            surf(real(U(:,:,1)),'LineStyle','none')
            title(['Time == ', num2str(t)])
            view([0 0])
            drawnow
        end
    end
    
end

function dt = getdt(dx,dy,U)

    global gamma
    [r,u,v,p]=getQ(U);
    a = sqrt(gamma.*p./r);
    max_u_eig = max(abs(u(:))+a(:));
    max_v_eig = max(abs(v(:))+a(:));
    dt = 0.1/(max_u_eig/dx+max_v_eig/dy);
    
end

function U = updateU(dx,dy,dt,U)

    Ix = 2:size(U(:,:,1),2)-1;
    Iy = 2:size(U(:,:,1),1)-1;

    Uip1 = U(Iy,Ix+1,:);
    Uim1 = U(Iy,Ix-1,:);
    Uij  = U(Iy,Ix,:);
    
    Eip1 = getE(Uip1);
    Eim1 = getE(Uim1);
    
    U(Iy,Ix,:) = Uij - dt/(2*dx)*(Eip1-Eim1);
    
    Ujp1 = U(Iy+1,Ix,:);
    Ujm1 = U(Iy-1,Ix,:);
    Uij  = U(Iy,Ix,:);
    
    Fjp1 = getF(Ujp1);
    Fjm1 = getF(Ujm1);    
    
    U(Iy,Ix,:) = Uij - dt/(2*dy)*(Fjp1-Fjm1);
    
    U = updateBoundary(U);

end

function U = updateBoundary(U)

    U(1,:,:) = U(2,:,:);
    U(end,:,:) = U(end-1,:,:);
    U(:,1,:) = U(:,2,:);
    U(:,end,:) = U(:,end-1,:);
    
end

function E = getE(U)

    [rho,u,v,P] = getQ(U);
    E(:,:,1) = rho.*u;
    E(:,:,2) = rho.*u.^2+P;
    E(:,:,3) = rho.*u.*v;
    E(:,:,4) = u.*(U(:,:,4)+P);
    
end

function F = getF(U)

    [rho,u,v,P] = getQ(U);
    F(:,:,1) = rho.*u;
    F(:,:,2) = rho.*u.*v;
    F(:,:,3) = rho.*v.^2+P;
    F(:,:,4) = v.*(U(:,:,4)+P);

end

function [rho,u,v,P] = getQ(U)
    
    global gamma
    rho = U(:,:,1);
    u   = U(:,:,2)./rho;
    v   = U(:,:,3)./rho;
    P   = (gamma-1).*(U(:,:,4)-rho.*(u.^2+v.^2)./2);
    
end