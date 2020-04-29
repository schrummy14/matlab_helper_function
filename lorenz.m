function dxdt = lorenz(~,a)    

    %     a = 10;
    %     c = 8/3;
    %     b = 28;


    sigma=10;
    beta=8/3;
    rho=28;
    dxdt = [
       -sigma*a(1) + sigma*a(2); 
        rho*a(1) - a(2) - a(1)*a(3); 
       -beta*a(3) + a(1)*a(2)];
   
end