function dt = get_dt_hertz(radi,Y,density,maxV)

    
    Yeff = 0.5 .* Y;
    meff = 4/3 .* pi .* radi.^3 .* density;
    reff = 0.5 .* radi;
    
    dt = 2.87 .* (meff.^2 ./ (reff.*Yeff.*Yeff.*maxV)).^0.2; 
    
end