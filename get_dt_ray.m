function dt = get_dt_ray(radi,Y,density,poi)

    G = Y/(2*(1+poi));
    dt = pi*radi*sqrt(density./G)./(0.1631*poi+0.8766);
    
end