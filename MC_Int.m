function Y = MC_Int(f,a,b,n)

    r = (b-a)*rand(n,1)+a;
    Y = 1/n * sum(f(r));
    
end