function dt = get_dt_bond(radi,Y,density,lb)

    if nargin < 4
        lb = 2*radi;
    end
    m1 = 4/3 * pi * radi.^3;
    m2 = m1;
    meff = m1.*m2./(m1+m2); % m*m/(2*m); m/2;
    k = pi.*Y*radi.^2 ./ lb;
    m = meff .* density;
    dt = sqrt(m/k);
    
end