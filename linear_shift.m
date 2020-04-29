function [out] = linear_shift(z,a,b,c,d)

    out = (z-a).*(d-c)./(b-a) + c;

end