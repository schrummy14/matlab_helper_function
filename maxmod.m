function maxV = maxmod(x,y)

    isSameSign = x.*y>=0;
    minValue = max(abs(x),abs(y));
    maxV = isSameSign.*sign(x).*minValue;

end