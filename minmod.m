function minV = minmod(x,y)

    isSameSign = x.*y>=0;
    minValue = min(abs(x),abs(y));
    minV = isSameSign.*sign(x).*minValue;

end
        