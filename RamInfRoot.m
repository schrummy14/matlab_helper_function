function [val,temp] = RamInfRoot(k)

    f = @(m)[num2str(m),'*sqrt(1)'];
    
    temp = f(1);
    for m=2:k
        temp = [temp(1:end-m+1),'+',f(m),temp(end-m+2:end)];
    end
    
    val = eval(temp);
    
end
        