function x = diffx(v,n,dim)

if nargin < 3
    dim = 1;
end
if dim == 2
    v = v.';
    x = diffx(v,n,1);
    x = x.';
    return
end

if size(v,1) < 2
    v = v.';
    flipFlag = true;
else
    flipFlag = false;
end

x = v(n+1:end,:) - v(1:end-n,:);

if flipFlag
    x = x.';
end

end
