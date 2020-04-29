function [nt,nx,t_set,x_set] = normalize_freq(t,x)

    ff = isnan(x);
    id = find(ff==false,1);
    x = x - x(id);
    [a,id_min] = min(x);
    [b,id_max] = max(x);
    if abs(a)>abs(b)
        id = id_min;
    else
        id = id_max;
    end

    x = x(id:end);
    t = t(id:end);
    t = t - t(1);
   
    [aa,bb] = get_freq(t,x);
    [~,id] = max(bb(2:end));
    id = id + 1;
    
    x_set = x(1);
    t_set = 1/aa(id);
    nx = x/x(1);
    nt = t*aa(id);

end