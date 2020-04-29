function [orientation, data] = get_fiber_orientation(num_par_per_fiber, filename)

    if nargin < 2
        [fn, fp] = uigetfile('*.*');
        filename = [fp,fn];
        if nargin < 1
            num_par_per_fiber = 4;
        end
    end
    
    data = get_positions(filename);
    
    Num_Par = size(data,1);
    Num_Fib = Num_Par/num_par_per_fiber;
    
    orientation = zeros(Num_Fib,2);
    
    for k = 1:Num_Fib
        id_1 = (k-1)*num_par_per_fiber+1;
        id_2 = (k-1)*num_par_per_fiber+num_par_per_fiber-1;
        val = data(id_2,:) - data(id_1,:);
        x = val(1);
        y = val(2);
        z = val(3);

        r = sqrt(x*x + y*y + z*z);
        t = atand(y/x);
        p = acosd(z/r);
        orientation(k,:) = [t, p];

    end

end


function data = get_positions(filename)
    
    fid = fopen(filename);
    
    fgetl(fid);fgetl(fid);fgetl(fid);
    line = fgetl(fid);
    num_atoms = eval(line);
    
    data = zeros(num_atoms, 3);
    
    fgetl(fid);fgetl(fid);fgetl(fid);fgetl(fid);
    line = fgetl(fid);
    line_split = strsplit(line,' ');
    id_x = find(strcmp(line_split,'x')) - 2;
    for k = 1:num_atoms
        line = fgetl(fid);
        line_split = eval(['[',line,']']);
        id = line_split(1);
        x = line_split(id_x);
        y = line_split(id_x+1);
        z = line_split(id_x+2);
        data(id,:) = [x y z];
    end
    
    fclose(fid);
end