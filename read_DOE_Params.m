function [Params, Param_Names] = read_DOE_Params(filename)

    if nargin < 1
        filename = 'DOE_Params.txt';
    end
    
    fid = fopen(filename);
    
    Num_Params = 0;
    Params = [];
    Param_Names = [];
    
    while true
        line = fgetl(fid);
        if line == -1
            break
        end
        spt_line = strsplit(line);
        if strcmp(spt_line{1},'variable')
            Num_Params = Num_Params + 1;
            Param_Names{Num_Params,1} = spt_line{2};
            count = 0;
            while true
                count = count + 1;
                line = fgetl(fid);
                spt_line = strsplit(line);
                Params(count,Num_Params) = eval(spt_line{1});
                if length(spt_line) < 2
                    break
                end
            end
        end
    end
    
    fclose(fid);
    
end