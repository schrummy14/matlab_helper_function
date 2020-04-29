function [data, data_name, sim_time] = read_liggghts_log(num_outs,file_name) %,stop_early)
    
    if nargin < 1, num_outs = 4; end
    
    CHUNK = 500;
    
    data = nan(CHUNK,num_outs);
    if nargin < 2
        [fn,pn] = uigetfile('*.liggghts');
        file_name = [pn,fn];
    end
    fid = fopen(file_name);
    fline = fgetl(fid);
    line_break = strsplit(fline);
    sim_time = 0;
    data_count = 1;
    
    while ~isnumeric(fline)
        
        done = false;
        while ~done
            if strcmp(['x',line_break{1}],'x') && numel(line_break)>1
                line_break = line_break(2:end);
            end
            if strcmp(['x',line_break{1}],'xStep')
                if size(line_break,2)==num_outs+1
                    done = true;
                else
                    fline = fgetl(fid);
                    line_break = strsplit(fline);
                end
            else
                fline = fgetl(fid);
                if fline == -1
                    fclose all;
                    data(data_count + 1:end,:) = []; % New
                    return;
                else
                    line_break = strsplit(fline);
                end
            end
        end
        
        data_name = strrep(line_break(1:end-1),'_',' ');
        fline = fgetl(fid);
%         line_break = str2num(fline);
        line_break = eval(['[',fline,']']);
        
        while ~isempty(line_break)
            if length(line_break) == num_outs
                data(data_count,:) = line_break;
                data_count = data_count + 1;
                if mod(data_count-1,CHUNK) == 0
                    data(end+1:end+CHUNK,:) = nan(CHUNK,num_outs);
                end
%                 data(end+1,:) = line_break;
                fline = fgetl(fid);
                if fline == -1
                    fclose all;
                    data(data_count + 1:end,:) = [];
                    return;
                else
                    try
                        line_break = eval(['[',fline,']']);
                    catch
                        line_break = []; % str2num(fline);
                    end
                end
            else
                fclose all;
                data(data_count + 1:end,:) = [];
                return;
            end
        end
        
        line_break = strsplit(fline);
        sim_time = sim_time + str2double(line_break{4});
        if isnan(sim_time)
            sim_time = old_sim_time;
        end
        old_sim_time = sim_time;
        
        fline = fgetl(fid);
        line_break = strsplit(fline);
        
    end
    
    fclose all;
    
end