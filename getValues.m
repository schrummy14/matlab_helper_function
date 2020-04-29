function [time,solar,uwind,vwind,wwind] = getValues(data)

    values = data.Sensor.Data;
    
    time  = values(:,1);
    time  = timefix(time);
    solar = values(:,2);
    uwind = values(:,3);
    vwind = values(:,4);
    wwind = values(:,5);

end

function tout = timefix(tin)

    done = false;

    while ~done
        A = [];
        n = 1;
        for i = 2 : (length(tin)-1)
            if tin(i+1)-tin(i)<0
                A(n) = i;
                n=n+1;
            end
        end

        if ~isempty(A)
            for i = 1:length(A)
                tin(A(i)) = tin(A(i))-1*10^(floor(log10(tin(A(i)))));
            end
        else
            done = true;
        end
    end

    tout = tin;
    
end