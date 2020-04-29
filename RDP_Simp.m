function out = RDP_Simp(points,myEps)

    outStart = points(1,:);
    outEnd = points(end,:);
    
    maxD = 0;
    for k = 2:(size(points,1)-1)
        a = [outStart-outEnd,0];
        b = [points(k,:)-outEnd,0];
        d = norm(cross(a,b)) / norm(a);
        if d > maxD 
            idVal = k;
            maxD = d; 
        end
    end
    
    if maxD > myEps
        out = [
            RDP_Simp(points(1:idVal,:),myEps);
            RDP_Simp(points(idVal+0:end,:),myEps)];
        out(find(points(idVal,:) == out,1),:) = [];
    else
        out = [outStart;outEnd];
    end            

end