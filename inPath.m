function out = inPath(wantedPath)

    
    out = false;
    allPaths = strsplit(path,';');
    
    for k = 1:length(allPaths)
        if strcmp(wantedPath,allPaths{k})
            out = true;
            break;
        end
    end
    
end