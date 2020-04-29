function [xmin, fxmin] = multiParForMin(myFun, numWorkers, numIts, LV, HV)

    curPool = gcp('nocreate');
    if isempty(curPool)
        curPool = parpool(numWorkers);
    else
        numWorkers = curPool.NumWorkers;
    end

    numVar = length(LV);
    minFx = inf(numWorkers,1);
    x = zeros(numWorkers,numVar);
    
    LV = (LV(:))';
    HV = (HV(:))';
    
    fxTry = myFun(LV);
    if length(fxTry) > 1
        opts = optimoptions('lsqnonlin');
        optFun = @(x0,opts)lsqnonlin(myFun,x0,LV,HV,opts);
    else
        opts = optimoptions('fmincon');
        optFun = @(x0,opts)fmincon(myFun,x0,[],[],[],[],LV,HV,[],opts);
    end
    opts.Display = 'None';
    opts.MaxFunctionEvaluations = 1000*length(LV);
    opts.MaxIterations = opts.MaxIterations*10;
    opts.FunctionTolerance = opts.FunctionTolerance/100;
    startTime = tic;
    ticBytes(curPool)
    parfor coreID = 1:numWorkers
        curTime = 0;
        oldRemainTime = 0;
        for k = 1:numIts
            if coreID == 1
                loopStart = tic;
            end
            x0 = (HV - LV).*rand(1,numVar) + LV;
            try
                [CurX, fx] = optFun(x0,opts);
            catch
                disp(['Failed at x0 = ', num2str(x0)]);
                fx = inf;
            end
            if fx < minFx(coreID)
                minFx(coreID) = fx;
                x(coreID,:) = CurX;
                fprintf('Core %i found fx = %0.3e with x = [', coreID, fx);
                fprintf(' %0.3f',CurX);
                fprintf(' ]\n');
                
            end
            if coreID == 1
                curTime = (k-1)*curTime + toc(loopStart);
                curTime = curTime/k;
                remainTime = (numIts-k+1)*curTime;
                if abs(remainTime - oldRemainTime) > 5
                    fprintf('Estimated time remaining = %0.3f\n',remainTime);
                    oldRemainTime = remainTime;
                end
            end
        end
    end
    tocBytes(curPool)
    disp(['Total Elapsed Time = ', num2str(toc(startTime))])
    [fxmin, minID] = min(minFx);
    xmin = x(minID,:);

end