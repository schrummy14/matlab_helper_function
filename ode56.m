function varargout = ode56(ode,tspan,y0,options,varargin)

    solver_name = 'ode56';
    
    if nargin < 4
        options = [];
        if nargin < 3
            y0 = [];
            if nargin < 2
                tspan = [];
                if nargin < 1
                    error('Not enough inputs...')
                end
            end
        end
    end
    
    % Stats 
    nsteps  = 0;
    nfailed = 0;
    nfevals = 0;
    
    % Output
    FcnHandlesUsed = isa(ode,'function_handle');
    output_sol = (FcnHandlesUsed && (nargout == 1));
    output_ty = (~output_sol && (nargout > 0));
    
    
    sol = []; f3d = [];
    if output_sol
        sol.solver = solver_name;
        sol.extdata.odefun = ode;
        sol.extdata.options = options;                       
        sol.extdata.varargin = varargin;  
    end  
    
    [neq, tspan, ntspan, next, t0, tfinal, tdir, y0, f0, odeArgs, odeFcn, ...
     options, threshold, rtol, normcontrol, normy, hmax, htry, htspan, dataType] = ...
        odearguments(FcnHandlesUsed, solver_name, ode, tspan, y0, options, varargin);
    nfevals = nfevals + 1;
    
    % Set Output
    if nargout > 0
        outputFcn = odeget(options,'OutputFcn',[],'fast');
    else
        outputFcn = odeget(options,'OutputFcn',@odeplot,'fast');
    end
    outputArgs = {};
    if isempty(outputFcn)
        haveOutputFcn = false;
    else
        haveOutputFcn = true;
        outputs = odeget(options,'OutputSel',1:neq,'fast');
        if isa(outputFcn,'function_handle')
            outputArgs = varargin;
        end
    end
    
    refine = max(1,odeget(options,'Refine',4,'fast'));
    if ntspan > 2
        outputAt = 'RequestedPoints';
    elseif refine <= 1
        outputAt = 'SolverSteps';
    else
        outputAt = 'RefinedSteps';
        S = (1:refine-1)/refine;
    end
    printstats = strcmp(odeget(options,'Stats','off','fast'),'on');
    
    [haveEventFcn,eventFcn,eventArgs,valt,teout,yeout,ieout] = ...
        odeevents(FcnHandlesUsed,odeFcn,t0,y0,options,varargin);
    
    % Handle mass matrix
    [Mtype, M, Mfun] =  odemass(FcnHandlesUsed,odeFcn,t0,y0,options,varargin);
    if Mtype > 0
        Msingular = odeget(options,'MassSingular','no','fast');
        if strcmp(Msingular,'maybe')
            warning('Assuming that the matrix is not singular...');
        elseif strcmp(Msingular,'yes')
            error('Mass matrix is singular')
        end
        
        % NEED TO UPDATE!!!!
        [odeFcn,odeArgs] = odemassexplicit(FcnHandlesUsed,Mtype,odeFcn,odeArgs,Mfun,M);
        f0 = feval(odeFcn,t0,y0,odeArgs{:});
        nfevals = nfevals + 1;
    end
    
    % Non-Negative solution components
    idxNonNegative = odeget(options,'NonNegative',[],'fast');
    nonNegative = ~isempty(idxNonNegative);
    
    
    % NEED TO UPDATE!!!
    if nonNegative  % modify the derivative function
        [odeFcn,thresholdNonNegative] = odenonnegative(odeFcn,y0,threshold,idxNonNegative);
        f0 = feval(odeFcn,t0,y0,odeArgs{:});
        nfevals = nfevals + 1;
    end
    
    t = t0;
    y = y0;
    
    % Generate Output if needed
    nout = 0;
    tout = [];
    yout = [];
    if nargout > 0
        if output_sol
            chunk = min(max(100,50*refine), refine+floor((2^11)/neq));      
            tout = zeros(1,chunk,dataType);
            yout = zeros(neq,chunk,dataType);
            f3d  = zeros(neq,10,chunk,dataType);
        else
            if ntspan > 2
                tout = zeros(1,ntspan,dataType);
                yout = zeros(neq,ntspan,dataType);
            else
                chunk = min(max(100,50*refine), refine+floor((2^13)/neq));
                tout = zeros(1,chunk,dataType);
                yout = zeros(neq,chunk,dataType);
            end
        end
        nout = 1;
        tout(nout) = t;
        yout(:,nout) = y;
    end
    
    % Initialize method parameters
    [pow,A,B,E] = getCoefs;
    f = zeros(neq,length(E),dataType);
    hmin = 16*eps(t);
    
    if isempty(htry)
        absh = min(hmax,htspan);
        if normcontrol 
            rh = (norm(f0) / max(normy,threshold)) / (0.8*rtol^pow);
        else
            rh = norm(f0 ./ max(abs(y),threshold),inf) / (0.8*rtol^pow);
        end
        if absh*rh>1
            absh = 1/rh;
        end
        absh = max(absh,hmin);
    else
        absh = min(hmax,max(hmin,htry));
    end
    f(:,1) = f0;
    
    % Initialize the output function
    if haveOutputFcn
        feval(outputFcn,[t,tfinal],y(outputs),'init',outputArgs{:});
    end
    
    % Main Loop
    
    done = false;
    while ~done
        hmin = 16*eps(t);
        absh = min(hmax,max(hmin,absh));
        h = tdir*absh;
        
        if absh > abs(tfinal - t)
            h = tfinal - t;
            absh = abs(h);
            done = true;
        end
        
        nofailed = true;
        while true
            hA = h*A;
            hB = h*B;
            for m = 2:length(E)
                f(:,m) = feval(odeFcn,t+hA(m),y+f*hB(:,m),odeArgs{:});
            end
            
            tnew = t + h;
            ynew = y + f*hB(:,m);
            nfevals = nfevals + m;
            
            % Estimate Error
            NNrejectStep = false;
            if normcontrol
                normynew = norm(ynew);
                errwt = max(max(normy,normynew),threshold);
                err = absh * (norm(f * E) / errwt);
                if nonNegative && (err <= rtol) && any(ynew(idxNonNegative)<0)
                    errNN = norm( max(0,-ynew(idxNonNegative)) ) / errwt ;
                        if errNN > rtol
                            err = errNN;
                            NNrejectStep = true;
                        end
                end      
            else
                err = absh * norm((f * E) ./ max(max(abs(y),abs(ynew)),threshold),inf);
                if nonNegative && (err <= rtol) && any(ynew(idxNonNegative)<0)
                    errNN = norm( max(0,-ynew(idxNonNegative)) ./ thresholdNonNegative, inf);      
                        if errNN > rtol
                            err = errNN;
                            NNrejectStep = true;
                        end
                end            
            end
            
            % Check if error is accaptable
            if err > rtol
                nfailed = nfailed + 1;
                if absh <= hmin
                    warning('Integration tollerance was not met...')
                    % NEED TO UPDATE!!!!!!!!!!!!!!!!!!!
                    solver_output = odefinalize(solver_name, sol,...
                                                printstats,...
                                                [nsteps,nfailed,nfevals],...
                                                nout,tout,yout,haveEventFcn,...
                                                teout,yeout,ieout,...
                                                {f3d,idxNonNegative});
                    if nargout > 0
                        varargout = solver_output;
                    end
                    return
                end
                
                if nofailed
                    nofailed = false;
                    if NNrejectStep
                        absh = max(hmin,0.5*absh);
                    else
                        absh = max(hmin, absh*max(0.1, 0.8 * (rtol/err)^pow));
                    end
                else
                    absh = max(hmin, 0.5*absh);
                end
                h = tdir * absh;
                done = false;
            else
                break;
            end
        end
        nsteps = nsteps + 1;
        
        if haveEventFcn
            [te,ye,ie,valt,stop] = ...
                odezero(@ntrp56,eventFcn,eventArgs,valt,t,y,tnew,ynew,t0,h,f,idxNonNegative);
            if ~isempty(te)
                if output_sol || (nargout > 2)
                    teout = [teout, te];
                    yeout = [yeout, ye];
                    ieout = [ieout, ie];
                end
                if stop               % Stop on a terminal event.               
                % Adjust the interpolation data to [t te(end)].   

                % Update the derivatives using the interpolating polynomial.
                    taux = t + (te(end) - t)*A;        
                    [~,f] = ntrp56(taux,t,y,[],[],h,f,idxNonNegative);        

                    tnew = te(end);
                    ynew = ye(:,end);
                    h = tnew - t;
                    done = true;
                end
            end
        end
        
        if output_sol
            nout = nout + 1;
            if nout > length(tout)
                tout = [tout, zeros(1,chunk,dataType)];  % requires chunk >= refine
                yout = [yout, zeros(neq,chunk,dataType)];
                f3d  = cat(3,f3d,zeros(neq,10,chunk,dataType)); 
            end
            tout(nout) = tnew;
            yout(:,nout) = ynew;
            f3d(:,:,nout) = f;
        end 
        
        if output_ty || haveOutputFcn
            switch outputAt
                case 'SolverSteps'
                    nout_new = 1;
                    tout_new = tnew;
                    yout_new = ynew;
                case 'RefinedSteps'
                    tref = t + (tnew-t)*S;
                    nout_new = refine;
                    tout_new = [tref,tnew];
                    yout_new = [ntrp56(tref,t,y,[],[],h,f,idxNonNegative),ynew];
                case 'RequestedPoints'
                    nout_new = 0;
                    tout_new = [];
                    yout_new = [];
                    while next <= ntspan
                        if tdir * (tnew - tspan(next)) < 0
                            if haveEventFcn && stop
                                nout_new = nout_new + 1;
                                tout_new = [tout_new, tnew];
                                yout_new = [yout_new, yout];
                            end
                            break;
                        end
                        nout_new = nout_new + 1;
                        tout_new = [tout_new, tspan(next)];
                        if tspan(next) == tnew
                            yout_new = [yout_new,ynew];
                        else
                            yout_new = [yout_new, ntrp56(tspan(next),t,y,[],[],h,f,idxNonNegative)];
                        end
                        next = next + 1;
                    end
            end
            
            if nout_new > 0
                if output_ty
                    oldnout = nout;
                    nout = nout + nout_new;
                    if nout > length(tout)
                        tout = [tout, zeros(1,chunk,dataType)];  % requires chunk >= refine
                        yout = [yout, zeros(neq,chunk,dataType)];
                    end
                    idx = oldnout+1:nout;        
                    tout(idx) = tout_new;
                    yout(:,idx) = yout_new;
                end
                if haveOutputFcn
                    stop = feval(outputFcn,tout_new,yout_new(outputs,:),'',outputArgs{:});
                        if stop
                            done = true;
                        end  
                end     
            end  
        end
        
        if done
            break
        end
        
        % If there were no failures compute a new step size (h)
        if nofailed
            temp = 1.25*(err/rtol)^pow;
            if temp > 0.2
                absh = absh / temp;
            else
                absh = 5*absh;
            end
        end
        
        % Advance the integration one step...
        t = tnew;
        y = ynew;
        if normcontrol
            normy = normynew;
        end
        f(:,1) = f(:,end);
    end
    solver_output = odefinalize(solver_name, sol,...
                            outputFcn, outputArgs,...
                            printstats, [nsteps, nfailed, nfevals],...
                            nout, tout, yout,...
                            haveEventFcn, teout, yeout, ieout,...
                            {f3d,idxNonNegative});
    if nargout > 0
        if nargout == 1 && length(tspan)>2
            solver_output{1}.y = deval56(solver_output{1},tspan);
            solver_output{1}.x = tspan.';
            solver_output{1}.stats.usedInterp = 'true';
        elseif nargout == 1 && length(tspan) == 2
            solver_output{1}.stats.usedInterp = 'false';
        end
        varargout = solver_output;
    end  
    
end

function [pow, A, B, E] = getCoefs
    
    pow = 1/6;
    
    B = zeros(10);
    B(2,1)   = 1/32;
    B(3,1:2) = 1/72*[1 2];
    B(4,1:3) = 1/64 * [1 0 3];
    B(5,1:4) = 1/125 * [53 0 -204 176];
    B(6,1:5) = [1/96 0 0 4/33 125/1056];
    B(7,1:6) = [-19/24 0 0 64/33 -875/264 8/3];
    B(8,1:7) = [987/320 0 0 -36/5 975/64 -459/40 177/160];
    B(9,1:8) = [-1238/105 0 0 33536/1155 -14900/231 1796/35 -21/5 8/7];
    B(10,1:9)= 1/90*[7 0 0 0 0 32 12 32 7];
    B = B';
    A = sum(B);
    b = 1/450*[
        35 0 0 0 0 160 60 160 28  7
        35 0 0 0 0 160 60 160 42 -7
        ];
    E = (b(1,:)-b(2,:)).';
end
    
function [neq, tspan, ntspan, next, t0, tfinal, tdir, y0, f0, args, odeFcn, ...
          options, threshold, rtol, normcontrol, normy, hmax, htry, htspan, ...
          dataType ] =   ...
    odearguments(FcnHandlesUsed, solver, ode, tspan, y0, options, extras)

    if FcnHandlesUsed  % function handles used
        if isempty(tspan) || isempty(y0) 
            error(['tspan or y0 was not supplied for ', solver]);
        end      
        if length(tspan) < 2
            error(['Size of tspan must be greater than one for ',solver]);%message('MATLAB:odearguments:SizeTspan', solver));
        end  
        htspan = abs(tspan(2) - tspan(1));  
        tspan = tspan(:);
        ntspan = length(tspan);
        t0 = tspan(1);  
        next = 2;       % next entry in tspan
        tfinal = tspan(end);     
        args = extras;                 % use f(t,y,p1,p2...) 
    else
        if isempty(tspan) || isempty(y0) 
            if exist(ode)==2 && ( nargout(ode)<3 && nargout(ode)~=-1 ) 
                error(message('MATLAB:odearguments:NoDefaultParams', funstring( ode ), solver, funstring( ode )));      
            end
            [def_tspan,def_y0,def_options] = feval(ode,[],[],'init',extras{:});
            if isempty(tspan)
                tspan = def_tspan;
            end
            if isempty(y0)
                y0 = def_y0;
            end
            options = odeset(def_options,options);
        end  
        tspan = tspan(:);
        ntspan = length(tspan);
        if ntspan == 1
            t0 = 0;
            next = 1;
        else
            t0 = tspan(1);
            next = 2;
        end
        htspan = abs(tspan(next)-t0);
        tfinal = tspan(end);
        
        % The input arguments of f determine the args to use to evaluate f
        if (exist(ode)==2)
            if (nargin(ode)==2)
                args = {};
            else
                args = [{''} extras];
            end
        else
            try
                args = [{''} extras];
                feval(ode,tspan(1),y0(:),args{:});
            catch %ME
                args = {};
            end
        end
    end
    
    y0 = y0(:);
    neq = length(y0);
    
    % Test that tspan is internally consistent
    if t0 == tfinal
        error('tspan endpoints are not distinct')
    end
    tdir = sign(tfinal-t0);
    if any(tdir*diff(tspan) <= 0)
        error('tspan is not monotonic')
    end
    
    f0 = feval(ode,t0,y0,args{:});
    [m,n] = size(f0);
    if n>1
        error('ode function must return a column vector')
    elseif m ~= neq
        error('Check Initial Conditions')
    end
    
    % Determine class types
    classT0 = class(t0);
    classY0 = class(y0);
    classF0 = class(f0);
    
    dataType = superiorfloat(t0,y0,f0);
    
    if ~(strcmp(classT0,dataType) && strcmp(classY0,dataType) && strcmp(classF0,dataType))
        error('Inconsistent datatypes...')
    end
    
    rtol = odeget(options,'RelTol',1e-3,'fast');
    if rtol < 100*eps(dataType)
        rtol = 100*eps(dataType);
        warning(['Increasing rtol, due to datatypes, to ', sprintf('%g',rtol)]);
    end
    atol = odeget(options,'AbsTol',1e-6,'fast');
    if any(atol <= 0)
        error('AbsTol not possible...')
    end
    normcontrol = strcmp(odeget(options,'NormControl','off','fast'),'on');
    if normcontrol
        if length(atol) ~= 1
            error('Non-Scalar AbsTol...')
        end
        normy = norm(y0);
    else
        if (length(atol) ~= 1) && (length(atol) ~= neq)
            error('size of AbsTol...');
        end
        atol = atol(:);
        normy = [];
    end
    threshold = atol/rtol;
    
    % Set hmax
    hmax = min(abs(tfinal-t0),abs(odeget(options,'MaxStep',0.1*(tfinal-t0),'fast')));
    if hmax <= 0
        error(' Max step size is zero....');
    end
    htry = abs(odeget(options,'InitialStep',[],'fast'));
    if ~isempty(htry) && (htry <= 0)
        error('Initial step size is zero...')
    end
    
    odeFcn = ode;
end

function [haveeventfun,eventFcn,eventArgs,eventValue,teout,yeout,ieout] =...
    odeevents(FcnHandlesUsed,ode,t0,y0,options,extras) 
    
    haveeventfun = 0;
    eventArgs = [];
    eventValue = [];
    teout = [];
    yeout = [];
    ieout = [];
    
    eventFcn = odeget(options,'Events',[],'fast');
    if isempty(eventFcn)
        return
    end
    
    if FcnHandlesUsed
        haveeventfun = 1;
        eventArgs = extras;
        eventValue = feval(eventFcn,t0,y0,eventArgs{:});
    else
        switch lower(eventFcn)
            case 'on'
                haveeventfun = 1;
                eventFcn = ode;
                eventArgs = [{'events'},extras];
                eventValue = feval(eventFcn,t0,y0,eventArgs{:});
            case 'off'
            otherwise
                error('Must set ode events either on or off...')
        end
    end
end

function [massType, massM, massFcn, massArgs, dMoptions] = ...
    odemass(FcnHandlesUsed,ode,t0,y0,options,extras)

    massType = 0;
    massFcn = [];
    massArgs = {};
    massM = speye(length(y0));
    dMoptions = [];
    
    if FcnHandlesUsed
        Moption = odeget(options,'Mass',[],'fast');
        if isempty(Moption)
            return
        elseif isnumeric(Moption)
            massType = 1;
            massM = Moption;
        else
            massFcn = Moption;
            massArgs = extras;
            Mstdep = odeget(options,'MStateDependence','weak','fast');
            switch lower(Mstdep)
                case 'none'
                    massType = 2;
                case 'weak' 
                    massType = 3;
                case 'strong'
                    massType = 4;
                    
                    dMoptions.diffvar = 3;
                    dMoptions.vectvars = [];
                    
                    atol = odeget(options,'AbsTol',1e-6,'fast');
                    dMoptions.thresh = zeros(size(y0)) + atol(:);
                    
                    dMoptions.fac = [];
                    
                    Mvs = odeget(options,'MvPattern',[],'fast');
                    if ~isempty(Mvs)
                        dMoptions.pattern = Mvs;
                        dMoptions.g = colgroup(Mvs);
                    end
                    
                otherwise
                    error('Check Mass Type Dependency')
            end
            
            if massType > 2
                massM = feval(massFcn,t0,y0,massArgs{:});
            else
                massM = feval(massFcn,t0,massArgs{:});
            end
        end
        
    else
        mass = lower(odeget(options,'Mass','none','fast'));
        
        switch mass
            case 'none', return;
            case 'm',      massType = 1;
            case 'm(t)',   massType = 2;
            case 'm(t,y)', massType = 3;
            otherwise
                error('Check Mass Properties')
        end
        
        massFcn = ode;
        massArgs = [{'mass'},extras];
        massM = feval(massFcn,t0,y0,massArgs{:});
    end
end
    
function solver_output = odefinalize(solver, sol,...
                                     outfun, outargs,...
                                     printstats, statvect,...
                                     nout, tout, yout,...
                                     haveeventfun, teout, yeout, ieout,...
                                     interp_data)
%ODEFINALIZE Helper function called by ODE solvers at the end of integration.
%
%   See also ODE113, ODE15I, ODE15S, ODE23, ODE23S, 
%            ODE23T, ODE23TB, ODE45, DDE23, DDESD.

%   Jacek Kierzenka
%   Copyright 1984-2005 The MathWorks, Inc.

    if ~isempty(outfun)
        feval(outfun,[],[],'done',outargs{:});
    end

% Return more stats for implicit solvers: ODE15i, ODE15s, ODE23s, ODE23t, ODE23tb
    fullstats = (length(statvect) > 3);  % faster than 'switch' or 'ismember'

    stats = struct('nsteps',statvect(1),'nfailed',statvect(2),'nfevals',statvect(3)); 
    if fullstats
        stats.npds     = statvect(4);
        stats.ndecomps = statvect(5);
        stats.nsolves  = statvect(6);  
    else 
        statvect(4:6) = 0;   % Backwards compatibility
    end  

    if printstats
        fprintf(getString(message('MATLAB:odefinalize:LogSuccessfulSteps', sprintf('%g',stats.nsteps))));
        fprintf(getString(message('MATLAB:odefinalize:LogFailedAttempts', sprintf('%g',stats.nfailed))));
        fprintf(getString(message('MATLAB:odefinalize:LogFunctionEvaluations', sprintf('%g',stats.nfevals))));
        if fullstats
            fprintf(getString(message('MATLAB:odefinalize:LogPartialDerivatives', sprintf('%g',stats.npds))));
            fprintf(getString(message('MATLAB:odefinalize:LogLUDecompositions', sprintf('%g',stats.ndecomps))));
            fprintf(getString(message('MATLAB:odefinalize:LogSolutionsOfLinearSystems', sprintf('%g',stats.nsolves))));
        end
    end

    solver_output = {};

    if (nout > 0) % produce output
        if isempty(sol) % output [t,y,...]
            solver_output{1} = tout(1:nout).';
            solver_output{2} = yout(:,1:nout).';
            if haveeventfun
                solver_output{3} = teout.';
                solver_output{4} = yeout.';
                solver_output{5} = ieout.';
            end
            solver_output{end+1} = statvect(:);  % Column vector
        else % output sol  
        % Add remaining fields
            sol.x = tout(1:nout);
            sol.y = yout(:,1:nout);
            if haveeventfun
                sol.xe = teout;
                sol.ye = yeout;
                sol.ie = ieout;
            end
            sol.stats = stats;
            switch solver
                case 'ode56'
                    [f3d,idxNonNegative] = deal(interp_data{:});
                    sol.idata.f3d = f3d(:,:,1:nout);      
                    sol.idata.idxNonNegative = idxNonNegative;
                otherwise
                    error(message('MATLAB:odefinalize:UnrecognizedSolver', solver));
            end  
            if strcmp(solver,'dde23') || strcmp(solver,'ddesd')
                solver_output = sol;
            else  
                solver_output{1} = sol;
            end  
        end
    end
        
end

function [yinterp,ypinterp] = ntrp56(tinterp,t,y,~,~,h,f,idxNonNegative)

    
    BI = zeros(10,5);
%     BI(1,1) = 1;
%     BI(:,2) = 1/6*[-25 0 0 0 0   48 -36   16 -1 -2];
%     BI(:,3) = 1/9*[ 70 0 0 0 0 -208 228 -112  1 21];
%     BI(:,4) = 2/3*[-10 0 0 0 0   36 -48   28  1 -7];
%     BI(:,5) = 8/15*[ 4 0 0 0 0  -16  24  -16 -1  5];
    BI(1,1) = 1;
    BI(:,2) = 0.166666666666666666666667*[-25 0 0 0 0   48 -36   16 -1 -2];
    BI(:,3) = 0.111111111111111111111111*[ 70 0 0 0 0 -208 228 -112  1 21];
    BI(:,4) = 0.666666666666666666666667*[-10 0 0 0 0   36 -48   28  1 -7];
    BI(:,5) = 0.533333333333333333333333*[  4 0 0 0 0  -16  24  -16 -1  5];
    
    s = (tinterp - t)/h;
    yinterp = y(:,ones(size(tinterp))) + f*(h*BI)*cumprod([s;s;s;s;s]);
    
    ypinterp = [];
    if nargout > 1
        ypinterp = f*BI*[ones(size(s));cumprod([2*s;3/2*s;4/3*s;5/4*s])];
    end
    
    if ~isempty(idxNonNegative)
        idx = find(yinterp(idxNonNegative,:)<0);
        if ~isempty(idx)
            w = yinterp(idxNonNegative,:);
            w(idx) = 0;
            yinterp(idxNonNegative,:) = w;
            if nargout > 1
                w = ypinterp(idxNonNegative,:);
                w(idx) = 0;
                ypinterp(idxNonNegative,:) = w;
            end
        end
    end
end

function [tout,yout,iout,vnew,stop] = ...
    odezero(ntrpfun,eventfun,eventargs,v,t,y,tnew,ynew,t0,varargin)
%ODEZERO Locate any zero-crossings of event functions in a time step.
%   ODEZERO is an event location helper function for the ODE Suite.  ODEZERO
%   uses Regula Falsi and information passed from the ODE solver to locate
%   any zeros in the half open time interval (T,TNEW] of the event functions
%   coded in eventfun.
%   
%   See also ODE45, ODE23, ODE113, ODE15S, ODE23S, ODE23T, ODE23TB.

%   Mark W. Reichelt, Lawrence F. Shampine, 6-13-94
%   Copyright 1984-2004 The MathWorks, Inc.

% Initialize.
tol = 128*max(eps(t),eps(tnew));
tol = min(tol, abs(tnew - t));
tout = [];
yout = [];
iout = [];
tdir = sign(tnew - t);
stop = 0;
rmin = realmin;

% Set up tL, tR, yL, yR, vL, vR, isterminal and direction.
tL = t;
yL = y;
vL = v;
[vnew,isterminal,direction] = feval(eventfun,tnew,ynew,eventargs{:});
if isempty(direction)
  direction = zeros(size(vnew));   % zeros crossings in any direction
end
tR = tnew;
yR = ynew;
vR = vnew;

% Initialize ttry so that we won't extrapolate if vL or vR is zero.
ttry = tR;

% Find all events before tnew or the first terminal event.
while true
  
  lastmoved = 0;
  while true
    % Events of interest shouldn't have disappeared, but new ones might
    % be found in other elements of the v vector.
    indzc = find((sign(vL) ~= sign(vR)) & (direction .* (vR - vL) >= 0));
    if isempty(indzc)
      if lastmoved ~= 0
        error(message('MATLAB:odezero:LostEvent'));
      end
      return;
    end
    
    % Check if the time interval is too short to continue looking.
    delta = tR - tL;
    if abs(delta) <= tol
      break;
    end
    
    if (tL == t) && any(vL(indzc) == 0 & vR(indzc) ~= 0)
      ttry = tL + tdir*0.5*tol;
      
    else
      % Compute Regula Falsi change, using leftmost possibility.
      change = 1;
      for j = indzc(:)'
        % If vL or vR is zero, try using old ttry to extrapolate.
        if vL(j) == 0
          if (tdir*ttry > tdir*tR) && (vtry(j) ~= vR(j))
            maybe = 1.0 - vR(j) * (ttry-tR) / ((vtry(j)-vR(j)) * delta);
            if (maybe < 0) || (maybe > 1)
              maybe = 0.5;
            end
          else
            maybe = 0.5;
          end
        elseif vR(j) == 0.0
          if (tdir*ttry < tdir*tL) && (vtry(j) ~= vL(j))
            maybe = vL(j) * (tL-ttry) / ((vtry(j)-vL(j)) * delta);
            if (maybe < 0) || (maybe > 1)
              maybe = 0.5;
            end
          else
            maybe = 0.5;
          end
        else
          maybe = -vL(j) / (vR(j) - vL(j)); % Note vR(j) ~= vL(j).
        end
        if maybe < change
          change = maybe;
        end
      end
      change = change * abs(delta);

      % Enforce minimum and maximum change.
      change = max(0.5*tol, min(change, abs(delta) - 0.5*tol));

      ttry = tL + tdir * change;
    end
    
    % Compute vtry.
    ytry = feval(ntrpfun,ttry,t,y,tnew,ynew,varargin{:});
    vtry = feval(eventfun,ttry,ytry,eventargs{:});

    % Check for any crossings between tL and ttry.
    indzc = find((sign(vL) ~= sign(vtry)) & (direction .* (vtry - vL) >= 0));
    if ~isempty(indzc)
      % Move right end of bracket leftward, remembering the old value.
      tswap = tR; tR = ttry; ttry = tswap;
      yswap = yR; yR = ytry; ytry = yswap;
      vswap = vR; vR = vtry; vtry = vswap;
      % Illinois method.  If we've moved leftward twice, halve
      % vL so we'll move closer next time.
      if lastmoved == 2
        % Watch out for underflow and signs disappearing.
        maybe = 0.5 * vL;
        i = find(abs(maybe) >= rmin);
        vL(i) = maybe(i);
      end
      lastmoved = 2;
    else
      % Move left end of bracket rightward, remembering the old value.
      tswap = tL; tL = ttry; ttry = tswap;
      yswap = yL; yL = ytry; ytry = yswap;
      vswap = vL; vL = vtry; vtry = vswap;
      % Illinois method.  If we've moved rightward twice, halve
      % vR so we'll move closer next time.
      if lastmoved == 1
        % Watch out for underflow and signs disappearing.
        maybe = 0.5 * vR;
        i = find(abs(maybe) >= rmin);
        vR(i) = maybe(i);
      end
      lastmoved = 1;
    end
  end

  j = ones(1,length(indzc));
  tout = [tout, tR(j)];
  yout = [yout, yR(:,j)];
  iout = [iout, indzc'];
  if any(isterminal(indzc))
    if tL ~= t0
      stop = 1;
    end
    break;
  elseif abs(tnew - tR) <= tol
    %  We're not going to find events closer than tol.
    break;
  else
    % Shift bracket rightward from [tL tR] to [tR+0.5*tol tnew].
    ttry = tR; ytry = yR; vtry = vR;
    tL = tR + tdir*0.5*tol;
    yL = feval(ntrpfun,tL,t,y,tnew,ynew,varargin{:});
    vL = feval(eventfun,tL,yL,eventargs{:});
    tR = tnew; yR = ynew; vR = vnew;
  end
end
end

function [odeFcn,odeArgs] = odemassexplicit( FcnHandlesUsed,massType,odeFcn,...
                                             odeArgs,massFcn,massM)  
%ODEMASSEXPLICIT  Helper function for handling the mass matrix
%   For explicit ODE solvers -- incorporate the mass matrix into the ODE
%   function.   
%
%   See also ODE113, ODE23, ODE45.

%   Jacek Kierzenka
%   Copyright 1984-2006 The MathWorks, Inc.

if FcnHandlesUsed
    switch massType
      case 1  % use LU factors of constant M 
        if issparse(massM) 
            [massL,massU,massP,massQ,massR] = lu(massM);
            odeArgs = [{odeFcn,massL,massU,massP,massQ,massR},odeArgs];    
            odeFcn = @ExplicitSolverHandleMass1sparse;       
        else % M full
            [massL,massU,massp] = lu(massM,'vector');
            odeArgs = [{odeFcn,massL,massU,massp},odeArgs];    
            odeFcn = @ExplicitSolverHandleMass1;
        end  
      case 2
        odeArgs = [{odeFcn,massFcn},odeArgs];    
        odeFcn = @ExplicitSolverHandleMass2;
      otherwise % case {3,4}
        odeArgs = [{odeFcn,massFcn},odeArgs];    
        odeFcn = @ExplicitSolverHandleMass34;
    end
else % ode-file:  F(t,y,'mass',p1,p2...)    
    if massType == 1   % use LU factors of constant M 
        if issparse(massM) 
            [massL,massU,massP,massQ,massR] = lu(massM);
            odeArgs = [{odeFcn,massL,massU,massP,massQ,massR},odeArgs];    
            odeFcn = @ExplicitSolverHandleMass1sparse;       
        else % M full
            [massL,massU,massp] = lu(massM,'vector');
            odeArgs = [{odeFcn,massL,massU,massp},odeArgs];    
            odeFcn = @ExplicitSolverHandleMass1;
        end  
    else  
        odeArgs = [{odeFcn},odeArgs];  
        odeFcn = @ExplicitSolverHandleMassOld;   
    end
end
end
% --------------------------------------------------------------------------

function yp = ExplicitSolverHandleMass1(t,y,odeFcn,L,U,p,varargin)
  ode = feval(odeFcn,t,y,varargin{:});
  yp = U \ (L \ ode(p));

% --------------------------------------------------------------------------
end
function yp = ExplicitSolverHandleMass1sparse(t,y,odeFcn,L,U,P,Q,R,varargin)
  yp = Q *( U \ (L \ (P * (R \ feval(odeFcn,t,y,varargin{:})))));
 
% --------------------------------------------------------------------------
end  
function yp = ExplicitSolverHandleMass2(t,y,odeFcn,massFcn,varargin)
  yp = feval(massFcn,t,varargin{:}) \ feval(odeFcn,t,y,varargin{:});
  
% --------------------------------------------------------------------------  
end
function yp = ExplicitSolverHandleMass34(t,y,odeFcn,massFcn,varargin)
  yp = feval(massFcn,t,y,varargin{:}) \ feval(odeFcn,t,y,varargin{:});

% --------------------------------------------------------------------------  
end  
function yp = ExplicitSolverHandleMassOld(t,y,odeFcn,varargin)
  yp = feval(odeFcn,t,y,'mass',varargin{2:end}) \ ...
       feval(odeFcn,t,y,varargin{:});
  
% --------------------------------------------------------------------------
end

function [Sxint,Spxint] = deval56(sol,xint,idx)
    
    if ~isa(sol,'struct')
        temp = sol;
        sol = xint;
        xint = temp;
    end
    
    try
        t = sol.x;
        y = sol.y;
    catch
        error('Solution not from diffeq solver')
    end
    
    if nargin < 3
        idx = 1:size(y,1);
    else
        if any(idx < 0) || any(idx > size(y,1))
            error(message('MATLAB:deval:IDXInvalidSolComp', inputname( 3 )));
        end  
    end  
    idx = idx(:);
        
    if isfield(sol,'solver')
        solver = sol.solver;
    else
        if isfield(sol,'yp')
            warning(message('MATLAB:deval:MissingSolverField', inputname(1), inputname(1)));  
            solver = 'bvp4c';
        else
            error(message('MATLAB:deval:NoSolverInStruct',inputname(1)));
        end
    end
    
    interpfcn = @ntrp56;
    dataType = superiorfloat(sol.x,xint);
    Spxint_requested = (nargout > 1);   
    
    if isfield(sol, 'idata') && isfield(sol.idata, 'idxNonNegative')
        idxNonNegative = sol.idata.idxNonNegative;
    else
        idxNonNegative = [];
    end
    
    n = length(idx);
    Nxint = length(xint);
    Sxint = zeros(n,Nxint,dataType);
    if Spxint_requested
        Spxint = zeros(n,Nxint,dataType);
    end

    % Make tint a row vector and if necessary, 
    % sort it to match the order of t.
    tint = xint(:).';  
    tdir = sign(t(end) - t(1));
    had2sort = any(tdir*diff(tint) < 0);
    if had2sort
        [tint,tint_order] = sort(tdir*tint);
        tint = tdir*tint;
    end  

    % Using the sorted version of tint, test for illegal values.
    if any(isnan(tint)) || (tdir*(tint(1) - t(1)) < 0) || (tdir*(tint(end) - t(end)) > 0)
        error(message('MATLAB:deval:SolOutsideInterval',sprintf('%e',t(1)),sprintf('%e',t(end))));
    end

    evaluated = 0;
    bottom = 1;
    while evaluated < Nxint
  
        % Find right-open subinterval [t(bottom), t(bottom+1)) containing the next entry of tint. 
        % Unless t(bottom) == t(end), a proper interval is returned: t(bottom+1) ~= t(bottom).
        Index = find( tdir*(t(bottom:end) - tint(evaluated+1)) <= 0, 1, 'last'); % one-based Index
        bottom = bottom + (Index - 1);  % convert to zero-based Index

      % Is it [t(end), t(end)]?
      at_tend = (t(bottom) == t(end));

      % Return solution already available at t(bottom)
      index1 = find(tint(evaluated+1:end) == t(bottom));

      % Interpolate solution inside (t(bottom), t(bottom+1))
      if at_tend
        index2 = [];
      else
        index2 = find( (tdir*(tint(evaluated+1:end) - t(bottom)) > 0) & ...
                       (tdir*(tint(evaluated+1:end) - t(bottom+1)) < 0) );
      end

      % Return the (adjusted) solution at t(bottom)
      if ~isempty(index1)
          if at_tend          
              yint1 = y(:,end);
              if Spxint_requested
                  % Extrapolate derivative from [t(bottom-1),t(bottom))
                  interpdata = extract_idata(solver,sol,t,bottom-1,idxNonNegative);
                  [~,ypint1] = interpfcn(t(bottom),t(bottom-1),y(:,bottom-1),...
                                       t(bottom),y(:,bottom),interpdata{:});    
              end

          elseif (bottom > 2) && (t(bottom) == t(bottom-1)) % Interface point
              % Average the solution (and its derivative) across the interface.
              yLeft  = y(:,bottom-1);       
              yRight = y(:,bottom);
              yint1 = (yLeft + yRight)/2;          
              if Spxint_requested
                  % Get the 'left' derivative by extrapolating in [t(bottom-2), t(bottom-1)).        
                  interpdata =  extract_idata(solver,sol,t,bottom-2,idxNonNegative); 
                  [~,ypLeft] = interpfcn(t(bottom-1),t(bottom-2),y(:,bottom-2),...
                                        t(bottom-1),y(:,bottom-1),interpdata{:});
                  % Get the 'right' derivative by interpolating in [t(bottom), t(bottom+1)).        
                  interpdata =  extract_idata(solver,sol,t,bottom,idxNonNegative); 
                  [~,ypRight] = interpfcn(t(bottom),t(bottom),y(:,bottom),...
                                        t(bottom+1),y(:,bottom+1),interpdata{:});
                  ypint1 = (ypLeft + ypRight)/2;
              end
              warning(message('MATLAB:deval:NonuniqueSolution',sprintf('%g',t(bottom))));
          else  
              % 'Regular' mesh point
              yint1 = y(:,bottom); 
              if Spxint_requested
                  % Interpolate derivative from [t(bottom),t(bottom+1))
                  interpdata = extract_idata(solver,sol,t,bottom,idxNonNegative);
                  [~,ypint1] = interpfcn(t(bottom),t(bottom),y(:,bottom),...
                                       t(bottom+1),y(:,bottom+1),interpdata{:});    
              end                    
          end

          % Accumulate the output.
          Sxint(:,evaluated+index1) = yint1(idx,ones(1,numel(index1)));  
          if Spxint_requested
              Spxint(:,evaluated+index1) = ypint1(idx,ones(1,numel(index1)));  
          end  
      end

      % Interpolate solution inside (t(bottom), t(bottom+1)).
      if ~isempty(index2) 
          % Get solver-dependent interpolation data for [t(bottom), t(bottom+1)).
          interpdata = extract_idata(solver,sol,t,bottom,idxNonNegative);

          % Evaluate the interpolant at all points from (t(bottom), t(bottom+1)).
          if Spxint_requested
              [yint2,ypint2] = interpfcn(tint(evaluated+index2),t(bottom),y(:,bottom),...
                                       t(bottom+1),y(:,bottom+1),interpdata{:});    
          else  
              yint2 = interpfcn(tint(evaluated+index2),t(bottom),y(:,bottom),...
                               t(bottom+1),y(:,bottom+1),interpdata{:});    
          end

          % Accumulate the output.
          Sxint(:,evaluated+index2) = yint2(idx,:);  
          if Spxint_requested
              Spxint(:,evaluated+index2) = ypint2(idx,:);  
          end  
      end

      evaluated = evaluated + length(index1) + length(index2);      
    end

    if had2sort     % Restore the order of tint in the output.
      Sxint(:,tint_order) = Sxint;  
      if Spxint_requested
        Spxint(:,tint_order) = Spxint;  
      end  
    end
end

function interpdata = extract_idata(~,sol,t,tidx,idxNonNegative)
% Data for interpolation in [t(tidx), t(tidx+1))

  interpdata = { t(tidx+1)-t(tidx), ...
                 sol.idata.f3d(:,:,tidx+1), ...
                 idxNonNegative };   
end    
    