function [x,y,z,combos] = BifurcationRe(sig,rho,beta,endTime,Ntrials)

    numPoints = 100;

    %% Set up trials
    % Conditions to be used
    switch numel(sig)
        case 1
            sigmaCon = sig;
        case 2 
            sigmaCon = linspace(sig(1),sig(end),Ntrials);
        otherwise
            sigmaCon = sig;
    end

    switch numel(rho)
        case 1
            rhoCon = rho;
        case 2 
            rhoCon = linspace(rho(1),rho(end),Ntrials);
        otherwise
            rhoCon = rho;
    end

    switch numel(beta)
        case 1
            betaCon = beta;
        case 2 
            betaCon = linspace(beta(1),beta(end),Ntrials);
        otherwise
            betaCon = beta;
    end

    % Set all possible comonations
    sets = {sigmaCon,rhoCon,betaCon};
    [s,r,b] = ndgrid(sets{:});
    combos = [s(:),r(:),b(:)];

    %% Integrate out to steady state starting a small pertibation
    % Only will look at the positive conditions
    
    if isempty(gcp('nocreate'))
        usePar = 0;
    else
        usePar = 1;
    end    
    
    X = cell(1,size(combos,1));
    options = optimoptions('fsolve','Display','off');
    sig = combos(:,1);
    rho = combos(:,2);
    bet = combos(:,3);
    if usePar
        parfor m = 1:size(combos,1)
            %% Steady State Condition
%             sig = combos(m,1);
%             rho = combos(m,2);
%             bet = combos(m,3);
            zs = rho(m) - 1;
            ys = bet(m)*(rho(m)-1);
            xs = ys;
            xstart = fsolve(@(x)Lorenz(0,x,rho(m),sig(m),bet(m)),[zs,ys,xs],options);
            
            [~,X{m}] = ode15s(@(t,x)Lorenz(t,x,rho(m),sig(m),bet(m)),[0 endTime],xstart);
        end
    else
        for m = 1:size(combos,1)
            %% Steady State Condition
%             sig = combos(m,1);
%             rho = combos(m,2);
%             bet = combos(m,3);
            zs = rho(m) - 1;
            ys = bet(m)*(rho(m)-1);
            xs = ys;
            xstart = fsolve(@(x)Lorenz(0,x,rho(m),sig(m),bet(m)),[zs,ys,xs],options);
            [~,X{m}] = ode15s(@(t,x)Lorenz(t,x,rho(m),sig(m),bet(m)),[0 endTime],xstart/1.1);
        end
    end
    x = cell(1,size(combos,1));
    y = x;
    z = x;
    for m = 1:size(combos,1)
        x{m} = X{m}(:,1);
        y{m} = X{m}(:,2);
        z{m} = X{m}(:,3);
    end
    %% Plot the figures
%     xs = zeros(size(combos,1),1);
%     ys = xs;
%     zs = xs;
%     for m = 1:size(combos,1)
%         rho = combos(m,2);
%         bet = combos(m,3);
%         zs(m) = rho - 1;
%         ys(m) = bet*(rho-1);
%         xs(m) = ys(m);
%     end
%     xs = 0;
%     ys = 0;
%     zs = 0;
%     
    % Note: This is only for the Rho Condition
    % Rho and x
    figure
    hold on
    for m = 1:size(combos,1)
        if size(x{m},1) < numPoints
            plot(combos(m,2)*ones(size(x{m},1)),x{m}(:,1),... -xs(m)*ones(size(x{m},1),1),...
                'LineStyle','none','Marker','.','Color','blue');
        else
            plot(combos(m,2)*ones(numPoints,1),x{m}((end-numPoints+1):end,1),...-xs(m)*ones(numPoints,1),...
                'LineStyle','none','Marker','.','Color','blue');
        end
    end
    hold off
    xlabel('Rho Value')
    ylabel('x')
    title('Bifurcation of x with respact to Rho')
    
%     % Rho and y
%     figure
%     hold on
%     for m = 1:size(combos,1)
%         if size(x{m},1) < numPoints
%             plot(combos(m,2)*ones(size(y{m},1)),y{m}(:,1),...-ys(m)*ones(size(y{m},1),1),...
%                 'LineStyle','none','Marker','.','Color','blue');
%         else
%             plot(combos(m,2)*ones(numPoints,1),y{m}((end-numPoints+1):end,1),...-ys(m)*ones(numPoints,1),...
%                 'LineStyle','none','Marker','.','Color','blue');
%         end
%     end
%     hold off
%     xlabel('Rho Value')
%     ylabel('y')
%     title('Bifurcation of y with respact to Rho')
%     
%     % Rho and z
%    figure
%     hold on
%     for m = 1:size(combos,1)
%         if size(x{m},1) < numPoints
%             plot(combos(m,2)*ones(size(z{m},1)),z{m}(:,1),...-zs(m)*ones(size(z{m},1),1),...
%                 'LineStyle','none','Marker','.','Color','blue');
%         else
%             plot(combos(m,2)*ones(numPoints,1),z{m}((end-numPoints+1):end,1),...-zs(m)*ones(numPoints,1),...
%                 'LineStyle','none','Marker','.','Color','blue');
%         end
%     end
%     hold off
%     xlabel('Rho Value')
%     ylabel('z')
%     title('Bifurcation of z with respact to Rho')
    
    
end

function dF = Lorenz(~,x,rho,sigma,beta)
    z = x(3);
    y = x(2);
    x = x(1);
    
    dx = sigma*(y-x);
    dy = x*(rho-z)-y;
    dz = x*y-beta*z;
    
    dF = [
        dx
        dy
        dz
        ];
end