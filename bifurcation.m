%% Program to produce the bifurcation of a single parameter of a function
% Written by Matt Schramm 

% Example bifurcation(linspace(1.4,2.1,100))
function bifurcation(param)

    % if param is only a single value, assume the parameter goes from eps
    % to param using 100 values.
    if length(param)<2
        param = linspace(eps,param,100);
    % If param only contains two values, assume the parameter goes form the
    % first value in param to the last value of param in 100 steps
    elseif length(param)<3
        param = linspace(param(1),param(2),100);
    end
    
    figure
    hold on
    x = zeros(51,length(param));
    parfor m = 1:length(param)
        temp = fun(param(m));
        x(:,m) = temp(end-50:end);
    end
    plot(param,x,'r.')
    xlim([min(param),max(param)]);
    ylim([0 1.5]);
    hold off
    
end


function x = fun(param)

    % Your function will go here. As an example I will investigate the
    % bifurcation due to the step size used by "The Runge Kutta" method in
    % solving the Logistic Map equation dxdt = r*x*(1-x) with r = 2
    
    dxdt = @(x)(2*x*(1-x));
    
    A = [ % A parameter matrix for "The Runge Kutta Method"
        0   0   0   0
        1/2 0   0   0
        0   1/2 0   0
        0   0   1   0
        ];
    b = [1/6 1/3 1/3 1/6];

    k = param;
    N = 200;
    x = zeros(N,1);
    x(1,:) = .99; % You want to start close to the stationary point 
                       % but not directly on it. 
    y = zeros(length(b),1);
    for n = 1:N
        for m = 1:length(b)
            y(m,:) = dxdt(x(n,:)+k*A(m,:)*y);
        end
        x(n+1,:) = x(n,:) + k*b*y;
    end
    
end