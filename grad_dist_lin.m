clear
clc

% rng(1)
% x_data = linspace(-1,1,15)';
% x_data = [x_data, x_data.^2, x_data.^3];
dim = 5;
data_length = 15000;
m0 = rand(1,dim);
x_data = 2*rand(data_length,1) - 1;

% x_data = zeros(data_length,dim);
xx = zeros(data_length,dim);
% x_data(:,1) = 2*rand(data_length,1) - 1;
xx(:,1) = linspace(-1,1,data_length)';
for k = 2:dim
    x_data(:,k) = x_data(:,1).^k;
    xx(:,k) = xx(:,1).^k;
end
y_data = x_data * m0' + 1 + 0.062500*(2*rand(size(x_data,1),1)-1);

[obs,deg] = size(x_data);

tol = 1e-10;
done = false;
m = zeros(1,deg);
b = 0;
step = 1e-2;
b_old = b;
m_old = m;

while ~done
    
    for k = 1:obs
        
        guess = x_data(k,:)*m' + b;
        error = y_data(k,:) - guess;
        b = b + step*error;
        m = m + step*(error * x_data(k,:));

    end
    
    res1 = abs(b - b_old)/b;
    res2 = abs(m - m_old)/norm(m);
    b_old = b;
    m_old = m;
    res = max([res1;res2(:)]);
    
    if res <= tol
        break
    end

end

y_hat = x_data*m' + b;
SS_tot = sum((y_data - mean(y_data)).^2);
SS_res = sum((y_data - y_hat).^2);

r_squared = 1 - SS_res/SS_tot;
adj_r_squared = 1 - ((1-r_squared)*(obs-1))/(obs-deg-1);
s = 1000;

yy = xx*m'+b;
plot(x_data,y_data,'k.',xx,(yy),'r--')
title(['RMS = ',num2str(rms(y_data-(y_hat))),...
       ', R^2_a_d_g = ',num2str(round(s*adj_r_squared)/s)])