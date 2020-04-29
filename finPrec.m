clc
clear
%% Finite Percision 
% https://en.wikipedia.org/wiki/Finite_difference_coefficient
x0 = .125;
derNumber = 1;
pow = linspace(-16/derNumber,0,80);
h(:,1) = 1*10.^(pow);
yx = @(x)exp(x);
dyT = exp(x0);
% Central Methods
err = zeros(length(pow),4);
for order = 2:2:8
    for n = 1:length(pow)
        x = [x0-4*h(n,1), x0-3*h(n,1), x0-2*h(n,1), x0-1*h(n,1), x0, x0+1*h(n,1), x0+2*h(n,1), x0+3*h(n,1), x0+4*h(n,1)];
        y = yx(x);
        dyN = numDer(y,h(n,1),derNumber,order,'central');
        err(n,order/2) = abs(dyT-dyN);
    end
end
figure
% loglog(h,err)
loglog(h,err(:,1),'k-o',h,err(:,2),'k-s',h,err(:,3),'k-d',h,err(:,4),'k-')
axis('tight')
legend('Order 2','Order 4','Order 6','Order 8','Location','best')
xlabel('Step Size (h)')
ylabel('Abs Error (abs(dy_t_r_u_e - dy_c_a_l_c))')
title('Central Numerical Differentiation Method Errors')
ax = gca;
xlim([h(1) h(end)])
ylim([1e-16 1e0])
ax.XTick = 10.^(-1/derNumber*[16 14 12 10 8 6 4 2 0]);
ax.YTick = [1e-16 1e-14 1e-12 1e-10 1e-8 1e-6 1e-4 1e-2 1e0];

% % Forward Methods
% err = zeros(length(pow),6);
% for order = 1:6
%     for n = 1:length(pow)
%         x = [x0 x0+h(n,1) x0+2*h(n,1) x0+3*h(n,1) x0+4*h(n,1) x0+5*h(n,1) x0+6*h(n,1) x0+7*h(n,1)];
%         y = yx(x);
%         dyN = numDer(y,h(n,1),derNumber,order,'forward');
%         err(n,order) = abs(dyT-dyN);
%     end
% end
% figure
% loglog(h,err)
% axis('tight')
% legend('Order 1','Order 2','Order 3','Order 4','Order 5','Order 6','Location','best')
% xlabel('Step Size (h)')
% ylabel('Abs Error (abs(dy_t_r_u_e - dy_c_a_l_c))')
% title('Forward Numerical Differentiation Method Errors')
% ax = gca;
% xlim([h(1) h(end)])
% ylim([1e-16 1e0])
% ax.XTick = 10.^(-1/derNumber*[16 14 12 10 8 6 4 2 0]);
% ax.YTick = [1e-16 1e-14 1e-12 1e-10 1e-8 1e-6 1e-4 1e-2 1e0];
% 
% % Backward Methods
% err = zeros(length(pow),6);
% for order = 1:6
%     for n = 1:length(pow)
%         x = [x0-7*h(n,1) x0-6*h(n,1),x0-5*h(n,1),x0-4*h(n,1),x0-3*h(n,1),x0-2*h(n,1),x0-h(n,1),x0];
%         y = yx(x);
%         dyN = numDer(y,h(n,1),derNumber,order,'backward');
%         err(n,order) = abs(dyT-dyN);
%     end
% end
% figure
% loglog(h,err)
% axis('tight')
% legend('Order 1','Order 2','Order 3','Order 4','Order 5','Order 6','Location','best')
% xlabel('Step Size (h)')
% ylabel('Abs Error (abs(dy_t_r_u_e - dy_c_a_l_c))')
% title('Backward Numerical Differentiation Method Errors')
% ax = gca;
% xlim([1e-16 1e0])
% ylim([1e-16 1e0])
% ax.XTick = [1e-16 1e-14 1e-12 1e-10 1e-8 1e-6 1e-4 1e-2 1e0];
% ax.YTick = ax.XTick;