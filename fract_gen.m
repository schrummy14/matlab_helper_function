clear
clc

Max_its = 1e6;
N = 3;
phi = linspace(0,2*pi,N+1)+pi/2;
e = exp(phi*1i);

Atrac(:,1) = real(e(1:N));
Atrac(:,2) = imag(e(1:N));

axis_lim = [min(Atrac(:,1)),max(Atrac(:,1)),min(Atrac(:,2)),max(Atrac(:,2))];

shift_dist = 0.5;

points = zeros(Max_its+1,2);
points(1,:) = Atrac(1,:);%0*rand(1,2);

for m = 1:Max_its
    rand_point = randi(N);
    points(m+1,:) = shift_dist*points(m,:) + (1-shift_dist)*Atrac(rand_point,:);
end

% hold on
% for m = 1:Max_its 
%     plot(points(m,1),points(m,2),'.')
%     axis(axis_lim)
% %     axis('equal')
%     drawnow
% end
% hold off

plot(points(:,1),points(:,2),'.')
axis(axis_lim)
axis('equal')