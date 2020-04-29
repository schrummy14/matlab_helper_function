clear
clc

Max_its = 1000000;
N = 3;
phi = linspace(0,2*pi,N+1)+pi/2;
e = exp(phi*1i);

Atrac(:,1) = real(e(1:N));
Atrac(:,2) = imag(e(1:N)) - min(imag(e(1:N)));

axis_lim = [min(Atrac(:,1)),max(Atrac(:,1)),min(Atrac(:,2)),max(Atrac(:,2))];

shift_dist = 0.5;
one_m_shift = 1.0 - shift_dist;

points = zeros(Max_its+1,2);
points(1,:) = Atrac(1,:);

rand_points = randi(N,Max_its,1);

for m = 1:Max_its
    rand_point = rand_points(m);
    points(m+1,:) = shift_dist*points(m,:) + one_m_shift*Atrac(rand_point,:);
end

clear('rand_points')

plot(points(:,1),points(:,2),'.')
axis(axis_lim)
axis('equal')