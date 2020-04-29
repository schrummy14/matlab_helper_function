% Sun Earth Moon System
clear
clc

SunMass = 1*1.9891e30;
MercuryMass = 1*3.3011e23;
VenusMass = 1*4.8675e24;
EarthMass = 1*5.972e24;
MoonMass = 1*7.34767309e22;
MarsMass = 1*6.4171e23;
JupiterMass = 1*1.8986e27;
SaturnMass = 1*5.6836e26;
UranusMass = 1*8.6810e25;
NeptuneMass= 1*1.0243e26; 

OrbitMass = [SunMass,MercuryMass,VenusMass,EarthMass,MoonMass,MarsMass,...
    JupiterMass,SaturnMass,UranusMass,NeptuneMass];

G = 6.67408e-11;

RSMe= 57909050e3; 
RSJ = 778.299e9;
RSV = (108939000000+107477000000)/2;
RSE = 92.96e6*1609.34;
REM = 238900*1609.34;
RSMa= 227.9392e9;
RSS = 1429.39e9;
RSU = 2875.04e9;
RSN = 4504.45e9;

MercuryVelocity = 47.362e3;
VenusVelocity = 35.02e3; % meters WRT Sun
JupiterVelocity = 13070; % WRT Sun
EarthVelocity = 30000; % WRT Sun
MoonVelocity  = 1023; % WRT Earth
MarsVelocity = 24.077e3; % WRT Sun
SaturnVelocity = 9.69e3; % WRT Sun
UranusVelocity = 6.80e3; % WRT Sun
NeptuneVelocity= 5.43e3;

tSpan = linspace(0,1000*365*24*3600,2);

% Set on sun
x0 = [
    0;0 % Sun Pos
    0;RSMe % Mercury
    0;RSV % Venus
    0;RSE % Earth Pos
    0;RSE+REM % Moon Pos
    0;RSMa % Mars
    0;RSJ % Jupiter Pos
    0;RSS
    0;RSU
    0;RSN
    0;0 % Sun Vel
    MercuryVelocity;0
    VenusVelocity;0
    EarthVelocity;0
    -MoonVelocity+EarthVelocity;0
    MarsVelocity;0
    JupiterVelocity;0
    SaturnVelocity;0
    UranusVelocity;0
    NeptuneVelocity;0
    ];

opts = odeset('RelTol',1e-6,'AbsTol',1e-10,'InitialStep',24*3600);
exopts = {'ProgressBar',true};
tic
Sol = ode56(@(t,x)n_body(t,x,OrbitMass,G,2),tSpan,x0.',opts,exopts{:});
toc

x = Sol.y';
t = Sol.x';

%% Plot it
rf = 5;
plot(...
    x(:, 1)-x(:,rf),x(:, 2)-x(:,rf+1),'-',...
    x(:, 3)-x(:,rf),x(:, 4)-x(:,rf+1),'-',...
    x(:, 5)-x(:,rf),x(:, 6)-x(:,rf+1),'-',...
    x(:, 7)-x(:,rf),x(:, 8)-x(:,rf+1),'-',...
    x(:, 9)-x(:,rf),x(:,10)-x(:,rf+1),'-',...
    x(:,11)-x(:,rf),x(:,12)-x(:,rf+1),'-',...
    x(:,13)-x(:,rf),x(:,14)-x(:,rf+1),'-',...
    x(:,15)-x(:,rf),x(:,16)-x(:,rf+1),'-',...
    x(:,17)-x(:,rf),x(:,18)-x(:,rf+1),'-',...
    x(:,19)-x(:,rf),x(:,20)-x(:,rf+1),'-')

axis('equal')
legend('Sun','Mercury','Venus','Earth','Moon','Mars')