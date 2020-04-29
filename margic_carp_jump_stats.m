clc
clear

Data = [
    1 0 2.51
    85980062 0 130.99
    11078387815 0 203.98
    25979076005 0 220.54
    31229190163 0 224.46
    34411545215 0 226.89
    ];

JP = Data(:,1).*(1+Data(:,2));
Height = Data(:,3);

% fun = @(b,x)b(1)+b(2).*x+b(3).*x.^2;
% b = [3 2 0.2];

% plot(log(JP),Height,'kx',log(JP),fun(b,log(JP)))

% b0 = fminunc(@(b)norm(Height-fun(b,JP)),[max(JP),-1]);
LM = fitlm(log(JP),log(Height),'poly3');
LM.disp;

LM.plot

% plot(JP,Height,'kx',JP,LM.Fitted,'r--')