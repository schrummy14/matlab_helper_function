time = [0 7  14  21 35 115]./240;
    
poss = ([2   1.66  1.33  1 0.5 0]/2);

xx = linspace(0,115/240,1001)';

[yy,rsq,beta] = LinReg(time,poss,2,xx);

plot(time,poss,'o-')
xlabel 'Time'
ylabel 'Normalized Height'

