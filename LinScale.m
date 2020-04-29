function out = LinScale(a,b,x,y,InFile)
    % a is the min value of the scaling you want
    % b is the max value of the scaling you want
    % x is the min value of the scaling you have
    % y is the max value of the scaling you have
    out = (b-a)/(y-x)*(InFile-x) + a;
    
end